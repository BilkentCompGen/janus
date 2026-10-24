#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>
#include <vector>
#include "cmd.h"
#include "io.h"

// `janus run`: the whole pipeline on one normal and one or more tumor samples of a patient, every step a child process
// (janus itself or minimap2) whose wall clock, CPU time and peak memory go to steps.tsv. A step whose marker done.<step>
// exists is skipped, so an interrupted run resumes and prepared inputs (packed reads) are reused.

struct Sample { std::string name, path; bool normal; };

class Run {
public:
    explicit Run(const Args& a);
    int go();

private:
    const Args& a;
    std::string dir, self, ref, threads;
    std::vector<Sample> smp;
    std::vector<std::string> tumors;
    std::string normal;
    std::string path(const std::string& f) const { return dir + "/" + f; }
    bool done(const std::string& step) const;
    void step(const std::string& name, const std::vector<std::vector<std::string>>& jobs);
    static int child(const std::vector<std::string>& argv, const std::string& log, bool append);
    std::string list(const std::string& suffix) const;
};

static double now() {
    timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + t.tv_nsec * 1e-9;
}

static bool exists(const std::string& f) {
    struct stat st;
    return stat(f.c_str(), &st) == 0;
}

static bool endswith(const std::string& s, const char* e) {
    size_t n = std::strlen(e);
    return s.size() >= n && s.compare(s.size() - n, n, e) == 0;
}

Run::Run(const Args& a) : a(a), dir(a.need(a.out, "--out")), ref(a.need(a.ref, "--ref")), threads(std::to_string(a.threads)) {
    char buf[4096];
    ssize_t n = readlink("/proc/self/exe", buf, sizeof buf - 1);
    if (n <= 0) { std::perror("janus run: /proc/self/exe"); std::exit(1); }
    self.assign(buf, n);
    std::string sheet = a.need(a.samples, "--samples"), base = sheet.find('/') == std::string::npos ? "" : sheet.substr(0, sheet.rfind('/') + 1);
    Lines in(sheet.c_str());
    while (in.next()) {
        char name[256], role[64], p[4096];
        if (std::sscanf(in.s, "%255s %63s %4095s", name, role, p) != 3 || !std::strcmp(name, "name")) continue;
        bool nm = !std::strcmp(role, "normal");
        if (!nm && std::strcmp(role, "tumor")) { std::fprintf(stderr, "janus run: role of %s must be normal or tumor\n", name); std::exit(2); }
        smp.push_back({name, p[0] == '/' ? std::string(p) : base + p, nm});   // relative paths are relative to the sheet
        (nm ? normal : tumors.emplace_back()) = name;
    }
    size_t normals = std::count_if(smp.begin(), smp.end(), [](const Sample& s) { return s.normal; });
    if (normals != 1 || tumors.empty()) { std::fprintf(stderr, "janus run: the sample sheet needs one normal and at least one tumor\n"); std::exit(2); }
    std::stable_partition(smp.begin(), smp.end(), [](const Sample& s) { return s.normal; });
    mkdir(dir.c_str(), 0755);
}

bool Run::done(const std::string& s) const {
    return exists(path("done." + s));
}

int Run::child(const std::vector<std::string>& argv, const std::string& log, bool append) {
    pid_t pid = fork();
    if (pid == 0) {
        FILE* f = std::fopen(log.c_str(), append ? "a" : "w");
        dup2(fileno(f), 1);
        dup2(fileno(f), 2);
        std::vector<char*> av;
        for (const std::string& s : argv) av.push_back((char*)s.c_str());
        av.push_back(nullptr);
        execvp(av[0], av.data());
        std::fprintf(stderr, "cannot run %s\n", av[0]);
        _exit(127);
    }
    int st;
    waitpid(pid, &st, 0);
    return WIFEXITED(st) ? WEXITSTATUS(st) : 128;
}

// one step: the jobs run concurrently, each job a chain of commands; the measurement covers all of them
void Run::step(const std::string& name, const std::vector<std::vector<std::string>>& jobs) {
    if (done(name)) return;
    std::string log = path("log." + name);
    std::fclose(std::fopen(log.c_str(), "w"));
    double t0 = now();
    std::vector<pid_t> pids;
    for (const std::vector<std::string>& job : jobs) {
        pid_t pid = fork();
        if (pid == 0) {
            std::vector<std::string> cmd;
            for (const std::string& s : job) {
                if (s != "&&") { cmd.push_back(s); continue; }
                if (child(cmd, log, true)) _exit(1);
                cmd.clear();
            }
            _exit(cmd.empty() ? 0 : child(cmd, log, true) != 0);
        }
        pids.push_back(pid);
    }
    rusage total = {};
    bool ok = true;
    for (pid_t pid : pids) {
        int st;
        rusage ru;
        wait4(pid, &st, 0, &ru);
        ok &= WIFEXITED(st) && WEXITSTATUS(st) == 0;
        total.ru_utime.tv_sec += ru.ru_utime.tv_sec; total.ru_utime.tv_usec += ru.ru_utime.tv_usec;
        total.ru_stime.tv_sec += ru.ru_stime.tv_sec; total.ru_stime.tv_usec += ru.ru_stime.tv_usec;
        total.ru_maxrss = std::max(total.ru_maxrss, ru.ru_maxrss);
    }
    double wall = now() - t0, user = total.ru_utime.tv_sec + total.ru_utime.tv_usec * 1e-6, sys = total.ru_stime.tv_sec + total.ru_stime.tv_usec * 1e-6;
    std::fprintf(stderr, "[run] %s %.2f s, %.2f CPU s, %ld MB\n", name.c_str(), wall, user + sys, total.ru_maxrss / 1024);
    if (!ok) {
        std::fprintf(stderr, "[run] %s failed; see %s\n", name.c_str(), log.c_str());
        std::exit(1);
    }
    FILE* f = std::fopen(path("steps.tsv").c_str(), "a");
    if (!std::ftell(f)) std::fprintf(f, "step\twall_s\tuser_s\tsys_s\tpeak_rss_kib\n");
    std::fprintf(f, "%s\t%.2f\t%.2f\t%.2f\t%ld\n", name.c_str(), wall, user, sys, total.ru_maxrss);
    std::fclose(f);
    std::fclose(std::fopen(path("done." + name).c_str(), "w"));
}

std::string Run::list(const std::string& suffix) const {
    std::string s;
    for (const Sample& x : smp) s += (s.empty() ? "" : ",") + path(x.name + suffix);
    return s;
}

int Run::go() {
    const std::string J = self, T = threads;
    for (const Sample& s : smp) {
        std::string svb = path(s.name + ".svb");
        if (endswith(s.path, ".svb")) {
            if (!exists(svb)) { symlink(s.path.c_str(), svb.c_str()); symlink((s.path + ".idx").c_str(), (svb + ".idx").c_str()); }
            continue;
        }
        step("pack_" + s.name, {{J, "pack", "--reads", s.path, "--out", svb, "--threads", T}});
    }
    for (const Sample& s : smp) step("extract_" + s.name, {{J, "extract", "--reads", path(s.name + ".svb"), "--threads", T}});
    step("count_" + normal, {{J, "count", "--reads", path(normal + ".svb"), "--out", path(normal + ".cnt"), "--threads", T, "--support", "0", "--shallow", std::to_string(a.shallow)}});
    for (const std::string& t : tumors) step("count_" + t, {{J, "count", "--reads", path(t + ".svb"), "--out", path(t + ".cnt"), "--threads", T, "--normal", path(normal + ".cnt"), "--clonal", std::to_string(a.clonal)}});
    std::vector<int> mins;
    for (const std::string& t : tumors) {
        Input c(path(t + ".cnt"));
        std::fclose(c.f);
        mins.push_back(c.h.resv[2]);
    }
    std::string M = std::to_string(*std::min_element(mins.begin(), mins.end())), svb = list(".svb");
    step("graph", {{J, "graph", "--counts", list(".cnt"), "--out", path("graph.gph"), "--support", M}});
    step("prune", {{J, "prune", "--graph", path("graph.gph"), "--reads", svb, "--out", path("prune"), "--support", M, "--threads", T}});
    step("select", {{J, "select", "--graph", path("graph.gph"), "--sites", path("prune.sites"), "--reads", svb, "--out", path("ev"), "--threads", T}});
    const std::string fa = ref + ".fa";
    if (tumors.size() == 1) {
        step("place", {{J, "place", "--pieces", path("ev"), "--index", ref + ".cidx", "--ref", fa, "--threads", T}});
        step("local", {{"minimap2", "-cx", "map-hifi", "--secondary=no", "-t", T, "-o", path("local.paf"), path("ev.win.fa"), path("ev.on.fa")}});
        step("lift", {{J, "lift", "--paf", path("local.paf"), "--ref", fa, "--out", path("aln.paf")}});
    } else if (!done("place")) {
        // the pieces of every tumor are placed and aligned on their own, up to four tumors at a time
        std::vector<FILE*> parts(tumors.size() + 1, nullptr);
        Lines in(path("ev.fa").c_str());
        FILE* cur = nullptr;
        while (in.next()) {
            if (in.s[0] == '>') {
                size_t s = std::strtoul(in.s + 2, nullptr, 10);
                if (!parts[s]) parts[s] = std::fopen(path("p" + std::to_string(s) + ".fa").c_str(), "w");
                cur = parts[s];
            }
            std::fputs(in.s, cur);
        }
        for (FILE* f : parts) if (f) std::fclose(f);
        size_t jobs = std::min<size_t>(tumors.size(), 4);
        std::string tj = std::to_string(std::max(1, a.threads / (int)jobs));
        std::vector<std::vector<std::string>> chains(jobs);
        for (size_t s = 1; s <= tumors.size(); s++) {
            if (!parts[s]) continue;
            std::string p = path("p" + std::to_string(s));
            std::vector<std::string>& c = chains[(s - 1) % jobs];
            if (!c.empty()) c.push_back("&&");
            std::vector<std::string> cmd = {J, "place", "--pieces", p, "--index", ref + ".cidx", "--ref", fa, "--threads", tj, "&&",
                                            "minimap2", "-cx", "map-hifi", "--secondary=no", "-t", tj, "-o", p + ".local.paf", p + ".win.fa", p + ".on.fa", "&&",
                                            J, "lift", "--paf", p + ".local.paf", "--ref", fa, "--out", p + ".aln.paf"};
            c.insert(c.end(), cmd.begin(), cmd.end());
        }
        step("place", chains);
        FILE* aln = std::fopen(path("aln.paf").c_str(), "w");
        for (size_t s = 1; s <= tumors.size(); s++) {
            if (!parts[s]) continue;
            Lines l(path("p" + std::to_string(s) + ".aln.paf").c_str());
            while (l.next()) std::fputs(l.s, aln);
        }
        std::fclose(aln);
    }
    step("clip", {{J, "clip", "--pieces", path("ev"), "--paf", path("aln.paf"), "--out", path("ev.gaps.fa")}});
    step("gaps", {{"minimap2", "-cx", "sr", "--secondary=no", "-t", T, "-o", path("ev.gaps.raw.paf"), ref + ".sr.mmi", path("ev.gaps.fa")}});
    step("patch", {{J, "patch", "--pieces", path("ev"), "--paf", path("ev.gaps.raw.paf"), "--aln", path("aln.paf"), "--out", path("ev.paf")}});
    step("call", {{J, "call", "--pieces", path("ev"), "--reads", svb, "--ref", fa, "--graph", path("graph.gph"), "--out", path("calls.vcf"), "--vntr", ref + ".trf.bed", "--threads", T, "--support", M}});
    if (!a.model) { std::fprintf(stderr, "[run] done without a model: %s\n", path("calls.vcf").c_str()); return 0; }
    std::string model = a.model;
    if (tumors.size() == 1) step("classify", {{J, "classify", "--calls", path("calls.vcf"), "--hist", path(tumors[0] + ".cnt.hist"), "--model", model, "--out", path("janus.vcf")}});
    else {
        std::string hists, ms;
        for (size_t k = 0; k < tumors.size(); k++) { hists += (k ? "," : "") + path(tumors[k] + ".cnt.hist"); ms += (k ? "," : "") + std::to_string(mins[k]); }
        step("joint", {{J, "joint", "--calls", path("calls.vcf"), "--hist", hists, "--minimum", ms, "--model", model, "--out", path("janus.vcf")}});
    }
    std::fprintf(stderr, "[run] done: %s\n", path("janus.vcf").c_str());
    return 0;
}

int cmd_run(const Args& a) {
    return Run(a).go();
}
