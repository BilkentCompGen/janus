#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "cmd.h"
#include "io.h"
#include "model.h"

static const int64_t DUPT = 500;

static void fail(const char* what, const char* path) {
    std::fprintf(stderr, "janus: %s: %s\n", what, path);
    std::exit(1);
}

Model::Model(const char* path) {
    Lines in(path);
    char tag[32];
    int seeds = 0, k, n;
    if (!in.next() || std::strncmp(in.s, "janus-model 1", 13)) fail("not a janus model", path);
    in.next();
    if (std::sscanf(in.s, "features %d seeds %d cut %lf", &features, &seeds, &cut) != 3) fail("bad model header", path);
    for (const char* p = std::strstr(in.s, "cuts "); p && *p && std::strncmp(p, "train", 5);) {
        double c;
        if (std::sscanf(p, "%d:%lf", &k, &c) == 2) cuts[k] = c;
        p = std::strchr(p, ' ');
        if (p) p++;
    }
    baseline.resize(seeds);
    trees.resize(seeds);
    for (int s = 0; s < seeds; s++) {
        int t = 0;
        in.next();
        if (std::sscanf(in.s, "seed %d baseline %lf trees %d", &k, &baseline[s], &t) != 3 || k != s) fail("bad seed line", path);
        trees[s].resize(t);
        for (std::vector<Node>& tree : trees[s]) {
            in.next();
            if (std::sscanf(in.s, "%s %d", tag, &n) != 2 || std::strcmp(tag, "tree")) fail("bad tree line", path);
            tree.resize(n);
            for (Node& d : tree) {
                int ml;
                in.next();
                if (std::sscanf(in.s, "%lf %d %lf %d %d %d", &d.value, &d.feature, &d.threshold, &ml, &d.left, &d.right) != 6) fail("bad node line", path);
                d.missing_left = ml;
            }
        }
    }
}

double Model::predict(const double* x) const {
    double sum = 0;
    for (size_t s = 0; s < trees.size(); s++) {
        double raw = baseline[s];
        for (const std::vector<Node>& tree : trees[s]) {
            const Node* d = &tree[0];
            while (d->feature >= 0) d = &tree[std::isnan(x[d->feature]) ? (d->missing_left ? d->left : d->right) : x[d->feature] <= d->threshold ? d->left : d->right];
            raw += d->value;
        }
        sum += 1 / (1 + std::exp(-raw));
    }
    return sum / trees.size();
}

double Model::threshold(int minimum) const {
    std::map<int, double>::const_iterator it = cuts.find(minimum);
    return it == cuts.end() ? cut : it->second;
}

static std::vector<int> ints(const std::string& s) {
    std::vector<int> v;
    for (const char* p = s.c_str(); *p; p += *p == ',') v.push_back(std::strtol(p, (char**)&p, 10));
    return v;
}

Vcf::Vcf(const char* path) {
    Lines in(path);
    while (in.next()) {
        if (in.s[0] == '#') { header.push_back(in.s); continue; }
        Sv q;
        q.line = in.s;
        std::vector<std::string> f;
        for (size_t a = 0, b; a <= q.line.size(); a = b + 1) { b = std::min(q.line.find('\t', a), q.line.size()); f.push_back(q.line.substr(a, b - a)); }
        std::map<std::string, std::string> info;
        for (size_t a = 0, b; a <= f[7].size(); a = b + 1) {
            b = std::min(f[7].find(';', a), f[7].size());
            size_t e = f[7].find('=', a);
            if (e < b) info[f[7].substr(a, e - a)] = f[7].substr(e + 1, b - e - 1);
        }
        q.c1 = f[0];
        q.p1 = std::strtoll(f[1].c_str(), nullptr, 10);
        size_t br = f[4].find_first_of("[]");
        if (br != std::string::npos) {
            size_t colon = f[4].find(':', br);
            q.c2 = f[4].substr(br + 1, colon - br - 1);
            q.p2 = std::strtoll(f[4].c_str() + colon + 1, nullptr, 10);
        } else {
            q.c2 = q.c1;
            q.p2 = q.p1 + std::strtoll(info["SVLEN"].c_str(), nullptr, 10);
        }
        q.type = info["SVTYPE"];
        q.o = std::strtoul(info["O"].c_str(), nullptr, 10);
        for (const char* p = info["F"].c_str(); *p; p += *p == ',') q.x.push_back(std::strtod(p, (char**)&p));
        q.own = 0;
        for (size_t k = 10; k < f.size(); k++)
            if (!f[k].compare(0, 3, "0/1")) { q.own = k - 10; break; }
        q.same = ints(info["SAME"]);
        q.span = ints(info["SPAN"]);
        q.same.erase(q.same.begin());
        q.span.erase(q.span.begin());
        r.push_back(q);
    }
}

std::vector<std::string> Vcf::tumors() const {
    std::vector<std::string> v;
    std::string h = header.back();
    h.erase(h.find_last_not_of("\r\n") + 1);
    size_t col = 0;
    for (size_t a = 0, b; a <= h.size(); a = b + 1, col++) {
        b = std::min(h.find('\t', a), h.size());
        if (col >= 10) v.push_back(h.substr(a, b - a));
    }
    return v;
}

int Vcf::minimum() const {
    double lo = r[0].x[1];
    for (const Sv& q : r) lo = std::min(lo, q.x[1]);
    return (int)std::nearbyint(std::pow(2.0, lo) - 1);
}

void Vcf::write(const char* path, const std::vector<double>& pr, const std::vector<uint8_t>& du, const std::vector<uint8_t>& pass, const std::vector<std::string>& extra) const {
    FILE* o = std::fopen(path, "w");
    for (const std::string& h : header) std::fputs(h.c_str(), o);
    for (size_t i = 0; i < r.size(); i++) {
        const std::string& l = r[i].line;
        size_t t[9], at = 0;
        for (int k = 0; k < 9; k++) { t[k] = l.find('\t', at); at = t[k] + 1; }
        std::fprintf(o, "%.*s\t%s\t%.*s;VP=%.4f%s%s", (int)t[5], l.c_str(), du[i] ? "DUPL" : pass[i] ? "PASS" : "LOW", (int)(t[7] - t[6] - 1), l.c_str() + t[6] + 1, pr[i], extra.empty() ? "" : extra[i].c_str(), l.c_str() + t[7]);
    }
    std::fclose(o);
}

static double pct(const std::vector<double>& h, double f, int lo) {
    double total = 0;
    for (size_t i = lo; i < h.size(); i++) total += h[i];
    total = std::max(total, 1.0);
    double c = 0;
    size_t i = lo;
    for (; i < h.size(); i++) {
        c += h[i];
        if (c / total >= f) break;
    }
    return i;
}

Context::Context(const char* path, int mt) {
    std::map<std::string, std::vector<double>> h;
    Lines in(path);
    while (in.next()) {
        char name[64];
        int k = 0;
        std::sscanf(in.s, "%63s%n", name, &k);
        std::vector<double>& v = h[name];
        for (const char* p = in.s + k; *p && *p != '\n'; v.push_back(std::strtod(p, (char**)&p)));
    }
    int mn = 0;
    while (h["shared_normal"][mn] == 0) mn++;
    tm = pct(h["shared_tumor"], 0.5, mt);
    nm = pct(h["shared_normal"], 0.5, mn);
    c[0] = std::log2(tm);
    c[1] = std::log2(nm);
    c[2] = pct(h["tumor_only"], 0.05, mt) / tm;
    c[3] = pct(h["tumor_only"], 0.95, mt) / tm;
}

void inputs(const Sv& r, const Context& c, int mt, std::vector<double>& x) {
    static const int T[] = {1, 3, 39}, N[] = {26, 27, 33, 34, 38, 41, 47, 51};
    x = r.x;
    for (int j : T) x[j] -= std::log2(c.tm + 1);
    for (int j : N) x[j] -= std::log2(c.nm + 1);
    x.push_back(mt);
    x.insert(x.end(), c.c, c.c + 4);
}

std::vector<uint8_t> duplicates(const std::vector<Sv>& r, const std::vector<double>& pr) {
    std::vector<size_t> ix(r.size());
    for (size_t i = 0; i < ix.size(); i++) ix[i] = i;
    std::sort(ix.begin(), ix.end(), [&](size_t a, size_t b) {
        const Sv &x = r[a], &y = r[b];
        return x.c1 != y.c1 ? x.c1 < y.c1 : x.c2 != y.c2 ? x.c2 < y.c2 : x.p1 != y.p1 ? x.p1 < y.p1 : x.p2 != y.p2 ? x.p2 < y.p2 : a < b;
    });
    std::vector<uint8_t> d(r.size(), 0);
    for (size_t a = 0; a < ix.size(); a++) {
        const Sv& x = r[ix[a]];
        for (size_t b = a + 1; b < ix.size(); b++) {
            const Sv& y = r[ix[b]];
            if (y.c1 != x.c1 || y.c2 != x.c2 || y.p1 > x.p1 + DUPT) break;
            if (std::llabs(y.p2 - x.p2) <= DUPT && (x.type != y.type || x.o == y.o)) d[pr[ix[a]] < pr[ix[b]] ? ix[a] : ix[b]] = 1;
        }
    }
    return d;
}

void classify(const Model& m, const std::vector<Sv>& r, const std::vector<size_t>& ix, const char* hist, int mt, std::vector<double>& pr) {
    Context c(hist, mt);
    std::vector<double> x;
    for (size_t i : ix) {
        inputs(r[i], c, mt, x);
        pr[i] = m.predict(x.data());
    }
}

int cmd_classify(const Args& a) {
    Vcf v(a.need(a.calls, "--calls"));
    Model m(a.need(a.model, "--model"));
    std::vector<size_t> ix(v.r.size());
    for (size_t i = 0; i < ix.size(); i++) ix[i] = i;
    int mt = v.minimum();
    std::vector<double> pr(v.r.size());
    classify(m, v.r, ix, a.need(a.hist, "--hist"), mt, pr);
    std::vector<uint8_t> du = duplicates(v.r, pr), pass(v.r.size());
    double cut = m.threshold(mt);
    size_t np = 0, nd = 0;
    for (size_t i = 0; i < pr.size(); i++) { pass[i] = pr[i] >= cut; np += pass[i] && !du[i]; nd += du[i]; }
    v.write(a.need(a.out, "--out"), pr, du, pass, {});
    std::fprintf(stderr, "[classify] %zu PASS, %zu LOW, %zu DUPL of %zu records, threshold %.2f\n", np, v.r.size() - np - nd, nd, v.r.size(), cut);
    return 0;
}
