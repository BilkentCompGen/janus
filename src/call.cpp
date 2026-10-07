#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "call.h"
#include "io.h"
#include "keys.h"
#include "reads.h"

static std::vector<uint32_t> nums(const char* s) {
    std::vector<uint32_t> v;
    for (char* p = (char*)s; *p && *p != '\t' && *p != '\n'; p += *p == ',') v.push_back(std::strtoul(p, &p, 10));
    return v;
}

void Spot::parse(const char* line, const Graph& g) {
    char name[256], cn[4096], ls[4096], rs[4096];
    unsigned long key, sample, index;
    std::sscanf(line, "%255s %u %u %*u %4095s %lu %u %4095s %4095s %u", name, &a, &b, cn, &key, &s, ls, rs, &lo);
    nm = name;
    n = nums(cn);
    t = n[s];
    tonly = nbr = 0;
    lk = 0;
    m.assign(g.S * KINDS, 0);
    dl.assign(g.S, std::vector<int32_t>());
    for (uint32_t i = 0; i < NA; i++)
        for (uint32_t j = 0; j < NA; j++) d[i][j] = UNSET;
    std::sscanf(name, "s%lu_%lu_", &sample, &index);
    read = read_id(sample, index);
    const char* p[2] = {ls, rs};
    for (int k = 0; k < 2; k++)
        for (const char* q = p[k]; *q && *q != '-' && *q != '\t' && *q != '\n';) {
            (k ? r : l).push_back(std::strtoull(q, (char**)&q, 10));
            std::strtoul(q + 1, (char**)&q, 10);
            if (*q == ',') q++;
        }
    uint64_t u = Graph::from(g.edges[Combos::in(key)]), w = Graph::to(g.edges[Combos::out(key)]);
    for (uint64_t k = 0; k < g.S; k++)
        dep.push_back(l.empty() || r.empty() ? std::max(g.reads(k, u), g.reads(k, w)) : std::min(g.reads(k, l[0]), g.reads(k, r[0])));
}

static const uint32_t PM = 20;
static const uint32_t PLIM = 15000;
static const uint32_t PDUP = 2000;
static const uint32_t PWIDE = 10000;
static const uint32_t PLOC = 1000;

class CallPass : public Walker {
public:
    CallPass(std::vector<Spot>& w, std::vector<Probe>& pr, const Graph& g, int t);
    void run(const ReadFiles& files);

private:
    std::vector<Spot>& w;
    std::vector<Probe>& pr;
    std::vector<uint32_t> av, ai, pv, pi;
    std::vector<std::pair<uint64_t, uint32_t>> byread;
    std::vector<uint8_t> isa, nz, br, gd, isp;
    std::vector<std::vector<uint32_t>> sb;
    std::vector<uint64_t> rs;
    std::vector<std::vector<uint64_t>> hb;
    bool measuring;
    void visit(Chain& c, const uint8_t* rec, uint32_t m, uint32_t rc, uint64_t rid) override;
    void hits(const Chain& c, uint32_t m, std::vector<uint64_t>& h) const;
    void own(Chain& c, uint32_t m, uint32_t rc, uint64_t rid);
    void place(Chain& c, uint32_t m, Spot& e);
    void measure(Chain& c, uint32_t m);
    void probe(Chain& c, uint32_t m);
};

CallPass::CallPass(std::vector<Spot>& w, std::vector<Probe>& pr, const Graph& g, int t) : Walker(t), w(w), pr(pr), isa(g.K, 0), nz(g.K, 0), br(g.K, 0), gd(g.K, 0), isp(g.K, 0), sb(t), hb(t), measuring(false) {
    for (uint64_t v = 0; v < g.K; v++) nz[v] = g.reads(0, v) == 0;
    std::vector<uint32_t> full;
    for (uint64_t v = 0; v < g.K; v++) {
        bool all = true;
        for (uint64_t s = 0; s < g.S; s++) all &= g.reads(s, v) > 0;
        if (all) full.push_back(g.reads(0, v));
    }
    uint32_t med = 0;
    if (!full.empty()) {
        std::nth_element(full.begin(), full.begin() + full.size() / 2, full.end());
        med = full[full.size() / 2];
    }
    for (uint64_t v = 0; v < g.K; v++) {
        gd[v] = g.reads(0, v) > 0 && g.reads(0, v) <= 2 * med;
        for (uint64_t s = 1; s < g.S; s++) gd[v] &= g.reads(s, v) > 0;
    }
    std::vector<uint8_t> in(g.K, 0), on(g.K, 0);
    for (uint64_t e : g.edges) {
        uint64_t u = Graph::from(e), v = Graph::to(e);
        if (on[u] < 2) on[u]++;
        if (in[v] < 2) in[v]++;
    }
    for (uint64_t v = 0; v < g.K; v++) br[v] = in[v] > 1 || on[v] > 1;
    std::vector<std::pair<uint32_t, uint32_t>> an;
    for (uint32_t i = 0; i < w.size(); i++) {
        for (uint32_t k = 0; k < w[i].l.size() && k < NA; k++) an.push_back({w[i].l[k], i << 3 | k});
        for (uint32_t k = 0; k < w[i].r.size() && k < NA; k++) an.push_back({w[i].r[k], i << 3 | 4 | k});
        rs.push_back(w[i].read);
    }
    std::sort(an.begin(), an.end());
    for (std::pair<uint32_t, uint32_t>& p : an) { av.push_back(p.first); ai.push_back(p.second); isa[p.first] = 1; }
    std::sort(rs.begin(), rs.end());
    rs.erase(std::unique(rs.begin(), rs.end()), rs.end());
    index(g.hashes);
}

void CallPass::visit(Chain& c, const uint8_t*, uint32_t m, uint32_t rc, uint64_t rid) {
    if (measuring) { measure(c, m); probe(c, m); }
    else own(c, m, rc, rid);
}

void CallPass::place(Chain& c, uint32_t m, Spot& e) {
    std::vector<uint32_t>& ids = sb[c.slot];
    ids.assign(c.id.begin(), c.id.begin() + m);
    std::sort(ids.begin(), ids.end());
    for (uint32_t x : e.pb)
        for (uint32_t pass = 0; pass < 2 && !pr[x].set; pass++) {
            Probe& b = pr[x];
            if (pass && b.qa == b.ta && b.qb == b.tb) continue;
            int64_t ca = (int64_t)e.lo + (pass ? b.ta : b.qa), cb = (int64_t)e.lo + (pass ? b.tb : b.qb);
            uint32_t nl = 0, nr = 0;
            b.absent = 0;
            for (uint32_t i = 0; i < NA; i++)
                for (uint32_t j = 0; j < NA; j++) b.d[i][j] = UNSET;
            for (int64_t i = (int64_t)m - 1; i >= 0 && nl < NA; i--) {
                if ((int64_t)c.end(i) > ca) continue;
                if (ca - (int64_t)c.start(i) > (int64_t)PLIM) break;
                uint32_t v = c.id[i];
                b.absent += nz[v];
                if (!gd[v] || std::upper_bound(ids.begin(), ids.end(), v) - std::lower_bound(ids.begin(), ids.end(), v) != 1) continue;
                b.v[0][nl] = v;
                b.at[0][nl++] = i;
            }
            for (uint32_t i = 0; i < m && nr < NA; i++) {
                if ((int64_t)c.start(i) < cb) continue;
                if ((int64_t)c.start(i) - cb > (int64_t)PLIM) break;
                uint32_t v = c.id[i];
                b.absent += nz[v];
                if (!gd[v] || std::upper_bound(ids.begin(), ids.end(), v) - std::lower_bound(ids.begin(), ids.end(), v) != 1) continue;
                b.v[1][nr] = v;
                b.at[1][nr++] = i;
            }
            if (!nl || !nr) continue;
            for (uint32_t i = 0; i < nl; i++)
                for (uint32_t j = 0; j < nr; j++) b.d[i][j] = (int32_t)c.start(b.at[1][j]) - (int32_t)c.start(b.at[0][i]);
            for (uint32_t i = nl; i < NA; i++) b.v[0][i] = NONE;
            for (uint32_t i = nr; i < NA; i++) b.v[1][i] = NONE;
            b.gap = ca - c.end(b.at[0][0]) + c.start(b.at[1][0]) - cb;
            for (uint32_t i = b.at[0][nl - 1] + 1; i < b.at[1][nr - 1]; i++) b.mid.push_back({c.id[i], i});
            std::sort(b.mid.begin(), b.mid.end());
            b.set = 1 + pass;
            break;
        }
}

void CallPass::probe(Chain& c, uint32_t m) {
    std::vector<uint64_t>& h = hb[c.slot];
    h.clear();
    for (uint32_t i = 0; i < m; i++) {
        uint32_t v = c.id[i];
        if (!isp[v]) continue;
        for (size_t j = std::lower_bound(pv.begin(), pv.end(), v) - pv.begin(); j < pv.size() && pv[j] == v; j++) h.push_back((uint64_t)pi[j] << 32 | i);
    }
    std::sort(h.begin(), h.end());
    for (size_t i = 0, j; i < h.size(); i = j) {
        uint32_t s = h[i] >> 35, ps[2][NA], oc[2][NA];
        for (uint32_t k = 0; k < NA; k++) oc[0][k] = oc[1][k] = 0;
        for (j = i; j < h.size() && h[j] >> 35 == s; j++) {
            uint32_t k = h[j] >> 32 & 7;
            ps[k >= 4][k & 3] = (uint32_t)h[j];
            oc[k >= 4][k & 3]++;
        }
        Probe& b = pr[s];
        for (uint32_t x = 0; x < NA; x++)
            for (uint32_t y = 0; y < NA; y++) {
                if (oc[0][x] != 1 || oc[1][y] != 1 || b.d[x][y] == UNSET) continue;
                int32_t dl = (int32_t)c.start(ps[1][y]) - (int32_t)c.start(ps[0][x]) - b.d[x][y], pd = 99;
                uint32_t i0 = ps[0][x], i1 = ps[1][y], o0 = b.at[0][x], o1 = b.at[1][y];
                if (i0 < i1) {
                    std::vector<uint32_t>& rd = sb[c.slot];
                    rd.assign(c.id.begin() + i0 + 1, c.id.begin() + i1);
                    std::sort(rd.begin(), rd.end());
                    pd = 0;
                    for (const std::pair<uint32_t, uint32_t>& z : b.mid)
                        if (z.second > o0 && z.second < o1 && !std::binary_search(rd.begin(), rd.end(), z.first)) pd++;
                    for (uint32_t v : rd) {
                        bool in = false;
                        for (std::vector<std::pair<uint32_t, uint32_t>>::const_iterator z = std::lower_bound(b.mid.begin(), b.mid.end(), std::make_pair(v, 0u)); z != b.mid.end() && z->first == v && !in; z++) in = z->second > o0 && z->second < o1;
                        pd += !in;
                    }
                }
                lock(b.lk);
                b.dl[sample * NA * NA + x * NA + y].push_back({dl, pd});
                unlock(b.lk);
            }
    }
}

void CallPass::hits(const Chain& c, uint32_t m, std::vector<uint64_t>& h) const {
    h.clear();
    for (uint32_t i = 0; i < m; i++) {
        uint32_t v = c.id[i];
        if (!isa[v]) continue;
        for (size_t j = std::lower_bound(av.begin(), av.end(), v) - av.begin(); j < av.size() && av[j] == v; j++) h.push_back((uint64_t)ai[j] << 32 | i);
    }
    std::sort(h.begin(), h.end());
}

void CallPass::own(Chain& c, uint32_t m, uint32_t rc, uint64_t rid) {
    if (!rc)
        for (std::vector<std::pair<uint64_t, uint32_t>>::const_iterator z = std::lower_bound(byread.begin(), byread.end(), std::make_pair(read_id(sample, rid), 0u)); z != byread.end() && z->first == read_id(sample, rid); z++)
            place(c, m, w[z->second]);
    std::vector<uint64_t>& h = hb[c.slot];
    hits(c, m, h);
    for (size_t i = 0, j; i < h.size(); i = j) {
        uint32_t s = h[i] >> 35;
        for (j = i; j < h.size() && h[j] >> 35 == s; j++);
        Spot* e = &w[s];
        if (e->read != read_id(sample, rid)) continue;
        for (size_t x = i; x < j; x++)
            for (size_t y = i; y < j; y++) {
                uint32_t la = h[x] >> 32 & 7, rb = h[y] >> 32 & 7, pa = (uint32_t)h[x], pb = (uint32_t)h[y];
                if (la >= 4 || rb < 4 || pa == pb) continue;
                uint32_t li = la, ri = rb - 4, al = c.start(pa), ar = c.start(pb);
                if (e->d[li][ri] != UNSET) continue;
                e->d[li][ri] = al < ar ? ar - al : al - ar;
                if (li || ri) continue;
                e->tonly = e->nbr = 0;
                for (uint32_t k = std::min(pa, pb); k <= std::max(pa, pb); k++) { e->tonly += nz[c.id[k]]; e->nbr += br[c.id[k]]; }
            }
    }
}

void CallPass::measure(Chain& c, uint32_t m) {
    std::vector<uint64_t>& h = hb[c.slot];
    hits(c, m, h);
    for (size_t i = 0, j; i < h.size(); i = j) {
        uint32_t s = h[i] >> 35, li = NONE, ri = NONE, pl = 0, pr = 0, ps[2][NA];
        for (uint32_t k = 0; k < NA; k++) ps[0][k] = ps[1][k] = NONE;
        for (j = i; j < h.size() && h[j] >> 35 == s; j++) {
            uint32_t k = h[j] >> 32 & 7, pos = (uint32_t)h[j];
            if (k < 4 && li == NONE) { li = k; pl = pos; }
            if (k >= 4 && ri == NONE) { ri = k - 4; pr = pos; }
            if (ps[k >= 4][k & 3] == NONE) ps[k >= 4][k & 3] = pos;
        }
        Spot* e = &w[s];
        if (li == NONE || ri == NONE || pl == pr || e->d[li][ri] == UNSET) continue;
        uint32_t al = c.start(pl), ar = c.start(pr), *cnt = e->m.data() + sample * KINDS;
        if (al > ar) std::swap(al, ar);
        int32_t dl = (int32_t)(ar - al) - e->d[li][ri], ad = dl < 0 ? -dl : dl;
        uint32_t ok[2] = {0, 0};
        for (uint32_t sd = 0; sd < 2; sd++)
            for (uint32_t k = 1; k < NA; k++) {
                if (ps[sd][0] == NONE || ps[sd][k] == NONE) continue;
                int32_t t0 = e->d[0][0], tk = sd ? e->d[0][k] : e->d[k][0];
                if (t0 == UNSET || tk == UNSET) continue;
                int32_t sp = std::labs(t0 - tk), got = std::labs((int32_t)c.start(ps[sd][k]) - (int32_t)c.start(ps[sd][0]));
                ok[sd] += std::labs(got - sp) <= (int32_t)DSAME;
            }
        if (ok[0] >= NA - 1 && ok[1] >= NA - 1) {
            __atomic_fetch_add(&cnt[CONSISTENT], 1, __ATOMIC_RELAXED);
            if (ad < (int32_t)DSAME) __atomic_fetch_add(&cnt[CONSISTENT_SAME], 1, __ATOMIC_RELAXED);
        }
        __atomic_fetch_add(&cnt[SPAN], 1, __ATOMIC_RELAXED);
        if (ad < (int32_t)DSAME) __atomic_fetch_add(&cnt[SAME], 1, __ATOMIC_RELAXED);
        lock(e->lk);
        e->dl[sample].push_back(dl);
        unlock(e->lk);
    }
}

void CallPass::run(const ReadFiles& files) {
    for (uint32_t i = 0; i < w.size(); i++)
        if (!w[i].pb.empty()) byread.push_back({w[i].read, i});
    std::sort(byread.begin(), byread.end());
    fetch(files, rs);
    std::fprintf(stderr, "[call] %zu sites on %zu tumor reads measured\n", w.size(), rs.size());
    std::vector<std::pair<uint32_t, uint32_t>> an;
    for (uint32_t i = 0; i < pr.size(); i++) {
        if (!pr[i].set) continue;
        for (uint32_t s = 0; s < 2; s++)
            for (uint32_t k = 0; k < NA; k++)
                if (pr[i].v[s][k] != NONE) an.push_back({pr[i].v[s][k], i << 3 | s << 2 | k});
    }
    std::sort(an.begin(), an.end());
    for (std::pair<uint32_t, uint32_t>& p : an) { pv.push_back(p.first); pi.push_back(p.second); isp[p.first] = 1; }
    measuring = true;
    stream(files.names, 0);
    for (Probe& b : pr)
        for (std::vector<std::pair<int32_t, int32_t>>& d : b.dl) std::sort(d.begin(), d.end());
    for (Spot& e : w)
        for (std::vector<int32_t>& d : e.dl) std::sort(d.begin(), d.end());
}

static const uint32_t ROW = 3000;
static const uint32_t TOL = 50;

void Call::sites() {
    Graph g(a.graph);
    {
        Lines in((pre + ".tsv").c_str());
        while (in.next()) {
            Spot r;
            r.parse(in.s, g);
            bypc[r.nm].push_back(w.size());
            w.push_back(r);
        }
    }
    alignments();
    probes(g.S);
    CallPass pass(w, pbs, g, a.threads);
    g = Graph();
    pass.run(ReadFiles(a.reads));
    std::vector<uint32_t> dep;
    for (const Spot& e : w) dep.push_back(e.dep[e.s]);
    std::sort(dep.begin(), dep.end());
    med = dep.empty() ? 0 : dep[dep.size() / 2];
}

void Call::probes(uint64_t S) {
    for (std::pair<const std::string, std::vector<uint32_t>>& pc : bypc) {
        if (!al.count(pc.first)) continue;
        for (Event& e : evs[pc.first] = events(al[pc.first])) {
            uint32_t u = NONE;
            for (uint32_t i : pc.second)
                if (e.q + ROW >= w[i].a && e.q <= w[i].b + ROW && (u == NONE || w[i].t > w[u].t)) u = i;
            e.pb = 0;
            if (u == NONE) continue;
            bool ins = e.ty == "INS", del = e.ty == "DEL", dup = e.ty == "DUP";
            uint32_t len = (ins || del || dup) && e.c1 == e.c2 ? e.p2 - e.p1 : 0, vs, ve, lx = 0, rx = 0, end = del || dup ? e.p2 : e.p1;
            if ((ins || del || dup) && bed.span(e.c1, e.p1, &vs, &ve)) {
                lx = e.p1 > vs ? e.p1 - vs : 0;
                rx = ve > end ? ve - end : 0;
            }
            if (lx + rx > PWIDE) { lx = std::min(lx, PLOC); rx = std::min(rx, PLOC); }
            if (dup && len <= PDUP) { lx += len; rx += len; }
            if (e.st == '-') std::swap(lx, rx);
            Probe b;
            b.qa = (int32_t)e.q - (int32_t)lx - (int32_t)PM;
            b.qb = e.q + (ins ? len : 0) + rx + PM;
            b.ta = (int32_t)e.q - (int32_t)PM;
            b.tb = e.q + (ins ? len : 0) + PM;
            b.lk = b.set = 0;
            b.gap = b.absent = 0;
            b.dl.assign(S * NA * NA, std::vector<std::pair<int32_t, int32_t>>());
            for (uint32_t i = 0; i < NA; i++) {
                b.v[0][i] = b.v[1][i] = NONE;
                for (uint32_t j = 0; j < NA; j++) b.d[i][j] = UNSET;
            }
            w[u].pb.push_back(pbs.size());
            pbs.push_back(b);
            e.pb = pbs.size();
        }
    }
}

void Call::alignments() {
    Lines in((pre + ".paf").c_str());
    Paf p;
    while (in.next()) {
        p.parse(in.s);
        const char *cg = std::strstr(in.s, "cg:Z:"), *de = std::strstr(in.s, "de:f:");
        al[p.q].push_back({p.qs, p.qe, p.ts, p.te, p.mq, p.st[0], p.t, cg ? std::string(cg + 5, std::strcspn(cg + 5, "\t\n")) : "", de ? std::strtof(de + 5, nullptr) : 0});
    }
}

void Call::sequences() {
    for (std::pair<const std::string, std::vector<uint32_t>>& pc : bypc)
        if (al.count(pc.first))
            for (const Event& e : evs[pc.first])
                if (e.ty == "INS" && e.p2 - e.p1 >= NEAR) sq[pc.first];
    Lines in((pre + ".fa").c_str());
    std::string* cur = nullptr;
    while (in.next())
        if (in.s[0] == '>') {
            std::unordered_map<std::string, std::string>::iterator it = sq.find(std::string(in.s + 1, in.n - 2));
            cur = it == sq.end() ? nullptr : &it->second;
        } else if (cur) cur->assign(in.s, in.n - 1);
}

static bool primary(const std::string& c) {
    return c.size() <= 5 && c.find('_') == std::string::npos && c != "chrM" && c != "M";
}

static bool by_event(const Record& x, const Record& y) {
    return x.e.ty != y.e.ty ? x.e.ty < y.e.ty : x.e.c1 != y.e.c1 ? x.e.c1 < y.e.c1 : x.e.c2 != y.e.c2 ? x.e.c2 < y.e.c2 : x.e.o != y.e.o ? x.e.o < y.e.o : x.e.p1 != y.e.p1 ? x.e.p1 < y.e.p1 : x.e.p2 != y.e.p2 ? x.e.p2 < y.e.p2 : x.e.q != y.e.q ? x.e.q < y.e.q : x.nm < y.nm;
}

static bool near(const Record& x, const Record& y) {
    uint32_t lx = x.e.p2 - x.e.p1, ly = y.e.p2 - y.e.p1;
    if ((x.e.ty == "INS" || x.e.ty == "DEL") && x.e.c1 == x.e.c2 && std::max(lx, ly) - std::min(lx, ly) > std::max<uint32_t>(10, std::max(lx, ly) / 10)) return false;
    return x.e.ty == y.e.ty && x.e.o == y.e.o && x.e.c1 == y.e.c1 && x.e.c2 == y.e.c2 && x.e.p1 + TOL >= y.e.p1 && y.e.p1 + TOL >= x.e.p1 && x.e.p2 + TOL >= y.e.p2 && y.e.p2 + TOL >= x.e.p2;
}

std::vector<Record> Call::attach(size_t* found) {
    std::vector<Record> som, res;
    for (std::pair<const std::string, std::vector<uint32_t>>& pc : bypc) {
        if (!al.count(pc.first)) continue;
        for (const Event& e : evs[pc.first]) {
            const Spot* r = nullptr;
            uint32_t u = 0;
            for (uint32_t i : pc.second)
                if (e.q + ROW >= w[i].a && e.q <= w[i].b + ROW && (!r || w[i].t > r->t)) { r = &w[i]; u = i; }
            if (!r || !primary(e.c1) || !primary(e.c2) || r->t < a.support) continue;
            Row row(*r, u);
            som.push_back({e, pc.first, 1, row, misfit(row, e), row.nsame, row.nspan, e.q < r->a ? r->a - e.q : e.q > r->b ? e.q - r->b : 0});
        }
    }
    std::sort(som.begin(), som.end(), by_event);
    for (const Record& v : som) {
        if (res.empty() || !near(res.back(), v)) { res.push_back(v); continue; }
        Record& q = res.back();
        bool more = top(v.r.n) > top(q.r.n);
        Row m = more ? v.r : q.r;
        if (more) { q.e = v.e; q.nm = v.nm; }
        if (v.nf < q.nf || (v.nf == q.nf && v.ms > q.ms)) { q.nf = v.nf; q.ms = v.ms; q.mn = v.mn; }
        m.nsame = std::max(v.r.nsame, q.r.nsame);
        m.nspan = std::max(v.r.nspan, q.r.nspan);
        m.nspc = std::max(v.r.nspc, q.r.nspc);
        m.n20c = std::max(v.r.n20c, q.r.n20c);
        m.ncl = std::max(v.r.ncl, q.r.ncl);
        q.r = m;
        q.pc++;
    }
    *found = som.size();
    return res;
}

int Call::run() {
    sites();
    sequences();
    size_t found;
    std::vector<Record> res = attach(&found);
    write(res, found);
    return 0;
}

int cmd_call(const Args& a) {
    return Call(a).run();
}
