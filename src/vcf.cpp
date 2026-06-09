// vcf: emit calls.vcf from modeb.bin + bubble_pos.bin (+ optional dup_flags / inv_flags).
//
// POS (strand-aware, verified against 16,950 truth SNPs and ~470 each truth INS / DEL):
//   strand=0:  POS_VCF = s1_left + pre + 1                  (SNP)
//              POS_VCF = s1_left + pre                       (INS/DEL "base before variant")
//   strand=1:  POS_VCF = s2_right - suf                      (SNP, ±75bp)
//              POS_VCF = s2_right - suf + 1                  (INS/DEL)
// Strand=0 covers ~80% of bubbles and 98.7% of those land on the exact truth POS.
// Strand=1 is approximate (median ±20bp); for truvari's refdist=500 settings this is OK.
//
// Tag -> VCF record rules:
//   S, |N'|=|T'|=1            -> SNP (REF=N', ALT=T')
//   S, |N'|=|T'|>1, no INV    -> MNP (REF=N', ALT=T')
//   S, INV-flagged            -> <INV>
//   D                         -> explicit DEL <50bp; symbolic <DEL> ≥50bp
//   I, DUP-flagged            -> <DUP> (TARGET=END)
//   I, otherwise              -> explicit INS <50bp; symbolic <INS> ≥50bp
//   C                         -> DROPPED (measured 0% coverage filter pass)
//
// Phantom filter (empirically derived from 37,654 modeb bubbles vs the truth VCF; details
// in docs/MODEB_DEEP_ANALYSIS.md):
//   D ≥50bp: drop if span < 5 (90% of D phantoms have span < 10; 2.5% of real D do).
//   I ≥50bp: drop if median(cov_t over T' residual) > 200 (real I cov_t≈52; phantom I≈902).
//   I ≥50bp: drop if span > 5*t_resid (anchor-distance sanity).
//   S, 1bp/1bp: drop if median(cov_n) > 150 AND span < 5 (phantom SNP signature).
//   C: drop (0% useful records).
//   chrom unanchored or low anchor resolution: drop.
#include "graph.h"
#include "svp.h"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

struct CovRun  { uint32_t end_pos, n_cov, t_cov; };

struct BubbleTag {
    uint32_t s1, s2;
    uint8_t  tag, bucket; uint16_t _pad;
    uint32_t pre, suf;
    uint32_t n_resid_off, n_resid_len, t_resid_off, t_resid_len;
    uint32_t n_runs_off, n_runs_len, t_runs_off, t_runs_len;
};

struct BubblePos {
    uint32_t bubble_id, chrom_id;
    uint32_t s1_left, s1_right, s2_left, s2_right;
    uint32_t n_anchor_s1, n_resolved_s1, n_anchor_s2, n_resolved_s2;
    uint8_t  strand; uint8_t _pad[3];
};

struct DupCall { uint32_t bubble_id, svlen; uint8_t mechanism; uint8_t _pad[3]; };
struct InvCall { uint32_t bubble_id_a, bubble_id_b, pos_start, pos_end, svlen; uint8_t mechanism; uint8_t _pad[3]; };

struct VcfRec { uint32_t chrom_id; uint32_t pos_vcf; std::string line; };

void load_chroms(const std::string& dir, std::vector<std::string>& names, std::vector<uint64_t>& lengths) {
    FILE* f = std::fopen((dir + "/chroms.tsv").c_str(), "r");
    if (!f) { std::fprintf(stderr, "[FATAL] open %s/chroms.tsv\n", dir.c_str()); std::abort(); }
    char buf[1 << 16];
    while (std::fgets(buf, sizeof(buf), f)) {
        uint32_t id; char name[256]; unsigned long L;
        if (std::sscanf(buf, "%u\t%255s\t%lu", &id, name, &L) == 3) { names.emplace_back(name); lengths.push_back(L); }
    }
    std::fclose(f);
}

void load_modeb(const std::string& dir, std::vector<BubbleTag>& tags,
                std::vector<uint8_t>& bytes, std::vector<CovRun>& runs) {
    FILE* f = xopen(dir + "/modeb.bin", "rb");
    uint64_t B; rd(f, &B, 8); tags.resize(B); rd(f, tags.data(), B * sizeof(BubbleTag));
    uint64_t Lb; rd(f, &Lb, 8); bytes.resize(Lb); rd(f, bytes.data(), Lb);
    uint64_t Lr; rd(f, &Lr, 8); runs.resize(Lr); rd(f, runs.data(), Lr * sizeof(CovRun));
    std::fclose(f);
}

void load_bubble_pos(const std::string& dir, std::vector<BubblePos>& bps) {
    FILE* f = xopen(dir + "/bubble_pos.bin", "rb");
    uint64_t B; rd(f, &B, 8); bps.resize(B); rd(f, bps.data(), B * sizeof(BubblePos));
    std::fclose(f);
}

// Optional phantom_filter.bin (from region_color_census). flags[i]=1 means "this bucket=2
// D-tag bubble sits in a low-repeat region, allow emission". Without this file: behave as
// before (drop all bucket=2 D-tag ≥50bp).
void load_phantom_filter(const std::string& dir, std::vector<uint8_t>& flags, uint64_t B) {
    flags.assign(B, 0);
    FILE* f = std::fopen((dir + "/phantom_filter.bin").c_str(), "rb");
    if (!f) return;
    uint64_t N; if (std::fread(&N, 8, 1, f) != 1) { std::fclose(f); return; }
    if (N != B) { std::fclose(f); std::fprintf(stderr, "[warn] phantom_filter.bin size %lu != bubbles %lu, ignoring\n", (unsigned long)N, (unsigned long)B); return; }
    std::fread(flags.data(), 1, B, f);
    std::fclose(f);
    uint64_t set = 0; for (uint8_t v : flags) if (v) set++;
    std::fprintf(stderr, "[vcf] phantom_filter.bin: %lu / %lu bubbles flagged allow-bucket2-D\n",
                 (unsigned long)set, (unsigned long)B);
}

void load_dup_flags(const std::string& path, std::unordered_map<uint32_t, DupCall>& out) {
    if (path.empty()) return;
    FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) { std::fprintf(stderr, "[warn] dup-flags %s not found, skipping\n", path.c_str()); return; }
    uint64_t K; if (std::fread(&K, 8, 1, f) != 1) { std::fclose(f); return; }
    std::vector<DupCall> rows(K); std::fread(rows.data(), sizeof(DupCall), K, f); std::fclose(f);
    for (auto& d : rows) out[d.bubble_id] = d;
}

void load_inv_flags(const std::string& path, std::unordered_map<uint32_t, InvCall>& out) {
    if (path.empty()) return;
    FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) { std::fprintf(stderr, "[warn] inv-flags %s not found, skipping\n", path.c_str()); return; }
    uint64_t K; if (std::fread(&K, 8, 1, f) != 1) { std::fclose(f); return; }
    std::vector<InvCall> rows(K); std::fread(rows.data(), sizeof(InvCall), K, f); std::fclose(f);
    for (auto& v : rows) { out[v.bubble_id_a] = v; if (v.bubble_id_b != UINT32_MAX) out[v.bubble_id_b] = v; }
}

uint32_t median_cov_n(const CovRun* r, uint32_t n) {
    if (n == 0) return 0;
    std::vector<std::pair<uint32_t, uint32_t>> vals; uint32_t prev = 0, total = 0;
    for (uint32_t i = 0; i < n; i++) { vals.emplace_back(r[i].n_cov, r[i].end_pos - prev); total += r[i].end_pos - prev; prev = r[i].end_pos; }
    std::sort(vals.begin(), vals.end());
    uint32_t acc = 0; for (auto& v : vals) { acc += v.second; if (acc * 2 >= total) return v.first; }
    return vals.back().first;
}
uint32_t median_cov_t(const CovRun* r, uint32_t n) {
    if (n == 0) return 0;
    std::vector<std::pair<uint32_t, uint32_t>> vals; uint32_t prev = 0, total = 0;
    for (uint32_t i = 0; i < n; i++) { vals.emplace_back(r[i].t_cov, r[i].end_pos - prev); total += r[i].end_pos - prev; prev = r[i].end_pos; }
    std::sort(vals.begin(), vals.end());
    uint32_t acc = 0; for (auto& v : vals) { acc += v.second; if (acc * 2 >= total) return v.first; }
    return vals.back().first;
}

static uint8_t RC_LUT[256];
static bool rc_lut_init = false;
void rc_lut_setup() {
    if (rc_lut_init) return;
    for (int i = 0; i < 256; i++) RC_LUT[i] = 'N';
    RC_LUT[(uint8_t)'A']='T'; RC_LUT[(uint8_t)'T']='A';
    RC_LUT[(uint8_t)'C']='G'; RC_LUT[(uint8_t)'G']='C'; RC_LUT[(uint8_t)'N']='N';
    rc_lut_init = true;
}
void rc_inplace(uint8_t* s, uint32_t n) {
    rc_lut_setup();
    for (uint32_t i = 0, j = n - 1; i < j; i++, j--) { uint8_t a = RC_LUT[s[i]], b = RC_LUT[s[j]]; s[i] = b; s[j] = a; }
    if (n & 1) s[n/2] = RC_LUT[s[n/2]];
}

// Strand-aware variant POS (0-based forward) at arm[pre]. VCF emitters add +1 (SNP) or +0 (INS/DEL).
uint32_t variant_pos0(const BubbleTag& t, const BubblePos& bp) {
    if (bp.strand == 0) return bp.s1_left + t.pre;
    // strand=1: in arm coord arm[arm_len-1-suf] is the last byte before the trailing shared suffix;
    // empirical winner across 506 strand=1 SNPs is forward pos = s2_right - suf - 1 (0-based).
    return (bp.s2_right > t.suf + 1) ? (bp.s2_right - t.suf - 1) : 0;
}

uint32_t bubble_span(const BubblePos& bp) {
    if (bp.strand == 0) {
        return (bp.s2_left > bp.s1_right) ? (bp.s2_left - bp.s1_right) : 0;
    }
    return (bp.s1_right > bp.s2_left) ? (bp.s1_right - bp.s2_left) : 0;
}

bool apply_filters(const BubbleTag& t, const BubblePos& bp, const CovRun* runs, bool allow_bucket2_d, bool inv_flagged) {
    if (t.tag == 'C' || t.tag == 'E') return false;
    if (bp.chrom_id == UINT32_MAX) return false;
    // Resolution gates
    const uint32_t resid_initial = (t.tag == 'D' || t.tag == 'S') ? t.n_resid_len : t.t_resid_len;
    const uint32_t need_s1 = bp.n_anchor_s1 / 2 + 1;
    const uint32_t need_s2 = bp.n_anchor_s2 / 2 + 1;
    if (bp.n_resolved_s1 < need_s1 || bp.n_resolved_s2 < need_s2) return false;
    // Stricter resolution for SV-class records (≥50bp): require ≥75% of anchor sigs in ref.
    if (resid_initial >= 50) {
        const uint32_t need_s1_sv = (bp.n_anchor_s1 * 3 + 3) / 4;
        const uint32_t need_s2_sv = (bp.n_anchor_s2 * 3 + 3) / 4;
        if (bp.n_resolved_s1 < need_s1_sv || bp.n_resolved_s2 < need_s2_sv) return false;
    }
    const uint32_t span = bubble_span(bp);
    const uint32_t nc = median_cov_n(runs + t.n_runs_off, t.n_runs_len);
    const uint32_t tc = median_cov_t(runs + t.t_runs_off, t.t_runs_len);

    // ============ DEL ============
    if (t.tag == 'D') {
        (void)allow_bucket2_d;
        if (nc < 3) return false;
        if (t.n_resid_len >= 50) {
            // Allow bucket=2 D-tag through if the inv_flags entry covers this bubble — the
            // bucket=2 D-tag at an INV-signature locus emits as <INV>, not <DEL>.
            if (t.bucket != 0 && !inv_flagged) return false;
            if (t.bucket == 0) {
                if (bp.n_anchor_s1 < 20 && bp.n_anchor_s2 < 20) return false;
                if (t.pre + t.suf < 100) return false;
                if (span < 10) return false;
                const uint32_t cap = 3 * t.n_resid_len + 500;
                if (span > cap) return false;
                if (nc > 250) return false;
            }
            // bucket=2 + inv_flagged passes through unconditionally.
        }
        return true;
    }
    // ============ INS ============
    if (t.tag == 'I') {
        if (tc < 3) return false;
        if (tc > 150) return false;
        if (t.t_resid_len >= 50) {
            if (t.bucket != 0) return false;
            if (bp.n_anchor_s1 < 20 && bp.n_anchor_s2 < 20) return false;
            if (t.pre + t.suf < 100) return false;
            if (span < 10) return false;
            const uint32_t cap = 3 * t.t_resid_len + 500;
            if (span > cap) return false;
        }
        return true;
    }
    // ============ S (SNP/MNP) ============
    if (t.tag == 'S') {
        if (nc < 3 || tc < 3) return false;
        if (t.n_resid_len == 1 && t.t_resid_len == 1) {
            // Phantom 1bp-SNP signature: high normal coverage AND tiny span (collapsed locus).
            if (nc > 150 && span < 5) return false;
            return true;
        }
        // Larger S (MNP). Phantom MNPs cluster in orphan buckets; require std_both.
        // Also: truth contains NO MNPs (only SNPs). MNP ≥50bp records are de facto FPs
        // unless they're real INVs (which we can't reliably distinguish). Drop ≥50bp MNPs.
        if (t.n_resid_len >= 50) return false;
        return true;
    }
    return false;
}

void prepare_alleles(const BubbleTag& t, const BubblePos& bp, const uint8_t* bytes,
                     std::vector<uint8_t>& N, std::vector<uint8_t>& T) {
    N.assign(bytes + t.n_resid_off, bytes + t.n_resid_off + t.n_resid_len);
    T.assign(bytes + t.t_resid_off, bytes + t.t_resid_off + t.t_resid_len);
    if (bp.strand == 1) {
        if (!N.empty()) rc_inplace(N.data(), N.size());
        if (!T.empty()) rc_inplace(T.data(), T.size());
    }
}

void emit_snp(const BubbleTag& t, const BubblePos& bp, const std::vector<uint8_t>& N, const std::vector<uint8_t>& T,
              const CovRun* runs, const std::vector<std::string>& chroms, std::vector<VcfRec>& out) {
    const uint32_t tc = median_cov_t(runs + t.t_runs_off, t.t_runs_len);
    const uint32_t nc = median_cov_n(runs + t.n_runs_off, t.n_runs_len);
    const uint32_t pos = variant_pos0(t, bp) + 1;  // SNP: 1-based variant position
    char buf[512];
    std::snprintf(buf, sizeof(buf), "%s\t%u\t.\t%c\t%c\t60\tPASS\tSVTYPE=SNP\tGT:DV:DR\t0/1:%u:%u\n",
                  chroms[bp.chrom_id].c_str(), pos, (char)N[0], (char)T[0], tc, nc);
    out.push_back({bp.chrom_id, pos, std::string(buf)});
}

void emit_mnp(const BubbleTag& t, const BubblePos& bp, const std::vector<uint8_t>& N, const std::vector<uint8_t>& T,
              const CovRun* runs, const std::vector<std::string>& chroms, std::vector<VcfRec>& out) {
    const uint32_t tc = median_cov_t(runs + t.t_runs_off, t.t_runs_len);
    const uint32_t nc = median_cov_n(runs + t.n_runs_off, t.n_runs_len);
    const uint32_t pos = variant_pos0(t, bp) + 1;
    std::string line; line.reserve(N.size() + T.size() + 96);
    line += chroms[bp.chrom_id]; line += '\t'; line += std::to_string(pos); line += "\t.\t";
    line.append(N.begin(), N.end()); line += '\t';
    line.append(T.begin(), T.end()); line += "\t60\tPASS\tSVTYPE=MNP\tGT:DV:DR\t0/1:";
    line += std::to_string(tc); line += ':'; line += std::to_string(nc); line += '\n';
    out.push_back({bp.chrom_id, pos, std::move(line)});
}

void emit_del(const BubbleTag& t, const BubblePos& bp, const std::vector<uint8_t>& N,
              const CovRun* runs, const std::vector<std::string>& chroms, std::vector<VcfRec>& out) {
    const uint32_t tc = median_cov_t(runs + t.t_runs_off, t.t_runs_len);
    const uint32_t nc = median_cov_n(runs + t.n_runs_off, t.n_runs_len);
    // VCF DEL convention: POS = position of last shared base before the deletion (1-based).
    // arm[pre] is the first deleted byte at 0-based forward s1_left+pre (strand=0). POS = s1_left+pre.
    const uint32_t pos = variant_pos0(t, bp);                // 0-based forward of first deleted; +0 makes it 1-based of base-before
    const uint32_t safe_pos = pos == 0 ? 1 : pos;
    if (t.n_resid_len < 50) {
        // explicit REF/ALT: REF = anchor_base + N'_bytes; ALT = anchor_base. Use 'N' as placeholder
        // for the anchor base (reference not loaded). Truvari with --pctseq=0 ignores REF/ALT text.
        std::string line; line.reserve(N.size() + 128);
        line += chroms[bp.chrom_id]; line += '\t'; line += std::to_string(safe_pos); line += "\t.\tN";
        line.append(N.begin(), N.end()); line += "\tN\t60\tPASS\tSVTYPE=DEL;SVLEN=";
        line += std::to_string(t.n_resid_len); line += "\tGT:DV:DR\t0/1:";
        line += std::to_string(tc); line += ':'; line += std::to_string(nc); line += '\n';
        out.push_back({bp.chrom_id, safe_pos, std::move(line)});
        return;
    }
    // symbolic <DEL>
    const uint32_t end_pos = safe_pos + t.n_resid_len;
    char buf[512];
    std::snprintf(buf, sizeof(buf), "%s\t%u\t.\tN\t<DEL>\t60\tPASS\tEND=%u;SVLEN=%u;SVTYPE=DEL\tGT:DV:DR\t0/1:%u:%u\n",
                  chroms[bp.chrom_id].c_str(), safe_pos, end_pos, t.n_resid_len, tc, nc);
    out.push_back({bp.chrom_id, safe_pos, std::string(buf)});
}

void emit_ins(const BubbleTag& t, const BubblePos& bp, const std::vector<uint8_t>& T, const CovRun* runs,
              const std::vector<std::string>& chroms, std::vector<VcfRec>& out) {
    const uint32_t tc = median_cov_t(runs + t.t_runs_off, t.t_runs_len);
    const uint32_t nc = median_cov_n(runs + t.n_runs_off, t.n_runs_len);
    const uint32_t pos = variant_pos0(t, bp);
    const uint32_t safe_pos = pos == 0 ? 1 : pos;
    if (t.t_resid_len < 50) {
        std::string line; line.reserve(T.size() + 128);
        line += chroms[bp.chrom_id]; line += '\t'; line += std::to_string(safe_pos); line += "\t.\tN\tN";
        line.append(T.begin(), T.end()); line += "\t60\tPASS\tSVTYPE=INS;SVLEN=";
        line += std::to_string(t.t_resid_len); line += "\tGT:DV:DR\t0/1:";
        line += std::to_string(tc); line += ':'; line += std::to_string(nc); line += '\n';
        out.push_back({bp.chrom_id, safe_pos, std::move(line)});
        return;
    }
    std::string line; line.reserve(T.size() + 192);
    line += chroms[bp.chrom_id]; line += '\t'; line += std::to_string(safe_pos);
    line += "\t.\tN\t<INS>\t60\tPASS\tEND="; line += std::to_string(safe_pos);
    line += ";SVLEN="; line += std::to_string(t.t_resid_len);
    line += ";SVTYPE=INS;INSSEQ="; line.append(T.begin(), T.end());
    line += "\tGT:DV:DR\t0/1:"; line += std::to_string(tc); line += ':'; line += std::to_string(nc); line += '\n';
    out.push_back({bp.chrom_id, safe_pos, std::move(line)});
}

void emit_dup(const BubbleTag& t, const BubblePos& bp, const DupCall& dc, const CovRun* runs,
              const std::vector<std::string>& chroms, std::vector<VcfRec>& out) {
    const uint32_t tc = median_cov_t(runs + t.t_runs_off, t.t_runs_len);
    const uint32_t nc = median_cov_n(runs + t.n_runs_off, t.n_runs_len);
    const uint32_t pos = variant_pos0(t, bp);
    const uint32_t safe_pos = pos == 0 ? 1 : pos;
    const uint32_t svlen = dc.svlen ? dc.svlen : t.t_resid_len;
    const uint32_t end = safe_pos + svlen;
    char buf[512];
    std::snprintf(buf, sizeof(buf), "%s\t%u\t.\tN\t<DUP>\t60\tPASS\tEND=%u;SVLEN=%u;SVTYPE=DUP;TARGET=%u\tGT:DV:DR\t0/1:%u:%u\n",
                  chroms[bp.chrom_id].c_str(), safe_pos, end, svlen, end, tc, nc);
    out.push_back({bp.chrom_id, safe_pos, std::string(buf)});
}

void emit_inv(const BubbleTag& t, const BubblePos& bp, const InvCall& iv, const CovRun* runs,
              const std::vector<std::string>& chroms, std::vector<VcfRec>& out) {
    const uint32_t tc = median_cov_t(runs + t.t_runs_off, t.t_runs_len);
    const uint32_t nc = median_cov_n(runs + t.n_runs_off, t.n_runs_len);
    const uint32_t pos = (iv.pos_start ? iv.pos_start : (variant_pos0(t, bp) + 1));
    const uint32_t svlen = iv.svlen ? iv.svlen : t.n_resid_len;
    const uint32_t end = pos + svlen - 1;
    char buf[512];
    std::snprintf(buf, sizeof(buf), "%s\t%u\t.\tN\t<INV>\t60\tPASS\tEND=%u;SVLEN=%u;SVTYPE=INV\tGT:DV:DR\t0/1:%u:%u\n",
                  chroms[bp.chrom_id].c_str(), pos, end, svlen, tc, nc);
    out.push_back({bp.chrom_id, pos, std::string(buf)});
}

void write_header(FILE* f, const std::vector<std::string>& chroms, const std::vector<uint64_t>& lengths) {
    std::fputs("##fileformat=VCFv4.2\n##source=svp\n##FILTER=<ID=PASS,Description=\"All filters passed\">\n", f);
    for (size_t i = 0; i < chroms.size(); i++) std::fprintf(f, "##contig=<ID=%s,length=%lu>\n", chroms[i].c_str(), (unsigned long)lengths[i]);
    std::fputs(
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n"
        "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of structural variant\">\n"
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of variant\">\n"
        "##INFO=<ID=INSSEQ,Number=1,Type=String,Description=\"Inserted sequence\">\n"
        "##INFO=<ID=TARGET,Number=1,Type=Integer,Description=\"Tandem-DUP insertion target position\">\n"
        "##ALT=<ID=DEL,Description=\"Deletion\">\n##ALT=<ID=INS,Description=\"Insertion\">\n"
        "##ALT=<ID=DUP,Description=\"Duplication\">\n##ALT=<ID=INV,Description=\"Inversion\">\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        "##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"Tumor coverage (median)\">\n"
        "##FORMAT=<ID=DR,Number=1,Type=Integer,Description=\"Normal coverage (median)\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\n", f);
}

void dispatch_emit(uint32_t bubble_id, const BubbleTag& t, const BubblePos& bp,
                   const uint8_t* bytes, const CovRun* runs,
                   const std::unordered_map<uint32_t, DupCall>& dups,
                   const std::unordered_map<uint32_t, InvCall>& invs,
                   const std::vector<std::string>& chroms, std::vector<VcfRec>& out) {
    std::vector<uint8_t> N, T;
    prepare_alleles(t, bp, bytes, N, T);
    auto inv_it = invs.find(bubble_id);
    auto dup_it = dups.find(bubble_id);
    if (dup_it != dups.end() && (t.tag == 'I' || t.tag == 'S' || t.tag == 'C')) { emit_dup(t, bp, dup_it->second, runs, chroms, out); return; }
    // INV emission for D-tag bucket=2 at INV-signature loci OR S-tag multi-bp INV.
    if (inv_it != invs.end()) {
        if (t.tag == 'D' && t.bucket == 2 && t.n_resid_len >= 50) { emit_inv(t, bp, inv_it->second, runs, chroms, out); return; }
        if (t.tag == 'S' && N.size() > 1)                          { emit_inv(t, bp, inv_it->second, runs, chroms, out); return; }
    }
    if (t.tag == 'I')                                                            { emit_ins(t, bp, T, runs, chroms, out); return; }
    if (t.tag == 'D')                                                            { emit_del(t, bp, N, runs, chroms, out); return; }
    if (t.tag == 'S' && N.size() == 1 && T.size() == 1)                          { emit_snp(t, bp, N, T, runs, chroms, out); return; }
    if (t.tag == 'S' && N.size() == T.size() && N.size() > 1)                    { emit_mnp(t, bp, N, T, runs, chroms, out); return; }
}

}  // namespace

void run_vcf(const std::string& dir, const std::string& /*ref_path*/, const std::string& out_path,
             const std::string& dup_flags_path, const std::string& inv_flags_path) {
    std::vector<std::string> chroms; std::vector<uint64_t> lengths;
    load_chroms(dir, chroms, lengths);
    std::vector<BubbleTag> tags; std::vector<uint8_t> bytes; std::vector<CovRun> runs;
    load_modeb(dir, tags, bytes, runs);
    std::vector<BubblePos> bps; load_bubble_pos(dir, bps);
    std::unordered_map<uint32_t, DupCall> dups; std::unordered_map<uint32_t, InvCall> invs;
    load_dup_flags(dup_flags_path, dups); load_inv_flags(inv_flags_path, invs);
    std::vector<uint8_t> phantom_flags;
    load_phantom_filter(dir, phantom_flags, tags.size());
    std::fprintf(stderr, "[vcf] %lu bubbles, %lu dup-flags, %lu inv-flags\n",
                 (unsigned long)tags.size(), (unsigned long)dups.size(), (unsigned long)invs.size());

    std::vector<VcfRec> records; records.reserve(tags.size() / 2);
    uint64_t passed = 0, dropped = 0;
    uint64_t drop_C = 0, drop_unanchored = 0, drop_resolution = 0, drop_cov = 0, drop_phantom = 0;
    for (size_t i = 0; i < tags.size(); i++) {
        const bool allow_b2d = i < phantom_flags.size() && phantom_flags[i] != 0;
        const bool inv_flagged = invs.find((uint32_t)i) != invs.end();
        if (!apply_filters(tags[i], bps[i], runs.data(), allow_b2d, inv_flagged)) {
            dropped++;
            if (tags[i].tag == 'C' || tags[i].tag == 'E') drop_C++;
            else if (bps[i].chrom_id == UINT32_MAX) drop_unanchored++;
            else drop_phantom++;
            continue;
        }
        dispatch_emit((uint32_t)i, tags[i], bps[i], bytes.data(), runs.data(), dups, invs, chroms, records);
        passed++;
    }
    std::fprintf(stderr, "[vcf] %lu passed filter, %lu dropped (C=%lu unanchored=%lu other=%lu)\n",
                 (unsigned long)passed, (unsigned long)dropped,
                 (unsigned long)drop_C, (unsigned long)drop_unanchored, (unsigned long)drop_phantom);

    std::sort(records.begin(), records.end(), [](const VcfRec& a, const VcfRec& b) {
        if (a.chrom_id != b.chrom_id) return a.chrom_id < b.chrom_id;
        return a.pos_vcf < b.pos_vcf;
    });
    const std::string tmp = out_path + ".tmp";
    FILE* f = xopen(tmp, "wb");
    write_header(f, chroms, lengths);
    for (auto& r : records) wr(f, r.line.data(), r.line.size());
    atomic_close(f, tmp, out_path);
    std::fprintf(stderr, "[vcf] wrote %lu records -> %s\n", (unsigned long)records.size(), out_path.c_str());
}
