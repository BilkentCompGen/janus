#include <algorithm>
#include <string>
#include <unordered_set>
#include <vector>
#include <htslib/vcf.h>
#include "ref.h"

std::vector<DeletionRegion> parseVCFDeletions(const std::string& path) {
    std::vector<DeletionRegion> dels;
    auto* fp = bcf_open(path.c_str(), "r"); if (!fp) return dels;
    auto* hdr = bcf_hdr_read(fp);
    auto* rec = bcf_init();

    while (bcf_read(fp, hdr, rec) == 0) {
        bcf_unpack(rec, BCF_UN_INFO);
        char *sv = nullptr; int n = 0;
        if (bcf_get_info_string(hdr, rec, "SVTYPE", &sv, &n) < 0) continue;
        if (std::string(sv) == "DEL") {
            int64_t *end = nullptr; int n2 = 0;
            if (bcf_get_info_int64(hdr, rec, "END", &end, &n2) >= 0) {
                std::string chr = bcf_hdr_id2name(hdr, rec->rid);
                if (chr.find("chr") == 0) chr = chr.substr(3);
                dels.push_back({chr, (uint64_t)rec->pos + 1, (uint64_t)end[0]});
                free(end);
            }
        }
        free(sv);
    }
    bcf_destroy(rec); bcf_hdr_destroy(hdr); bcf_close(fp);
    return dels;
}

std::unordered_set<node_hash> findDeletionNodeHashes(const std::vector<DeletionRegion>& dels, const LCPGraph& g, uint32_t rsid) {
    std::unordered_set<node_hash> res;
    for (const auto& [h, n] : g.nodes) {
        for (const auto& o : n.occurrences) {
            if (o.seq_id != rsid) continue;
            for (const auto& d : dels) {
                if (o.start < d.end && o.end > d.start) { res.insert(h); break; }
            }
            if (res.count(h)) break;
        }
    }
    return res;
}

void exportAlignedVCF(const std::string& path, std::vector<VariantCall>& vars) {
    if (vars.empty()) return;
    std::sort(vars.begin(), vars.end(), [](const auto& a, const auto& b) { return a.pos < b.pos; });

    auto* fp = bcf_open(path.c_str(), "w");
    auto* hdr = bcf_hdr_init("w");
    bcf_hdr_append(hdr, "##fileformat=VCFv4.2");
    bcf_hdr_append(hdr, "##contig=<ID=19,length=58617616>");
    bcf_hdr_append(hdr, "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"SV Type\">");
    bcf_hdr_append(hdr, "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"SV Length\">");
    bcf_hdr_append(hdr, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End Position\">");
    bcf_hdr_write(fp, hdr);

    auto* rec = bcf_init();
    for (const auto& v : vars) {
        bcf_clear(rec);
        rec->rid = bcf_hdr_name2id(hdr, v.chrom.c_str());
        rec->pos = v.pos - 1;
        std::string alleles = v.ref_allele + "," + v.alt_allele;
        bcf_update_alleles_str(hdr, rec, alleles.c_str());
        bcf_update_info_string(hdr, rec, "SVTYPE", v.svtype.c_str());
        int32_t slen = (int32_t)v.svlen;
        bcf_update_info_int32(hdr, rec, "SVLEN", &slen, 1);
        int32_t end_pos = (int32_t)(v.pos + std::abs(v.svlen));
        bcf_update_info_int32(hdr, rec, "END", &end_pos, 1);
        bcf_write(fp, hdr, rec);
    }
    bcf_destroy(rec); bcf_hdr_destroy(hdr); bcf_close(fp);
}
