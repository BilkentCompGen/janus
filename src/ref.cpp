#include <algorithm>
#include <string>
#include <vector>
#include <htslib/vcf.h>
#include "ref.h"

std::vector<VariantRegion> parseVCF(const std::string& path) {
    std::vector<VariantRegion> res;
    auto* fp = bcf_open(path.c_str(), "r"); if (!fp) return res;
    auto* hdr = bcf_hdr_read(fp);
    auto* rec = bcf_init();
    while (bcf_read(fp, hdr, rec) == 0) {
        bcf_unpack(rec, BCF_UN_INFO);
        char *sv = nullptr; int n = 0;
        std::string type = "SNP";
        if (bcf_get_info_string(hdr, rec, "SVTYPE", &sv, &n) >= 0) { type = sv; free(sv); }
        int64_t end = rec->pos + 1;
        int64_t *end_ptr = nullptr; int n2 = 0;
        if (bcf_get_info_int64(hdr, rec, "END", &end_ptr, &n2) >= 0) { end = end_ptr[0]; free(end_ptr); }
        std::string chr = bcf_hdr_id2name(hdr, rec->rid);
        if (chr.find("chr") == 0) chr = chr.substr(3);
        res.push_back({chr, type, (uint64_t)rec->pos + 1, (uint64_t)end});
    }
    bcf_destroy(rec); bcf_hdr_destroy(hdr); bcf_close(fp);
    return res;
}

void exportAlignedVCF(const std::string& path, std::vector<VariantCall>& vars) {
    if (vars.empty()) return;
    std::sort(vars.begin(), vars.end(), [](const auto& a, const auto& b) {
        if (a.pos != b.pos) return a.pos < b.pos;
        if (a.svtype != b.svtype) return a.svtype < b.svtype;
        return a.ref_allele < b.ref_allele;
    });
    auto* fp = bcf_open(path.c_str(), "w");
    auto* hdr = bcf_hdr_init("w");
    bcf_hdr_append(hdr, "##fileformat=VCFv4.2");
    bcf_hdr_append(hdr, "##contig=<ID=19,length=58617616>");
    bcf_hdr_append(hdr, "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type\">");
    bcf_hdr_append(hdr, "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Len\">");
    bcf_hdr_append(hdr, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">");
    bcf_hdr_append(hdr, "##ALT=<ID=DEL,Description=\"Deletion\">");
    bcf_hdr_append(hdr, "##ALT=<ID=INS,Description=\"Insertion\">");
    if (bcf_hdr_write(fp, hdr) < 0) { bcf_hdr_destroy(hdr); bcf_close(fp); return; }
    auto* rec = bcf_init();
    for (size_t i = 0; i < vars.size(); ++i) {
        if (i > 0 && vars[i].pos == vars[i-1].pos && vars[i].svtype == vars[i-1].svtype && vars[i].ref_allele == vars[i-1].ref_allele && vars[i].alt_allele == vars[i-1].alt_allele) continue;
        const auto& v = vars[i];
        bcf_clear(rec);
        rec->rid = bcf_hdr_name2id(hdr, v.chrom.c_str());
        rec->pos = v.pos - 1;
        bcf_update_alleles_str(hdr, rec, (v.ref_allele + "," + v.alt_allele).c_str());
        bcf_update_info_string(hdr, rec, "SVTYPE", v.svtype.c_str());
        if (v.svtype != "SNP") {
            int32_t slen = (int32_t)v.svlen, end = (int32_t)(v.pos + std::abs(v.svlen));
            bcf_update_info_int32(hdr, rec, "SVLEN", &slen, 1);
            bcf_update_info_int32(hdr, rec, "END", &end, 1);
        }
        if (bcf_write(fp, hdr, rec) < 0) break;
    }
    bcf_destroy(rec); bcf_hdr_destroy(hdr); bcf_close(fp);
}
