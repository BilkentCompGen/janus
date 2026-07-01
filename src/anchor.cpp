// svp anchor <dir> <ref.fasta>  (run_anchor): coordinate-reporting stage AFTER call.
// The graph/refine/call stages are reference-free; this final stage maps the coordinate-free
// somatic calls onto a REFERENCE genome so they can be written as a VCF. ref.fasta is the only
// reference allowed here (a coordinate oracle) — NEVER the tumor/normal assemblies.
//
// How: a somatic alt is all colour-2 (tumor-only, absent from the reference). Its graph neighbours
// with colour in {1,3} (normal-only / shared) ARE present in the reference, so their core_sig is
// locatable in ref.fasta. We map those anchor sigs to forward-genome positions (streaming ref
// fwd+rc, needed sigs only, unique positions) and bracket the variant between the min/max anchor.
// SVLEN: caller's spliced-Δ size for SNP/INS/DEL; the anchored span for BND/DUP/INV (caller leaves
// those unsized). Reads dir/refined.bin (else graph.bin) + dir/calls.{tsv,altsigs.bin}; writes dir/calls.vcf.
#define XXH_INLINE_ALL
#include <xxhash.h>
#include "svp.h"
#include "graph.h"
#include "unordered_dense.h"
#include <fcntl.h>
#include <unistd.h>
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

namespace {
constexpr int LEVEL = 5, MARGIN = 10;   // MARGIN must match the graph build (blacklist merge operates on margin-kept cores)
uint8_t COMP[256];
void initc(){ for(int i=0;i<256;i++)COMP[i]='N'; COMP['A']='T';COMP['T']='A';COMP['C']='G';COMP['G']='C'; }
[[noreturn]] void die(const char* m){ std::fprintf(stderr,"[anchor] FATAL %s\n",m); std::abort(); }
bool exists(const std::string& p){ FILE* f=std::fopen(p.c_str(),"rb"); if(f){std::fclose(f);return true;} return false; }

struct G{
  uint64_t K=0,E=0;
  std::vector<uint8_t> color; std::vector<uint32_t> ncov,tcov; std::vector<uint16_t> seqlen;
  std::vector<uint8_t> payload; std::vector<uint64_t> off;
  std::vector<uint32_t> out_off,out_nbr,in_off,in_nbr;
};
void load(const char* gp, G& g){
  FILE* f=std::fopen(gp,"rb"); if(!f)die("open graph");
  std::fread(&g.K,8,1,f); const uint64_t K=g.K;
  g.color.resize(K); std::fread(g.color.data(),1,K,f);
  g.ncov.resize(K); g.tcov.resize(K); std::fread(g.ncov.data(),4,K,f); std::fread(g.tcov.data(),4,K,f);
  std::fread(&g.E,8,1,f); std::vector<std::pair<uint32_t,uint32_t>> ed(g.E); std::fread(ed.data(),8,g.E,f);
  uint64_t total; std::fread(&total,8,1,f); g.seqlen.resize(K); std::fread(g.seqlen.data(),2,K,f);
  g.payload.resize(total); std::fread(g.payload.data(),1,total,f); std::fclose(f);
  g.off.assign(K+1,0); for(uint64_t i=0;i<K;i++) g.off[i+1]=g.off[i]+g.seqlen[i];
  const uint64_t E=g.E; g.out_off.assign(K+1,0); g.in_off.assign(K+1,0);
  for(auto&e:ed){g.out_off[e.first+1]++;g.in_off[e.second+1]++;}
  for(uint64_t i=0;i<K;i++){g.out_off[i+1]+=g.out_off[i];g.in_off[i+1]+=g.in_off[i];}
  g.out_nbr.resize(E); g.in_nbr.resize(E);
  std::vector<uint32_t> oc(g.out_off.begin(),g.out_off.end()-1), ic(g.in_off.begin(),g.in_off.end()-1);
  for(auto&e:ed){g.out_nbr[oc[e.first]++]=e.second; g.in_nbr[ic[e.second]++]=e.first;}
  std::fprintf(stderr,"[anchor] graph K=%lu E=%lu payload=%.1fGB\n",(unsigned long)K,(unsigned long)E,total/1e9);
}
inline const char* P(const G&g,uint32_t x){ return (const char*)g.payload.data()+g.off[x]; }
inline uint64_t csig(const G&g,uint32_t x){ return XXH3_64bits(P(g,x),g.seqlen[x]); }

// faidx single-base fetch for the VCF REF column + contig lengths/order for the header.
struct FaEntry{uint64_t len,offset;uint32_t lb,lw;}; ankerl::unordered_dense::map<std::string,FaEntry> IX;
std::vector<std::string> FAORDER; int FD=-1;
void load_faidx(const std::string& fa){ FD=open(fa.c_str(),O_RDONLY); if(FD<0)die("open ref.fasta");
  FILE* fi=std::fopen((fa+".fai").c_str(),"r"); if(!fi)die("open .fai (run: samtools faidx ref.fasta)"); char nm[512]; uint64_t l,o; uint32_t a,b;
  while(std::fscanf(fi,"%511s\t%lu\t%lu\t%u\t%u\n",nm,&l,&o,&a,&b)==5){ IX[nm]={l,o,a,b}; FAORDER.push_back(nm); } std::fclose(fi); }
char fetch1(const std::string& c,long s0){ auto it=IX.find(c); if(it==IX.end()||s0<0)return 'N';
  const FaEntry&e=it->second; if((uint64_t)s0>=e.len)return 'N';
  uint64_t fs=e.offset+(s0/e.lb)*e.lw+(s0%e.lb); char ch='N'; if(pread(FD,&ch,1,(off_t)fs)!=1)return 'N';
  return (ch>='a'&&ch<='z')?(char)(ch-32):ch; }
std::string fetch_range(const std::string& c,long s0,long want){ auto it=IX.find(c); if(it==IX.end()||want<=0||s0<0)return"";
  const FaEntry&e=it->second; if(s0>=(long)e.len)return""; if((uint64_t)(s0+want)>e.len)want=e.len-s0;
  uint64_t fs=e.offset+(s0/e.lb)*e.lw+(s0%e.lb),nl=e.lw-e.lb,fb=(uint64_t)want+((uint64_t)want/e.lb+2)*nl;
  std::string raw(fb,'\0'); ssize_t g=pread(FD,&raw[0],fb,(off_t)fs); if(g<=0)return""; raw.resize((size_t)g);
  std::string out; out.reserve(want); for(char ch:raw){ if(ch=='\n'||ch=='\r')continue; out.push_back(ch>='a'&&ch<='z'?(char)(ch-32):ch); if((long)out.size()==want)break; } return out; }

// caller `type` column -> VCF SVTYPE class.  0 SNP 1 INS 2 DEL 3 BND 4 BND_CNgain 5 TIP
uint8_t type_code(const char* t){
  if(!std::strcmp(t,"SNP"))return 0;
  if(!std::strcmp(t,"INS"))return 1;
  if(!std::strcmp(t,"DEL"))return 2;
  if(!std::strcmp(t,"BND_CNgain"))return 4;
  if(!std::strcmp(t,"TIP"))return 5;
  return 3;  // BND
}
const char* SVT[6]={"SNP","INS","DEL","BND","BND","BND"};
const char* SVALT[6]={"<SNP>","<INS>","<DEL>","<BND>","<BND>","<BND>"};

struct Call{ uint64_t id; uint8_t tc; long size,jlen,sup; float vaf; std::vector<uint32_t> anchors; };
struct Rec{ uint32_t cidx; uint64_t pos,end; long svlen,sup; uint8_t tc,local; uint64_t id; uint32_t nanc; char rb; };
struct Hit{ uint16_t c, sl; uint64_t s; };   // located anchor: contig, core seqlen, fwd start
bool hit_lt(const Hit&a,const Hit&b){ return a.c!=b.c ? a.c<b.c : a.s<b.s; }
// dedup class: SNP/INS/DEL distinct; BND(3) & BND_CNgain(4) together; TIP(5) own.
int cls(uint8_t tc){ return tc<=2 ? tc : (tc==5 ? 5 : 3); }
bool rec_locus_lt(const Rec&a,const Rec&b){ if(a.cidx!=b.cidx)return a.cidx<b.cidx; int ca=cls(a.tc),cb=cls(b.tc); if(ca!=cb)return ca<cb; return a.pos<b.pos; }
bool rec_pos_lt(const Rec&a,const Rec&b){ return a.cidx!=b.cidx ? a.cidx<b.cidx : a.pos<b.pos; }

// record sig->packed-pos; a second DISTINCT position marks the sig ambiguous (repeat) and drops it.
void rec_pos(uint64_t sig,uint64_t v,ankerl::unordered_dense::map<uint64_t,uint64_t>&rpos,ankerl::unordered_dense::set<uint64_t>&ambig){
  if(ambig.count(sig))return;
  auto it=rpos.find(sig);
  if(it==rpos.end())rpos[sig]=v; else if(it->second!=v){ambig.insert(sig);rpos.erase(it);}
}
// record needed sigs of ONE strand's cores. rc=true maps a core [st,en) back to its forward-genome
// start (len-en). When a blacklist is active, the reference is parsed through the IDENTICAL
// (A+X)->(X+B) merge the graph used (margin-kept cores -> blacklist_merge -> core_sig per span), so the
// merged reference cores match the graph's merged nodes.
void index_pass(const char* data,uint64_t len,bool rc,uint16_t cur,const ankerl::unordered_dense::set<uint64_t>&needed,
                ankerl::unordered_dense::map<uint64_t,uint64_t>&rpos,ankerl::unordered_dense::set<uint64_t>&ambig){
  int n; const uint64_t *st,*en,*sg;
  extract(data,(int)len,0,LEVEL,&n,&st,&en,&sg);
  if(blacklist_on()){
    std::vector<int> keep; std::vector<std::pair<uint32_t,uint32_t>> eff;
    for(int i=0;i<n;i++) if(st[i]>=(uint64_t)MARGIN && en[i]+MARGIN<=len) keep.push_back(i);
    blacklist_merge(st,en,sg,keep,65535,eff);
    for(const auto&e:eff){ uint64_t sig=core_sig(data+e.first,(int)(e.second-e.first));
      if(needed.count(sig)){ uint64_t pos=rc?(len-e.second):e.first; rec_pos(sig,((uint64_t)cur<<40)|pos,rpos,ambig); } }
    return;
  }
  for(int i=0;i<n;i++) if(needed.count(sg[i])){ uint64_t pos=rc?(len-en[i]):st[i]; rec_pos(sg[i],((uint64_t)cur<<40)|pos,rpos,ambig); }
}
// extract one contig's cores (fwd + rc->fwd) and record positions of the needed anchor sigs. Clears seq.
void index_seq(std::string&seq,uint16_t cur,const ankerl::unordered_dense::set<uint64_t>&needed,
               ankerl::unordered_dense::map<uint64_t,uint64_t>&rpos,ankerl::unordered_dense::set<uint64_t>&ambig){
  if(seq.empty())return;
  index_pass(seq.data(),seq.size(),false,cur,needed,rpos,ambig);
  const uint64_t SL=seq.size(); std::string r(SL,'N'); for(uint64_t j=0;j<SL;j++)r[j]=(char)COMP[(uint8_t)seq[SL-1-j]];
  index_pass(r.data(),SL,true,cur,needed,rpos,ambig);
  seq.clear();
}

// parse one calls.tsv line -> id, type-code, size, vaf, junction-len (cols 0,3,4,9,10). false on header/blank.
bool parse_call(char* ln, uint64_t& id, uint8_t& tc, long& size, float& vaf, long& jlen, long& sup){
  if(ln[0]=='c'||ln[0]=='\0'||ln[0]=='\n') return false;
  char* f[11]; int nf=0; char* s=ln; for(char* p=ln; *p&&nf<11; p++){ if(*p=='\t'){*p='\0'; f[nf++]=s; s=p+1;} }
  if(nf<11) f[nf++]=s;
  if(nf<11) return false;
  id=std::strtoull(f[0],nullptr,10); tc=type_code(f[3]); size=std::atol(f[4]); sup=std::atol(f[5]);
  vaf=(float)std::atof(f[9]); jlen=std::atol(f[10]); return true;
}
}  // namespace

void run_anchor(const std::string& dir, const std::string& ref){
  initc(); load_faidx(ref);
  if(const char* blp=std::getenv("SVP_BLACKLIST")) load_blacklist(blp);   // same merge as the graph, applied to the reference
  const std::string gp = exists(dir+"/refined.bin") ? dir+"/refined.bin" : dir+"/graph.bin";
  G g; load(gp.c_str(),g); const uint64_t K=g.K;

  ankerl::unordered_dense::map<uint64_t,uint32_t> sig2node; sig2node.reserve(K);
  for(uint64_t i=0;i<K;i++) sig2node[csig(g,(uint32_t)i)]=(uint32_t)i;
  std::fprintf(stderr,"[anchor] sig2node built (%zu)\n",sig2node.size());

  // read calls.altsigs.bin (id -> alt nodes), calls.tsv (id -> type/size/vaf), in lockstep.
  FILE* af=std::fopen((dir+"/calls.altsigs.bin").c_str(),"rb"); if(!af)die("open calls.altsigs.bin");
  FILE* tf=std::fopen((dir+"/calls.tsv").c_str(),"rb"); if(!tf)die("open calls.tsv");
  std::vector<Call> calls; ankerl::unordered_dense::set<uint64_t> needed; needed.reserve(1u<<20);
  char* line=nullptr; size_t cap=0; uint64_t id; uint32_t ns; std::vector<uint64_t> buf;
  while(std::fread(&id,8,1,af)==1){ std::fread(&ns,4,1,af); buf.resize(ns); if(ns)std::fread(buf.data(),8,ns,af);
    uint64_t tid; uint8_t tc; long size,jlen,sup; float vaf;
    do{ if(getline(&line,&cap,tf)<=0)die("calls.tsv shorter than altsigs"); }while(!parse_call(line,tid,tc,size,vaf,jlen,sup));
    if(tid!=id)die("calls.tsv / altsigs id mismatch (regenerate both with one `svp call`)");
    // alt nodes; anchor candidates = colour {1,3} neighbours not in this call's alt set
    std::vector<uint32_t> alt; alt.reserve(ns);
    for(uint32_t i=0;i<ns;i++){ auto it=sig2node.find(buf[i]); if(it!=sig2node.end())alt.push_back(it->second); }
    ankerl::unordered_dense::set<uint32_t> altset(alt.begin(),alt.end()), as;
    for(uint32_t n:alt){
      for(uint32_t k=g.out_off[n];k<g.out_off[n+1];k++){ uint32_t w=g.out_nbr[k]; if((g.color[w]==1||g.color[w]==3)&&!altset.count(w))as.insert(w); }
      for(uint32_t k=g.in_off[n];k<g.in_off[n+1];k++){ uint32_t w=g.in_nbr[k]; if((g.color[w]==1||g.color[w]==3)&&!altset.count(w))as.insert(w); }
    }
    Call c; c.id=id; c.tc=tc; c.size=size; c.jlen=jlen; c.sup=sup; c.vaf=vaf;
    for(uint32_t w:as){ c.anchors.push_back(w); needed.insert(csig(g,w)); }
    calls.push_back(std::move(c));
  }
  std::fclose(af); std::fclose(tf); std::free(line);
  std::fprintf(stderr,"[anchor] calls=%zu anchor-sigs-needed=%zu\n",calls.size(),needed.size());

  // index ref: needed sig -> (contig<<40 | fwdpos). 2nd distinct pos => ambiguous (dropped).
  ankerl::unordered_dense::map<uint64_t,uint64_t> rpos; ankerl::unordered_dense::set<uint64_t> ambig;
  ankerl::unordered_dense::map<std::string,uint32_t> cidx; for(uint32_t i=0;i<FAORDER.size();i++)cidx[FAORDER[i]]=i;
  {
    FILE* fa=std::fopen(ref.c_str(),"rb"); if(!fa)die("open ref.fasta"); char* l=nullptr; size_t cp=0; ssize_t L;
    std::string seq,name; uint16_t cur=0;
    while((L=getline(&l,&cp,fa))>0){
      if(l[0]=='>'){ index_seq(seq,cur,needed,rpos,ambig); name.assign(l+1); size_t sp=name.find_first_of(" \t\r\n"); if(sp!=std::string::npos)name.resize(sp);
        auto it=cidx.find(name); cur=(it==cidx.end())?0:it->second; }
      else for(ssize_t i=0;i<L;i++){ char ch=l[i]; if(ch!='\n'&&ch!='\r')seq.push_back(ch>='a'&&ch<='z'?(char)(ch-32):ch); }
    }
    index_seq(seq,cur,needed,rpos,ambig); std::free(l); std::fclose(fa);
  }
  std::fprintf(stderr,"[anchor] ref indexed: located=%zu ambiguous=%zu\n",rpos.size(),ambig.size());

  // bracket each call -> VCF record on the consensus contig.
  std::vector<Rec> recs; recs.reserve(calls.size()); uint64_t unanch=0;
  std::vector<Hit> hits;
  for(const Call& c:calls){
    hits.clear();
    for(uint32_t w:c.anchors){ auto it=rpos.find(csig(g,w)); if(it==rpos.end())continue;
      Hit h; h.c=(uint16_t)(it->second>>40); h.s=it->second&((1ull<<40)-1); h.sl=g.seqlen[w]; hits.push_back(h); }
    if(hits.empty()){ unanch++; continue; }
    std::sort(hits.begin(),hits.end(),hit_lt);
    uint16_t bestc=hits[0].c; size_t bi=0,bj=0,run=0;
    for(size_t i=0;i<hits.size();){ size_t j=i; while(j<hits.size()&&hits[j].c==hits[i].c)j++;
      if(j-i>run){run=j-i;bestc=hits[i].c;bi=i;bj=j;} i=j; }
    // sorted within bestc: bi=leftmost anchor (min start), bj-1=rightmost. The breakpoint sits just
    // PAST the left anchor's core, so the true coordinate is start+seqlen (not start) — fixes the
    // ~1-core upstream offset. Single anchor: keep its start (side unknown).
    uint64_t lo=hits[bi].s, hi=hits[bj-1].s; uint16_t losl=hits[bi].sl;
    uint8_t tc=c.tc; bool retyped=false;
    // Retype BND via anchor geometry: the caller routes large deletions AND insertions through the
    // junction path (they don't form clean bubbles).
    //  DEL: spanned reference (50-50kb) is mostly NORMAL-ONLY (colour-1) — present in normal, gone in
    //       tumor (not a DUP=shared or artifact). Strand-agnostic.
    //  INS: the two anchors are ADJACENT (span<50) with a long novel junction (jlen>=50) = a point
    //       insertion of novel sequence. (Must precede the LOCALSUB test, which else eats it.)
    if(tc>=3 && run>=2){ long sspan=(long)(hi-lo);
      if(sspan>=50 && sspan<=50000 && !(2*c.jlen>=sspan)){
        std::string sp=fetch_range(FAORDER[bestc],(long)lo,sspan);
        if((int)sp.size()>=24){ int nn; const uint64_t *st,*en,*sg; extract(sp.data(),(int)sp.size(),0,LEVEL,&nn,&st,&en,&sg);
          int n1=0,nt=0; for(int i=0;i<nn;i++){ auto it=sig2node.find(sg[i]); if(it==sig2node.end())continue; nt++; if(g.color[it->second]==1)n1++; }
          if(nt>=3 && n1*10>=nt*6){ tc=2; retyped=true; } }
      } else if(sspan<50 && c.jlen>=50){ tc=1; retyped=true; }
    }
    // bubbles (SNP/INS/DEL, tc<=2) are bracketed by two flank cores: the breakpoint is just past the
    // LEFT flank, so start+seqlen (not start). BND is a point junction (unknown side) -> anchor start.
    uint64_t bp = (tc<=2 && run>=2) ? lo+losl : lo;
    long span=(long)(hi>bp ? hi-bp : 0);
    long dsize = !retyped ? c.size : (tc==2 ? span : std::max<long>(1,c.jlen-span));   // DEL=gap, INS=novel junction, bubble=Δ
    long svlen = (tc==1)?dsize : (tc==2)?-dsize : (tc==0)?1 : span;      // INS +, DEL -, SNP 1, BND span
    uint64_t pos1=bp+1, end = (tc==2)?(pos1+dsize) : (tc==1||tc==0)?pos1 : (pos1+(uint64_t)span);
    // A BND whose tumor junction REPRODUCES the spanned reference (jlen >= span/2) is not a novel
    // adjacency — it is a local substitution (a SNP/MNP the bubble path already calls), flag it out.
    uint8_t local = (tc>=3 && span>0 && 2*c.jlen >= span) ? 1 : 0;
    Rec r; r.cidx=bestc; r.pos=pos1; r.end=end; r.svlen=svlen; r.sup=c.sup; r.tc=tc; r.local=local; r.id=c.id; r.nanc=(uint32_t)run;
    r.rb=fetch1(FAORDER[bestc],(long)bp); recs.push_back(r);
  }
  // Dedup: a real SV gets several calls (bubble + BND + fwd/RC). Collapse same-class calls at the same
  // locus, keeping the highest-support representative (the redundant copies are the bulk of the FP).
  std::sort(recs.begin(),recs.end(),rec_locus_lt);
  std::vector<Rec> ded; ded.reserve(recs.size()); uint64_t before=recs.size();
  for(const Rec& r:recs){
    if(!ded.empty()){ Rec& k=ded.back();
      long W = (k.tc==0) ? 10 : std::max<long>(300,std::labs(k.svlen));
      if(k.cidx==r.cidx && cls(k.tc)==cls(r.tc) && (long)(r.pos-k.pos)<=W){ if(r.sup>k.sup)k=r; continue; }
    }
    ded.push_back(r);
  }
  recs.swap(ded);
  std::fprintf(stderr,"[anchor] dedup: %lu -> %lu records (dropped %lu redundant)\n",(unsigned long)before,(unsigned long)recs.size(),(unsigned long)(before-recs.size()));
  std::sort(recs.begin(),recs.end(),rec_pos_lt);

  // write VCF
  const std::string out=dir+"/calls.vcf", tmp=out+".tmp"; FILE* o=std::fopen(tmp.c_str(),"wb"); if(!o)die("open calls.vcf.tmp");
  std::fprintf(o,"##fileformat=VCFv4.2\n##source=svp-anchor\n##reference=%s\n",ref.c_str());
  for(const std::string& nm:FAORDER) std::fprintf(o,"##contig=<ID=%s,length=%lu>\n",nm.c_str(),(unsigned long)IX[nm].len);
  std::fprintf(o,"##ALT=<ID=DEL,Description=\"Deletion\">\n##ALT=<ID=INS,Description=\"Insertion\">\n");
  std::fprintf(o,"##ALT=<ID=BND,Description=\"Breakend (DUP/INV/translocation/large indel; subtype unresolved)\">\n");
  std::fprintf(o,"##ALT=<ID=SNP,Description=\"Single-nucleotide somatic variant (allele not base-resolved)\">\n");
  std::fprintf(o,"##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n");
  std::fprintf(o,"##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Signed length; spliced-delta for SNP/INS/DEL, anchored span for BND\">\n");
  std::fprintf(o,"##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">\n");
  std::fprintf(o,"##INFO=<ID=NANCHOR,Number=1,Type=Integer,Description=\"Reference anchor cores used to place this call\">\n");
  std::fprintf(o,"##INFO=<ID=VAF,Number=1,Type=Float,Description=\"Tumor variant allele fraction estimate\">\n");
  std::fprintf(o,"##INFO=<ID=CALLSRC,Number=1,Type=String,Description=\"svp caller source: BUBBLE/BND/TIP\">\n");
  std::fprintf(o,"##INFO=<ID=SUBTYPE,Number=1,Type=String,Description=\"BND subclass hint (CNgain = copy-number gain)\">\n");
  std::fprintf(o,"##FILTER=<ID=LOCALSUB,Description=\"BND junction reproduces the spanned reference (jlen>=span/2): a local SNP/MNP misfire, not a novel adjacency\">\n");
  std::fprintf(o,"##FILTER=<ID=TIP,Description=\"Tumor-only dead-end tip: insertion-end/breakend candidate, not a confirmed SV\">\n");
  std::fprintf(o,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
  for(const Rec& r:recs){ const char* src = r.tc<=2?"BUBBLE":(r.tc==5?"TIP":"BND");
    std::fprintf(o,"%s\t%lu\tsvp_%lu\t%c\t%s\t.\t%s\tSVTYPE=%s;SVLEN=%ld;END=%lu;NANCHOR=%u;CALLSRC=%s%s\n",
      FAORDER[r.cidx].c_str(),(unsigned long)r.pos,(unsigned long)r.id,r.rb,SVALT[r.tc],r.local?"LOCALSUB":(r.tc==5?"TIP":"PASS"),SVT[r.tc],r.svlen,
      (unsigned long)r.end,r.nanc,src, r.tc==4?";SUBTYPE=CNgain":""); }
  std::fclose(o); if(std::rename(tmp.c_str(),out.c_str())!=0)die("rename calls.vcf");
  std::fprintf(stderr,"[anchor] DONE wrote %zu records (unanchored=%lu) -> %s\n",recs.size(),(unsigned long)unanch,out.c_str());
}
