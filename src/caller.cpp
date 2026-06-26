// svp call <dir>  (run_call): reference-free somatic variant caller — the final pipeline stage.
// Reference-free somatic variant caller over the 2-colour REFINED core graph (no simplify needed).
// HIGH-RECALL by design: call EVERY structure that carries a tumor-only (colour-2) allele; the
// false positives are characterised/filtered downstream. Uses ONLY graph-intrinsic signals.
//   colour: 1=normal-only(N) 2=tumor-only(T) 3=shared(S). per-node ncov/tcov. topology. payload.
//
// Call sources (a somatic allele MUST involve >=1 colour-2 arm/core; N/N,N/S,S/S carry 0 variants):
//   BUBBLE  : (S,T) two-arm bubble with >=1 tumor (colour-2) arm -> T/N, T/S, T/T. alt=tumor arm,
//             ref=other arm. type/size from signed spliced Δ (SNP |Δ|<=eps, INS Δ>0, DEL Δ<0).
//   BND     : tumor-only junction core (colour-2, ncov==0) in a linear chain between two shared
//             cores = tumor-private adjacency (no reference counterpart) -> DUP/INV/transloc breakend.
//   TIP     : tumor-only dead-end chain (colour-2) -> candidate insertion-end / breakend.
// Outputs:  <prefix>.tsv (features+sequences, call_id keyed) and <prefix>.altsigs.bin
//           (per call: u64 id, u32 nsig, u64 sig[] of the alt-allele cores) for recall/precision scoring.
// fwd+RC strand duplicates collapsed by canonical (min seq/revcomp) haplotype hash.
#define XXH_INLINE_ALL
#include <xxhash.h>
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <unordered_set>

namespace {
int MINSUP = 3, SNP_EPS = 2;
uint8_t COMP[256];
void initc(){ for(int i=0;i<256;i++)COMP[i]='N'; COMP['A']='T';COMP['T']='A';COMP['C']='G';COMP['G']='C'; }
bool exists(const std::string& p){ FILE* f=std::fopen(p.c_str(),"rb"); if(f){std::fclose(f);return true;} return false; }

struct G{
  uint64_t K=0,E=0;
  std::vector<uint8_t> color; std::vector<uint32_t> ncov,tcov; std::vector<uint16_t> seqlen;
  std::vector<uint8_t> payload; std::vector<uint64_t> off;
  std::vector<uint32_t> out_off,out_nbr,in_off,in_nbr,live_out,live_in;
};
void load(const char* gp, G& g){
  FILE* f=std::fopen(gp,"rb"); if(!f){std::fprintf(stderr,"[caller] FATAL open %s\n",gp);std::exit(1);}
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
  g.live_out.resize(K); g.live_in.resize(K);
  for(uint64_t i=0;i<K;i++){g.live_out[i]=g.out_off[i+1]-g.out_off[i]; g.live_in[i]=g.in_off[i+1]-g.in_off[i];}
  std::fprintf(stderr,"[caller] K=%lu E=%lu payload=%.1fGB\n",(unsigned long)K,(unsigned long)E,total/1e9);
}
inline const uint8_t* P(const G&g,uint32_t x){ return g.payload.data()+g.off[x]; }
inline uint32_t one_out(const G&g,uint32_t v){ return g.out_nbr[g.out_off[v]]; }
inline uint32_t one_in (const G&g,uint32_t v){ return g.in_nbr[g.in_off[v]]; }
inline uint64_t csig(const G&g,uint32_t x){ return XXH3_64bits(P(g,x),g.seqlen[x]); }

uint32_t trace_arm(const G&g,uint32_t start,std::vector<uint32_t>&arm){
  arm.clear(); uint32_t cur=start;
  for(int s=0;s<1000000;s++){ arm.push_back(cur);
    if(g.live_out[cur]!=1) return UINT32_MAX; uint32_t nxt=one_out(g,cur);
    if(g.live_in[nxt]>=2) return nxt; if(g.live_in[nxt]!=1) return UINT32_MAX; cur=nxt; }
  return UINT32_MAX;
}
uint32_t trace_tip(const G&g,uint32_t v,int dir,std::vector<uint32_t>&chain){
  chain.clear(); uint32_t u=v;
  for(int s=0;s<100000;s++){ chain.push_back(u);
    uint32_t deg=dir>0?g.live_in[u]:g.live_out[u]; if(deg==0)return UINT32_MAX; if(deg!=1)return u;
    uint32_t p=dir>0?one_in(g,u):one_out(g,u); uint32_t pb=dir>0?g.live_out[p]:g.live_in[p];
    if(pb>=2)return p; u=p; }
  return UINT32_MAX;
}
int arm_color(const G&g,const std::vector<uint32_t>&a){ uint32_t c[4]={0,0,0,0}; for(uint32_t x:a)c[g.color[x]]++;
  if(c[2]>=c[1]&&c[2]>=c[3])return 2; if(c[1]>=c[3])return 1; return 3; }
uint32_t sp(const uint8_t*A,uint32_t la,const uint8_t*B,uint32_t lb){
  uint32_t m=std::min(la,lb); std::vector<uint32_t> fail(lb+1,0);
  for(uint32_t i=1,k=0;i<lb;i++){ while(k&&B[i]!=B[k])k=fail[k]; if(B[i]==B[k])k++; fail[i+1]=k; }
  uint32_t k=0; for(uint32_t i=la-m;i<la;i++){ while(k&&A[i]!=B[k])k=fail[k]; if(A[i]==B[k])k++; }
  return k>m?m:k;
}
std::string splice(const G&g,const std::vector<uint32_t>&nodes){
  std::string out; if(nodes.empty())return out;
  out.assign((const char*)P(g,nodes[0]), g.seqlen[nodes[0]]);
  for(size_t i=1;i<nodes.size();i++){ uint32_t a=nodes[i-1],b=nodes[i];
    uint32_t ov=sp(P(g,a),g.seqlen[a],P(g,b),g.seqlen[b]); out.append((const char*)P(g,b)+ov, g.seqlen[b]-ov); }
  return out;
}
std::string rc(const std::string&s){ std::string r(s.size(),'N'); for(size_t j=0;j<s.size();j++)r[j]=(char)COMP[(uint8_t)s[s.size()-1-j]]; return r; }
uint64_t canon_key(const std::string&s){ std::string r=rc(s); const std::string&c=(s<=r)?s:r; return XXH3_64bits(c.data(),c.size()); }
long arm_pathlen(const G&g,uint32_t S,const std::vector<uint32_t>&arm,uint32_t T){
  long L=0; uint32_t prev=S; for(uint32_t c:arm){ uint32_t ov=sp(P(g,prev),g.seqlen[prev],P(g,c),g.seqlen[c]); L+=(long)g.seqlen[c]-ov; prev=c; }
  uint32_t ov=sp(P(g,prev),g.seqlen[prev],P(g,T),g.seqlen[T]); L+=(long)g.seqlen[T]-ov; return L;
}
uint32_t arm_min_tcov(const G&g,const std::vector<uint32_t>&a){ uint32_t m=UINT32_MAX; for(uint32_t x:a)m=std::min(m,g.tcov[x]); return a.empty()?0:m; }
uint32_t arm_mean_tcov(const G&g,const std::vector<uint32_t>&a){ uint64_t s=0; for(uint32_t x:a)s+=g.tcov[x]; return a.empty()?0:(uint32_t)(s/a.size()); }

// From a colour-2 seed, walk `dir` (-1 via in-edges, +1 via out-edges) through the highest-tcov
// tumor-only cores until a SHARED (colour-3) anchor; collect the tumor-only run + its min tcov.
// Returns the shared anchor, or UINT32_MAX if none within jmax (not a clean junction to a real locus).
uint32_t reach_shared(const G&g,const std::vector<uint8_t>&consumed,uint32_t seed,int dir,int jmax,
                      std::vector<uint32_t>&run,uint32_t&mintc){
  run.clear(); uint32_t cur=seed;
  for(int s=0;s<jmax;s++){
    uint32_t a0=dir<0?g.in_off[cur]:g.out_off[cur], a1=dir<0?g.in_off[cur+1]:g.out_off[cur+1];
    const std::vector<uint32_t>& nbr=dir<0?g.in_nbr:g.out_nbr;
    uint32_t bs=UINT32_MAX,bst=0, bt=UINT32_MAX,btt=0;
    for(uint32_t k=a0;k<a1;k++){ uint32_t w=nbr[k];
      if(g.color[w]==3){ if(g.tcov[w]>=bst){bs=w;bst=g.tcov[w];} }
      else if(g.color[w]==2&&!consumed[w]&&w!=seed){ if(g.tcov[w]>=btt){bt=w;btt=g.tcov[w];} } }
    if(bs!=UINT32_MAX) return bs;            // reached a shared (real) anchor
    if(bt==UINT32_MAX) return UINT32_MAX;    // no tumor continuation -> not a clean bridge
    run.push_back(bt); mintc=std::min(mintc,g.tcov[bt]); cur=bt;
  }
  return UINT32_MAX;
}
FILE* SIGF;
void emit_sigs(uint64_t id,const G&g,const std::vector<uint32_t>&nodes){
  uint32_t n=nodes.size(); std::fwrite(&id,8,1,SIGF); std::fwrite(&n,4,1,SIGF);
  for(uint32_t x:nodes){ uint64_t s=csig(g,x); std::fwrite(&s,8,1,SIGF); }
}
}  // namespace

// dir/refined.bin (else graph.bin) -> dir/calls.tsv + dir/calls.altsigs.bin. minsup via CALL_MINSUP (default 3).
void run_call(const std::string& dir){
  if(const char* e=std::getenv("CALL_MINSUP")) MINSUP=std::atoi(e);
  const std::string in = exists(dir+"/refined.bin") ? dir+"/refined.bin" : dir+"/graph.bin";
  initc(); G g; load(in.c_str(),g); const uint64_t K=g.K;
  std::vector<uint8_t> consumed(K,0);
  std::string pre=dir+"/calls";
  FILE* o=std::fopen((pre+".tsv").c_str(),"wb"); std::setvbuf(o,nullptr,_IOFBF,1<<22);
  SIGF=std::fopen((pre+".altsigs.bin").c_str(),"wb"); std::setvbuf(SIGF,nullptr,_IOFBF,1<<22);
  std::fprintf(o,"call_id\tsource\tsubtype\ttype\tsize\talt_tcov\talt_meantcov\tanchor_tcov\tanchor_ncov\tvaf\talt_len\tref_len\tdelta\tlflank\talt_seq\tref_seq\trflank\n");
  std::unordered_set<uint64_t> seen; seen.reserve(1u<<22);
  uint64_t id=0; uint64_t cnt_bub=0,cnt_bnd=0,cnt_tip=0, raw_bub=0;
  uint64_t by_sub[3]={0}; const char* SUB[3]={"T/N","T/S","T/T"};

  // ---- Pass A: bubbles with >=1 tumor (colour-2) arm ----
  std::vector<uint32_t> a1,a2,succ;
  for(uint64_t s=0;s<K;s++){ if(g.live_out[s]<2)continue;
    succ.clear(); for(uint32_t k=g.out_off[s];k<g.out_off[s+1];k++)succ.push_back(g.out_nbr[k]);
    for(size_t i=0;i<succ.size();i++)for(size_t j=i+1;j<succ.size();j++){
      uint32_t T1=trace_arm(g,succ[i],a1); if(T1==UINT32_MAX)continue;
      uint32_t T2=trace_arm(g,succ[j],a2); if(T2==UINT32_MAX||T2!=T1)continue;
      int c1=arm_color(g,a1),c2=arm_color(g,a2);
      int nt=(c1==2)+(c2==2); if(nt==0) continue;                    // need >=1 tumor arm (else 0 variants)
      // alt = a tumor arm; ref = the other arm. T/T: alt=a1, ref=a2.
      std::vector<uint32_t> *alt,*ref; int sub;
      if(c1==2&&c2==2){ alt=&a1; ref=&a2; sub=2; }                   // T/T
      else if(c1==2){ alt=&a1; ref=&a2; sub=(c2==1)?0:1; }          // c2 in {1->T/N,3->T/S}
      else { alt=&a2; ref=&a1; sub=(c1==1)?0:1; }
      uint32_t sup=arm_min_tcov(g,*alt); if(sup<(uint32_t)MINSUP) continue;
      raw_bub++;
      long d = arm_pathlen(g,(uint32_t)s,*alt,T1) - arm_pathlen(g,(uint32_t)s,*ref,T1);
      std::vector<uint32_t> hap; hap.push_back((uint32_t)s); for(uint32_t x:*alt)hap.push_back(x); hap.push_back(T1);
      std::string hs=splice(g,hap); uint64_t key=canon_key(hs);
      if(!seen.insert(key).second) continue;
      for(uint32_t x:*alt)consumed[x]=1; for(uint32_t x:*ref)consumed[x]=1;
      const char* ty; long size;
      if(std::labs(d)<=SNP_EPS){ ty="SNP"; size=1; } else if(d>0){ ty="INS"; size=d; } else { ty="DEL"; size=-d; }
      uint32_t anc_t=std::min(g.tcov[s],g.tcov[T1]), anc_n=std::min(g.ncov[s],g.ncov[T1]);
      double vaf = anc_t>0 ? (double)arm_mean_tcov(g,*alt)/anc_t : 0;
      std::string altseq=splice(g,*alt), refseq=splice(g,*ref);
      std::string lf((const char*)P(g,(uint32_t)s),std::min<uint32_t>(g.seqlen[s],31));
      std::string rf((const char*)P(g,T1),std::min<uint32_t>(g.seqlen[T1],31));
      std::fprintf(o,"%lu\tBUBBLE\t%s\t%s\t%ld\t%u\t%u\t%u\t%u\t%.3f\t%ld\t%ld\t%ld\t%s\t%.200s\t%.200s\t%s\n",
        (unsigned long)id,SUB[sub],ty,size,sup,arm_mean_tcov(g,*alt),anc_t,anc_n,vaf,
        arm_pathlen(g,(uint32_t)s,*alt,T1),arm_pathlen(g,(uint32_t)s,*ref,T1),d,lf.c_str(),altseq.c_str(),refseq.c_str(),rf.c_str());
      emit_sigs(id,g,*alt); id++; cnt_bub++; by_sub[sub]++;
    }
  }
  // ---- Pass B: breakend junctions (DUP/INV/transloc/large-indel) ----
  // A breakpoint = a tumor-only run that BRIDGES two shared (real) loci = a novel adjacency with no
  // reference counterpart (so it never forms a bubble). For each unconsumed tumor core, trace the
  // highest-depth tumor path to the nearest shared anchor on each side; if both land on shared cores
  // AND the junction runs at ~full tumor depth (real homozygous junction, not a low-cov error spur),
  // emit it. Works inside the repeat tangle (no degree-1 requirement) — this is what recovers DUP/INV.
  const int JMAX = std::getenv("CALL_JMAX") ? std::atoi(std::getenv("CALL_JMAX")) : 12;
  const double DEPTH_FRAC = std::getenv("CALL_DEPTH_FRAC") ? std::atof(std::getenv("CALL_DEPTH_FRAC")) : 0.30;
  std::vector<uint32_t> up,dn,jrun;
  for(uint64_t v=0;v<K;v++){
    if(consumed[v]||g.color[v]!=2||g.tcov[v]<(uint32_t)MINSUP) continue;
    uint32_t mintc=g.tcov[v];
    uint32_t L=reach_shared(g,consumed,(uint32_t)v,-1,JMAX,up,mintc);
    uint32_t R=reach_shared(g,consumed,(uint32_t)v,+1,JMAX,dn,mintc);
    if(L==UINT32_MAX||R==UINT32_MAX) continue;                    // not a clean two-anchor bridge
    jrun.clear(); for(auto it=up.rbegin();it!=up.rend();++it)jrun.push_back(*it); jrun.push_back((uint32_t)v); for(uint32_t x:dn)jrun.push_back(x);
    if(jrun.size()>(size_t)JMAX) continue;
    uint32_t anc=std::min(g.tcov[L],g.tcov[R]);
    if((double)mintc < DEPTH_FRAC*anc) continue;                  // junction must run at ~full depth (reject error spurs)
    bool dup=false; for(uint32_t x:jrun) if(consumed[x])dup=true; if(dup)continue;
    std::string js=splice(g,jrun); uint64_t key=canon_key(js); if(!seen.insert(key).second) continue;
    bool cn=(g.ncov[L]>0&&(double)g.tcov[L]/g.ncov[L]>=1.5)||(g.ncov[R]>0&&(double)g.tcov[R]/g.ncov[R]>=1.5);
    std::string lf((const char*)P(g,L),std::min<uint32_t>(g.seqlen[L],31)), rf((const char*)P(g,R),std::min<uint32_t>(g.seqlen[R],31));
    std::fprintf(o,"%lu\tBND\t-\t%s\t%zu\t%u\t%u\t%u\t%u\t0\t%zu\t0\t0\t%s\t%.200s\t-\t%s\n",
      (unsigned long)id, cn?"BND_CNgain":"BND", jrun.size(), mintc, arm_mean_tcov(g,jrun), anc, std::min(g.ncov[L],g.ncov[R]), js.size(), lf.c_str(), js.c_str(), rf.c_str());
    emit_sigs(id,g,jrun); for(uint32_t x:jrun)consumed[x]=1; id++; cnt_bnd++;
  }
  // ---- Pass C: tumor-only dead-end tips (insertion-end / breakend candidates) ----
  std::vector<uint32_t> chain;
  for(uint64_t v=0;v<K;v++){ if(consumed[v]||g.color[v]!=2)continue;
    for(int dir=+1;dir>=-1;dir-=2){ uint32_t ed=dir>0?g.live_out[v]:g.live_in[v]; if(ed!=0)continue;
      uint32_t jx=trace_tip(g,(uint32_t)v,dir,chain); if(chain.size()>16)continue;
      // all-tumor tip
      bool allt=true; uint32_t sup=UINT32_MAX; for(uint32_t x:chain){ if(g.color[x]!=2)allt=false; sup=std::min(sup,g.tcov[x]); }
      if(!allt||sup<(uint32_t)MINSUP) continue;
      bool dup=false; for(uint32_t x:chain) if(consumed[x])dup=true; if(dup)continue;
      std::string ts=splice(g,chain); uint64_t key=canon_key(ts); if(!seen.insert(key).second) continue;
      const char* anchctx = (jx==UINT32_MAX)?"-":""; (void)anchctx;
      std::string lf="-"; if(jx!=UINT32_MAX){ lf.assign((const char*)P(g,jx),std::min<uint32_t>(g.seqlen[jx],31)); }
      std::fprintf(o,"%lu\tTIP\t-\tTIP\t%zu\t%u\t%u\t0\t0\t0\t%zu\t0\t0\t%s\t%.200s\t-\t-\n",
        (unsigned long)id, ts.size(), sup, arm_mean_tcov(g,chain), ts.size(), lf.c_str(), ts.c_str());
      emit_sigs(id,g,chain); for(uint32_t x:chain)consumed[x]=1; id++; cnt_tip++;
    }
  }
  std::fclose(o); std::fclose(SIGF);
  std::fprintf(stderr,"[caller] BUBBLE=%lu (T/N=%lu T/S=%lu T/T=%lu) BND=%lu TIP=%lu  total=%lu  (raw bubbles pre-dedup=%lu) minsup=%d -> %s.{tsv,altsigs.bin}\n",
    (unsigned long)cnt_bub,(unsigned long)by_sub[0],(unsigned long)by_sub[1],(unsigned long)by_sub[2],
    (unsigned long)cnt_bnd,(unsigned long)cnt_tip,(unsigned long)id,(unsigned long)raw_bub,MINSUP,pre.c_str());
}
