#!/usr/bin/env python3
"""Turn `trace probeseq` output (independent re-derivation of each hand-built allele's node
colours) + verify_targets.json into claim-level verdicts for the report.

probeseq columns: label n_cores n_present n_color1 n_color2 n_color3 n_outonly n_error
Usage: verify_interpret.py <probeseq_out.tsv> <verify_targets.json> <verify_findings.json>
"""
import sys, csv, json, collections

PROBE, TARGETS, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
tg = json.load(open(TARGETS))
P = {}
with open(PROBE) as f:
    for r in csv.DictReader(f, delimiter='\t'):
        P[r['label']] = {k:int(r[k]) for k in r if k!='label'}

def supports(meta, p):
    """does this single probed sequence support its expectation?"""
    pres = p['n_present']; exp = meta['expect']
    if exp == 'c2':   # a novel somatic seam -> at least one out-only core present
        return p['n_outonly'] > 0
    if exp == 'c3':   # shared body/interior -> colour-3 dominates present, little somatic leak
        return pres > 0 and p['n_color3'] >= 0.8*pres and p['n_outonly'] <= max(1, 0.1*pres)
    if exp == 'c1':   # deleted segment -> normal-only colour present and dominant over somatic
        return p['n_color1'] > 0 and p['n_color1'] >= p['n_outonly']
    return False

# group labels by claim
CLAIMS = [
 ('snp',   'SNP alt alleles create colour-2 somatic cores at the substitution',
           lambda m: m['type']=='SNP' and m['role']=='seam'),
 ('indel_s','Short DEL/INS junctions are colour-2 seams (type not encoded — look insertion-like)',
           lambda m: m['role']=='seam' and m['type'] in ('DEL','INS') and m['svlen']<25),
 ('del_l', 'Large DEL exposes a colour-1 (normal-only) deleted-segment arm',
           lambda m: m['type']=='DEL' and m['role']=='delbody'),
 ('ins_l', 'Large INS body is colour-2 (tumor-only) inserted sequence',
           lambda m: m['type']=='INS' and m['role']=='seam' and m['svlen']>=1000),
 ('dup_seam','Long DUP junction is a colour-2 seam (the only novel signal)',
           lambda m: m['type']=='DUP' and m['role']=='seam'),
 ('dup_body','Long DUP body is colour-3 SHARED — collapses onto the backbone, NOT a novel cycle',
           lambda m: m['type']=='DUP' and m['role']=='body'),
 ('inv_bp','INV breakpoints are colour-2 seams',
           lambda m: m['type']=='INV' and m['role']=='seam'),
 ('inv_int','INV interior is colour-3 SHARED (inverted body, same up to strand)',
           lambda m: m['type']=='INV' and m['role']=='interior'),
]

findings = []
print(f"{'claim':10} {'n':>3} {'ok':>3}  verdict     evidence")
for key, text, sel in CLAIMS:
    labs = [l for l,m in tg.items() if sel(m) and l in P]
    if not labs: continue
    ok = sum(1 for l in labs if supports(tg[l], P[l]))
    n = len(labs); frac = ok/n
    verdict = 'supported' if frac>=0.8 else ('partial' if frac>=0.5 else 'refuted')
    # evidence: median colour mix
    def med(field): vals=sorted(P[l][field] for l in labs); return vals[len(vals)//2]
    if tg[labs[0]]['expect']=='c2':
        ev = f"{ok}/{n} probed alleles place a colour-2/out-only core at the junction (median out-only cores/seq = {med('n_outonly')})."
    elif tg[labs[0]]['expect']=='c3':
        ev = f"{ok}/{n} bodies are colour-3 dominant (median present={med('n_present')}, colour-3={med('n_color3')}, out-only leak={med('n_outonly')})."
    else:
        ev = f"{ok}/{n} carry a colour-1 normal-only arm (median colour-1 cores={med('n_color1')}, out-only={med('n_outonly')})."
    findings.append({'key':key,'claim':text,'verdict':verdict,'evidence':ev,'n':n,'ok':ok})
    print(f"{key:10} {n:>3} {ok:>3}  {verdict:10}  {ev}")

json.dump(findings, open(OUT,'w'), indent=1)
print(f"\nwrote {OUT}  ({len(findings)} claims)")
