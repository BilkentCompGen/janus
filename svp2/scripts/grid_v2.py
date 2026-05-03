#!/usr/bin/env python3
"""Grid v2 — deeper exploration around the most promising configs from v1.

Focus:
  - DUP threshold, pctseq, refdist sweeps for dup + mixed sims
  - INV detection threshold (stage 5 knob) — but that requires code tweak,
    so we sweep by varying `aligner` variants (wfa vs wfa_rc vs mappy).
  - dedupe refdist/pctsize/pctseq to trade precision vs recall
  - chain_max_span sweep
"""
from __future__ import annotations
import itertools
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from grid_search import Config, run_config, ProcessPoolExecutor, as_completed  # type: ignore

import argparse
import json
import time
from dataclasses import asdict


def v2_configs() -> list[Config]:
    out: list[Config] = []
    # DUP tuning sweep — bubble over extractor × dup_thr × pctseq
    for ext in ("chain_rt", "chain_3w"):
        for aln in ("wfa", "wfa_rc"):
            for dup_thr in (0.45, 0.55, 0.60, 0.65, 0.75, 0.85):
                for pctseq in (0.50, 0.70):
                    out.append(Config(
                        sim="dup", extractor=ext, aligner=aln,
                        dup_thr=dup_thr, dedupe_pctseq=pctseq,
                        tag=f"dup_sw_d{dup_thr}_ps{pctseq}",
                    ))
    # MIXED tuning
    for ext in ("chain_rt", "chain_3w"):
        for aln in ("wfa", "wfa_rc", "edlib"):
            for dup_thr in (0.55, 0.70, 0.85):
                for dedupe_pctsize in (0.50, 0.70, 0.85):
                    out.append(Config(
                        sim="mixed_v2", extractor=ext, aligner=aln,
                        dup_thr=dup_thr, dedupe_pctsize=dedupe_pctsize,
                        tag=f"mix_d{dup_thr}_ds{dedupe_pctsize}",
                    ))
    # INV — sweep aligners at depth variations
    for depth in (4, 5, 6):
        for ext in ("chain_rt", "chain_3w"):
            for aln in ("wfa_rc", "mappy"):
                out.append(Config(
                    sim="inv", extractor=ext, aligner=aln, depth=depth,
                    tag=f"inv_d{depth}",
                ))
    # refdist sweep on all sims
    for sim in ("del", "ins", "dup", "inv", "mixed_v2"):
        for refdist in (100, 250, 1000, 2000):
            out.append(Config(
                sim=sim, extractor="chain_rt", aligner="wfa_rc",
                dedupe_refdist=refdist, tag=f"rd{refdist}",
            ))
    # chain_max_span sweep on mixed
    for cms in (50_000, 500_000, 1_000_000):
        out.append(Config(
            sim="mixed_v2", extractor="chain_rt", aligner="wfa_rc",
            chain_max_span=cms, tag=f"cms{cms}",
        ))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--parallel", type=int, default=6)
    ap.add_argument("--workers-stage5", type=int, default=4)
    ap.add_argument("--log", default="/home/yuio/svp/svp2/outputs/metrics/grid_v2.jsonl")
    ap.add_argument("--out-dir", default="/home/yuio/svp/svp2/outputs")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    cfgs = v2_configs()
    for c in cfgs:
        c.threads = args.workers_stage5

    already: set[str] = set()
    if os.path.exists(args.log):
        with open(args.log) as f:
            for line in f:
                try:
                    r = json.loads(line)
                    if "key" in r and r.get("summary"):
                        already.add(r["key"])
                except Exception:
                    pass
    cfgs = [c for c in cfgs if c.key() not in already]
    if args.limit:
        cfgs = cfgs[: args.limit]

    os.makedirs(os.path.dirname(args.log), exist_ok=True)
    print(f"[gridv2] running {len(cfgs)} new configs (skipping {len(already)})")
    if args.dry_run:
        for c in cfgs[:20]:
            print(c.key())
        if len(cfgs) > 20:
            print(f"... and {len(cfgs)-20} more")
        return 0

    t0 = time.time()
    with ProcessPoolExecutor(max_workers=args.parallel) as pool:
        futures = {pool.submit(run_config, c, args.out_dir, False): c for c in cfgs}
        done = 0
        with open(args.log, "a") as logf:
            for fut in as_completed(futures):
                cfg = futures[fut]
                try:
                    rec = fut.result()
                except Exception as e:
                    rec = {"config": asdict(cfg), "error": f"exception:{e!r}"}
                logf.write(json.dumps(rec) + "\n")
                logf.flush()
                done += 1
                s = rec.get("summary", {}) or {}
                f1 = float(s.get("f1") or 0)
                err = rec.get("error", "")
                print(f"[gridv2] {done}/{len(cfgs)} {cfg.key()[:100]} "
                      f"f1={f1:.3f} err={err}")
    print(f"[gridv2] done in {(time.time()-t0)/60:.1f} min")
    return 0


if __name__ == "__main__":
    sys.exit(main())
