# Benchmark Results

All benchmarks: depth=5, chr19, simulation reference = `simulation/ref.fasta` (~58MB human chr19). Evaluated with [Truvari](https://github.com/ACEnglish/truvari) (SV-level, >=50bp).

## 2026-04-16

**Config:** `chr19_mixed_v2.yaml` — `min_intersv_dist: 200`, `allow_sv_overlap: True` for SNPs/short indels, `max_tries: 300`.

### Truvari (mixed, per-type)

| Type  |  TP |  FP |  FN | Precision | Recall |     F1 |
|-------|-----|-----|-----|-----------|--------|--------|
| DEL   | 118 |  79 |   9 |    0.5990 | 0.9291 | 0.7284 |
| INS   | 121 |  60 |   6 |    0.6685 | 0.9528 | 0.7857 |
| DUP   |  88 |  43 |  43 |    0.6718 | 0.6718 | 0.6718 |
| INV   | 129 |  58 |   7 |    0.6898 | 0.9485 | 0.7988 |
| ALL   | 456 | 272 |  65 |    0.6264 | 0.8752 | 0.7302 |

---

## 2026-04-02

### Truvari (mixed, per-type)

| Type  |  TP |  FP |  FN | Precision | Recall |     F1 |
|-------|-----|-----|-----|-----------|--------|--------|
| DEL   | 133 |   1 |  10 |    0.9925 | 0.9301 | 0.9603 |
| INS   | 140 |   1 |   7 |    0.9929 | 0.9524 | 0.9722 |
| DUP   |  77 |  37 |  65 |    0.6754 | 0.5423 | 0.6016 |
| INV   | 115 |  65 |  36 |    0.6389 | 0.7616 | 0.6949 |
| ALL   | 473 | 149 |  68 |    0.7605 | 0.8743 | 0.8134 |

---

## 2026-03-26

### Truvari (mixed, per-type)

| Type  |  TP |  FP |  FN | Precision | Recall |     F1 |
|-------|-----|-----|-----|-----------|--------|--------|
| DEL   | 135 |   1 |   8 |    0.9926 | 0.9441 | 0.9677 |
| INS   | 145 |   2 |   2 |    0.9864 | 0.9864 | 0.9864 |
| DUP   |  28 |  94 | 114 |    0.2295 | 0.1972 | 0.2121 |
| ALL   | 288 | 105 | 106 |    0.7328 | 0.7310 | 0.7319 |
