# Rewrite plan: compare on cores, then place

Measured on the five HiFi pairs (fair minda: snape as the fifth HiFi voter; NIST for hg008), AA2 loses F1 in three places that
thresholds cannot fix: it judges every site on one tumor read and one normal read, it decides "absent from the normal" on exact
core paths, and it never checks that the aligned event is the object the combination count belongs to. The graph detector and its
counts are sound (full-BAM support / combination count has median 1.00 for true and false calls alike). The rewrite keeps
extract, count, graph and prune, and changes what happens after detection.

## Decisions

1. **Asymmetric.** A combination is interesting when a tumor has it on at least 3 reads and the normal on none. Normal-only
   combinations are never evidence. Only tumor reads are cut and aligned.
2. **Shadow counts instead of a lower mm.** count also writes the cores it saw once or twice (`.cnt.low`). graph marks, for every
   node another sample kept, whether this sample holds it below the threshold (`sh`, one byte per node per sample). The graph
   itself is unchanged; the information enters the call as a feature.
3. **Sites and anchors.** A site is one interesting combination in one tumor. Its anchors are the nodes on each side of the
   combination on the chosen tumor read that all samples share and that occur once in that read, up to 48 per side, innermost
   first, with their positions on that read.
4. **Measure every read on cores, cut the minimum.** One pass over all reads of all samples: a read that contains an anchor on
   each side of a site (innermost present-and-unique pair) yields its anchor distance minus the tumor read's distance for the
   same pair. Tumor reads on the route are clustered by that difference (3 clusters, 20 bp); the best read per cluster, by flank,
   is the one cut. Normal reads are counted (spanning), their differences kept (up to 32), and none are cut. Reads off the route
   are counted as spanning.
5. **Site statistics carried as numbers.** Per site: tumor reads on the route (the combination count), tumor reads spanning,
   clusters and their sizes, normal reads spanning, normal reads whose difference is within max(20, |median|/2) of the tumor
   (germline), median / min / max normal difference, shadow presence and tumor-only node density between the anchors, branch
   density between the anchors, anchor span.
6. **call attaches events by anchor span**, cross-checks the aligned size against the measured difference, orders BND ends so
   reciprocal duplicates cluster, keeps only primary chromosomes, and emits every record with an error probability from a logistic
   model over the site statistics plus alignment features (MAPQ, divergence, reference repeat score, VNTR flag, type, size,
   relative depth). PASS is a cut on that probability. No CHREND, no hard NORMAL / NOREF / LOWSUP.
7. **Weights** are fitted leave-one-pair-out on the five pairs against the no-snape consensus (NIST for hg008), evaluated with
   the fair ensemble (snape as the fifth HiFi voter) and, for the Severus head-to-head, `minda_4_2_ws`. The shipped weights are
   the fit on all five.
8. **Not now:** the single streaming pass with paths kept in memory, graph-derived anchors, random-access cutting. F1 first.

## Stages and files

| stage | file | change |
|---|---|---|
| count | `count.cpp` | second output `<out>.low`: sorted hashes seen 1 or 2 times |
| graph | `graph.cpp`, `files.*` | `sh` plane appended to `.gph` |
| select | `select.cpp` | asymmetric gather, anchor pass with positions, measure pass, cut one read per cluster, site rows |
| select | `measure.*` | per-read anchor-distance measurement and cluster bookkeeping (replaces `mates.*`) |
| call | `call.cpp`, `main.cpp` | rows without normal pieces, features, logistic verdict, `-w weights`, `-v vntr.bed` |

## Evidence rows (`ev.tsv`)

`piece a b site counts depth key tumor clusters cluster_size tumor_span normal_span normal_same lest dmin dmax shadow tonly branches span`

## Order of work

count → graph/files → measure → select → call/main → build → one pair end to end → all five → fit → evaluate → commit.
