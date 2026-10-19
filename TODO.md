# TODO

- Site detection on long chains: a marked chain gives a site only if a read traverses the whole chain, so a chain that no read covers end to end is dropped. Consider splitting it into two sites instead (A -> XY... on entry and ...Z -> B on exit). Tested on the five CASTLE pairs (2026-10-12): 1,550 to 2,170 such chains per pair, 491 to 842 of them entered or left by enough tumor reads; with the two-site rule no consensus SV gained a record in any pair, and the PASS calls matching the consensus changed by 0, 0, -1, 0, -1 (HCC1954, NCI-H2009, HCC1937, NCI-H1437, HCC1395). Not adopted for now.

- Anchor uniqueness on the reference: during the alignment step, check whether the anchors of a site are also unique on the reference, and use this as an additional signal. Not tested yet.

- ONT reads with homopolymer compression: LCP cores on ONT reads were tried once without homopolymer compression and failed. Try again with the compression that janus applies, on recent ONT data, and at lower LCP levels than 5, since the share of error-free cores falls quickly with core length. No ONT reads are on this machine yet (only the ONT call sets of the CASTLE panel). Not tested yet.
