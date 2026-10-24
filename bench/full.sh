#!/bin/bash
# Fig. 2a, 2c, Fig. 4 and the tables: janus on the seven pairs (five CASTLE pairs, HG008, COLO829), measured, into
# runs/full/<pair>. The calls are classified afterwards by classifier.sh (each CASTLE pair held out).
#   bench/full.sh [pair ...]
source "$(dirname "$0")/env.sh"
for p in ${*:-$LABELED}; do
    mkdir -p $R/full; sheet $R/full/$p.tsv $(src $p n) $(src $p t)
    janus_run full $R/full/$p $p $R/full/$p.tsv ""
done
