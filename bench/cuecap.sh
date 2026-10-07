#!/usr/bin/env bash
while true; do
    for c in $(docker ps -q --filter ancestor=snape-bench/cue2:latest); do
        [ "$(docker inspect -f '{{.HostConfig.Memory}}' $c)" -lt 165000000000 ] && docker update --memory 165g --memory-swap 165g $c > /dev/null
    done
    sleep 30
done
