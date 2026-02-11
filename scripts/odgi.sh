odgi() {
    docker run --rm -v "$PWD":/data -w /data quay.io/biocontainers/odgi:0.9.4--h077b44d_0 odgi "$@"
}

GRAPH_NAME="${1:-test}"

echo "Building graph..."
odgi build -g outputs/gfa/${GRAPH_NAME}_graph.gfa -o outputs/og/${GRAPH_NAME}_graph.og
echo "Unchopping graph..."
odgi unchop -i outputs/og/${GRAPH_NAME}_graph.og -o outputs/og/${GRAPH_NAME}_compressed.og
echo "Viewing graph..."
odgi view -i outputs/og/${GRAPH_NAME}_compressed.og -g > outputs/gfa/${GRAPH_NAME}_compressed.gfa
echo "Done!"