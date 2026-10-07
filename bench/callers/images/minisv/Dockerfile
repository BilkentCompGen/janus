FROM debian:bookworm-slim
RUN apt-get update && apt-get install -y --no-install-recommends curl ca-certificates bzip2 samtools git coreutils
WORKDIR /opt
RUN curl -sL https://github.com/lh3/minimap2/releases/download/v2.31/minimap2-2.31_x64-linux.tar.bz2 | tar --no-same-owner -jxf - \
 && cp minimap2-2.31_x64-linux/minimap2 /usr/local/bin/ \
 && curl -sL https://github.com/lh3/minigraph/releases/download/v0.21/minigraph-0.21_x64-linux.tar.bz2 | tar --no-same-owner -jxf - \
 && cp minigraph-0.21_x64-linux/minigraph /usr/local/bin/ \
 && curl -sL https://github.com/attractivechaos/k8/releases/download/v1.2/k8-1.2.tar.bz2 | tar --no-same-owner -jxf - \
 && cp k8-1.2/k8-x86_64-Linux /usr/local/bin/k8 \
 && git clone https://github.com/lh3/minisv.git && cd minisv && git checkout d07d02e6 \
 && printf '#!/bin/sh\nexec k8 /opt/minisv/minisv.js "$@"\n' > /usr/local/bin/minisv.js && chmod +x /usr/local/bin/minisv.js
