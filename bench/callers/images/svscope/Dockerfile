FROM mambaorg/micromamba:1.5.10
USER root
RUN apt-get update && apt-get install -y --no-install-recommends git ca-certificates
RUN git clone https://github.com/Goatofmountain/SVscope.git /opt/SVScope && cd /opt/SVScope && git checkout 0ad3686
RUN micromamba install -y -n base -f /opt/SVScope/environment.yml && micromamba clean -a -y
ENV PATH=/opt/conda/bin:$PATH
WORKDIR /opt/SVScope
