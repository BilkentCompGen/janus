FROM python:3.11-slim-bookworm
RUN apt-get update && apt-get install -y --no-install-recommends git samtools minimap2 ca-certificates
RUN git clone https://github.com/eioyuou/SomaSV.git /opt/SomaSV && cd /opt/SomaSV && git checkout 92f1560 && pip install --no-cache-dir -e .
WORKDIR /opt/SomaSV
