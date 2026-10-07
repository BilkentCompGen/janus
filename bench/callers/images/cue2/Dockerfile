FROM pytorch/pytorch:2.0.1-cuda11.7-cudnn8-runtime
RUN apt-get update && apt-get install -y --no-install-recommends git libgl1 libglib2.0-0 build-essential && apt-get clean
COPY cue2 /opt/cue2
RUN cd /opt/cue2 && pip install --no-cache-dir . && mkdir -p /opt/models && \
    python -c "import urllib.request; urllib.request.urlretrieve('https://storage.googleapis.com/cue-models/latest/cue2.v2.pt', '/opt/models/cue2.v2.pt')"
ENV PYTHONPATH=/opt/cue2
