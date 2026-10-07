FROM nvidia/cuda:11.8.0-cudnn8-runtime-ubuntu22.04
RUN apt-get update && apt-get install -y --no-install-recommends python3-pip git && apt-get clean
RUN pip install --no-cache-dir tensorflow==2.12.1 pandas pysam scikit-learn==1.2.2 scipy
RUN git clone https://github.com/coopsor/SVDF.git /opt/SVDF && git -C /opt/SVDF checkout -q 41c1607
RUN sed -i "s#'model/#'/opt/SVDF/model/#" /opt/SVDF/*.py
