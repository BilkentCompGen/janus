FROM snape-bench/svision-pro:2.5-pysam0.21
RUN pip install --no-cache-dir torch==1.10.1+cu113 -f https://download.pytorch.org/whl/cu113/torch_stable.html
RUN sed -i "s/set_device(options.gpu_id)/set_device(int(options.gpu_id))/" /usr/local/lib/python3.7/site-packages/src/inheritype_op.py
