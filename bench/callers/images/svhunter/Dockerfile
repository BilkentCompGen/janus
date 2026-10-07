FROM nvidia/cuda:11.8.0-cudnn8-runtime-ubuntu22.04
RUN apt-get update && apt-get install -y --no-install-recommends python3-pip git && apt-get clean
RUN pip install --no-cache-dir tensorflow==2.12.1 pandas pysam scikit-learn==1.5.1 scipy
RUN git clone https://github.com/eioyuou/SVHunter.git /opt/SVHunter && git -C /opt/SVHunter checkout -q 28d1ea5 && sed -i '/CUDA_VISIBLE_DEVICES/d' /opt/SVHunter/*.py
RUN sed -i 's/np.save(filename_data,all_data)/np.savez_compressed(filename_data,all_data)/' /opt/SVHunter/SVHunter_generate_data.py && sed -i '/def model_predict/,/def average_read_coverage/s/x_t = np.load(x_train_name1)/x_t = np.load(x_train_name1 + ".npz")["arr_0"]/' /opt/SVHunter/SVHunter_detect.py
