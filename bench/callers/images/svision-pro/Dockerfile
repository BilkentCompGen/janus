# SVision-pro 2.5: the biocontainer pins pysam 0.9.1, whose VariantRecordInfo has no .pop(), so SVision-pro 2.5 crashes
# while writing breakend records (src/output_vcf.py:194). Same image with pysam 0.21.0 (last cp37 wheel).
FROM quay.io/biocontainers/svision-pro:2.5--pyhdfd78af_0
RUN pip install --no-cache-dir pysam==0.21.0
