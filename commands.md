# Simulation

I followed akmami's tutorial but with few changes along the path. Instead of 20 sublines this creates 10 sublines.

https://github.com/akmami/docs/blob/main/PSITE.md

## MS Tool
Instead of installing ms, you can use msprime which offers the `mspms` command.

### Installation
To install msprime, use the following command:

```bash
pip install msprime
```

### Example Usage
You can run the following command to generate a tree and save it to a file:

```bash
mspms 10 1 -T -G 1 | tail -n1 > ms_tree.txt
```

This command simulates a coalescent tree with 10 samples and saves the output to `ms_tree.txt`.

## PSITE Usage (for me)

Note: cfg_template_female.yaml file was in hg19 from the original repository. I converted the chr lengths to match them to hg38

```bash
cd tool/PSiTE/
python3 psite.py vcf2fa -r /home/yuio/playground/data/processed/hg38/human_v38_no_chr.fasta -v ../../data/ref/HG001_GRCh38_1_22_v4.2.1_all.phased.vcf.gz --autosomes 1..22 --sex_chr X,X -o ../../data/psite/normal_fa
python3 psite.py phylovar -t ../../data/ref/ms_tree.txt --config ../../data/psite/cfg_template_female.yaml --purity 0.8 --sex_chr X,X --trunk_length 2.0 --prune 0.05 --chain tumor_chain --map map
python3 psite.py chain2fa -c tumor_chain -n ../../data/psite/normal_fa/normal.parental_0.fa,../../data/psite/normal_fa/normal.parental_1.fa -o ../../data/psite/tumor_fa --cores 16
```

## Merge Diploid Files

```bash
cd data/psite
nodes=$(ls ./tumor_fa | grep -oP 'node\d+' | sort | uniq); for node in $nodes; do parental_0="./tumor_fa/${node}.parental_0.fa"; parental_1="./tumor_fa/${node}.parental_1.fa"; if [[ -f $parental_0 && -f $parental_1 ]]; then cat $parental_0 $parental_1 > ./tumor_fa/${node}.fa; echo "Merged $parental_0 and $parental_1 into ./tumor_fa/${node}.fa"; fi done
```

## Move simulation data

```bash
mkdir -p ./simulation/tumors/
mkdir -p ./simulation/normal/
mv ./data/psite/tumor_fa/node{11,12,15,16,18,19,3,5,7,8}.fa ./simulation/tumors/
mv ./data/psite/normal_fa/normal.parental_0.fa ./simulation/normal/normal.p0.fa
mv ./data/psite/normal_fa/normal.parental_1.fa ./simulation/normal/normal.p1.fa
```

## Cleanup

```bash
rm ./data/psite/tumor_fa -rf
rm ./data/psite/normal_fa -rf
rm ./tool/PSiTE/{tumor_chain,phylovar_cnvs,phylovar_snvs,map} -rf
```

# LCP Tools

```bash
git clone https://github.com/BilkentCompGen/lcptools.git
cd lcptools

sudo make install
```

Then I needed to change the line 81 in lcptools.c from

```cpp
LCP_INIT(0);
```

to

```cpp
LCP_INIT2(0);
```

I also needed to run these for library to work.

```bash
echo "/usr/local/lib" | sudo tee /etc/ld.so.conf.d/local.conf
sudo ldconfig
```