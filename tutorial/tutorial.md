# Documentation on running FLAMES end-to-end
- In this document, we are going to run FLAMES on the twinning GWAS to replicate the findings in Schipper et al. 2025

## Setting up
1. Obtain PoPS
	```
	git clone https://github.com/FinucaneLab/pops.git
	```

2. Obtain PolyFun (used for finemapping)
	```
	git clone https://github.com/omerwe/polyfun.git
	```

3. Obtain FLAMES
	```
	git clone https://github.com/Marijn-Schipper/FLAMES.git
	```

4. Installation of packages
- We are going to create 2 conda environments (in progress: compling all requirements into one file)
	- one for FLAMES and PoPS (since the `requirement.txt` file from FLAMES includes the requirements in PoPS so the FLAMES conda environment should run both PoPS and FLAMES)
	- one for polyfun

- To install the flames conda enviroment: 
	```
	conda env create --file environment.yml
	conda activate FLAMES
	```

- To install the polyfun conda environment (source: https://github.com/omerwe/polyfun)
	```
	git clone https://github.com/omerwe/polyfun
	cd polyfun
	conda env create -f polyfun.yml
	conda activate polyfun
	```

- Make sure to also have tabix installed: 
	```
	wget https://github.com/samtools/htslib/releases/download/1.23/htslib-1.23.tar.bz2
	```

	or try: https://anaconda.org/channels/bioconda/packages/tabix/overview

- Make sure to install ensemble-vep if using the cache ensemble
	- I used a docker image. Example docker file: 
	```
	FROM ensemblorg/ensembl-vep:latest

	# Set working directory
	WORKDIR /data

	# Default command
	CMD ["vep", "--help"]
	```

	- build docker image
	```
	docker build -t vep .
	```

	- Alternatively, commands to install from source: 
	```
	git clone https://github.com/Ensembl/ensembl-vep.git
	cd ensembl-vep
	perl INSTALL.pl

	```


5. Download annotation data from Zenodo
- For compatibility with FUMA, you need to download `Annotation_data.tar.gz`, `pops_features_full_FUMA_compatible.tar.gz`, and `gtex_v8_ts_avg_log2TPM.txt` from: https://zenodo.org/records/12635505
- Untar:
```
tar -xzvf Annotation_data.tar.gz
tar -xzvf pops_features_full_FUMA_compatible.tar.gz
```

8. (optional) CADD and VEP
- FLAMES annotation step queries the API for CADD score and VEP. The query can be timed-out if there are too many requests. If preferred, you can utilize the cached functionality that is built in to FLAMES. To use this, you need to download the following: 
- CADD: 
```
wget https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh37/whole_genome_SNVs.tsv.gz
wget https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh37/whole_genome_SNVs.tsv.gz.tbi
```

- ensemble-vep
```
wget https://ftp.ensembl.org/pub/release-115/variation/indexed_vep_cache/homo_sapiens_vep_115_GRCh37.tar.gz
tar -xzf homo_sapiens_vep_115_GRCh37.tar.gz
```

## Obtain GWAS summary statistics
- Download the GWAS summary statistics from Mbarek et al. 2025
	- Paper: https://pubmed.ncbi.nlm.nih.gov/38052102/
	- Download GWAS sumstat from https://www.ebi.ac.uk/gwas/publications/38052102

	```
	https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90244001-GCST90245000/GCST90244755/GCST90244755_buildGRCh37.tsv
	```

- Run FUMA SNP2GENE (https://fuma.ctglab.nl/snp2gene), making sure to select MAGMA at step 6.
	- Download these files: 
		- `input.snps`
			- Please note that this file is **not** available to be downloaded from FUMA. To simplify, you can download the already bgzipped file from the `tutorial` folder. 
		- `GenomicRiskLoci.txt`
		- `magma.genes.out`
		- `magma.genes.raw`
		- `magma_exp_gtex_v8_ts_general_avg_log2TPM.gsa.out`
- Bgzip `input.snps`
	- Add the character `#` to the header
	- The header of the file `input.snps` from FUMA is `chr     bp      non_effect_allele       effect_allele   rsID    p       beta`. To successfully run polyfun for the finemapping step, update the header to `#chr    bp      A2       A1   rsID    p       beta`
	```
	bgzip input.snps
	```
	-Note that since the bgzipped `input.snps.gz` has already been prepared for you, you can **skip this step**. 

- Tabix
```
tabix -p vcf input.snps.gz
```

- In sum, make sure that you have these files in the working directory for each GWAS before proceeding
	- `input.snps.gz`
	- `input.snps.gz.tbi`
	- `GenomicRiskLoci.txt`
	- `magma.genes.out`
	- `magma.genes.raw`
	- `magma_exp_gtex_v8_ts_general_avg_log2TPM.gsa.out`


## Run PoPS

- Below is the example command, please replace with correct path before running. 

```
python {/PATH/}/pops/pops.py \
  --gene_annot_path {/PATH/}/pops_features_full_FUMA_compatible/gene_annots.txt \
  --feature_mat_prefix {/PATH/}/pops_features_full_FUMA_compatible/features_munged/pops_features \
  --num_feature_chunks 116 \
  --magma_prefix magma \
  --control_features {/PATH}/pops_features_full_FUMA_compatible/control.features \
  --out_prefix twinning_Mbareketal2024
```

- Outputs: 
	- `twinning_Mbareketal2024.preds`
	- `twinning_Mbareketal2024.coefs`
	- `twinning_Mbareketal2024.marginals`

## Run finemapping with polyfun
- Use the provided script `run_polyfun.py`: 
	```
	python run_polyfun.py --filedir {/path/to/working/directory/} --sample_size {/sample/size/}
	```
- Before running the above script, please make sure: 
	- on line 49, make sure that the path to header.txt file is correct
	- on line 64, make sure to update the path to polyfun `munge_polyfun_sumstats.py` script
	- on line 70, make sure to update the path to polyfun `munge_polyfun_sumstats.py` script

- Outputs: `locus_{chr}:{start}-{end}.susie.finemapped`
- Additional formating steps: `format_polyfun_out.py`
```
python format_polyfun_out.py locus_{chr}:{start}-{end}.susie.finemapped locus_{X}.cred1
```

## Make the index file
- Use the provided script `make_indexfile.py`
	```
	python make_indexfile.py --filedir {/path/to/working/directory/}
	```

### Run FLAMES annotate

```
python {/PATH/}/FLAMES/FLAMES.py annotate \
-o {/PATH/} \
-a {/PATH/}/Annotation_data \
-p {/PATH/}twinning_Mbareketal2024.preds \
-m {/PATH/}magma.genes.out \
-mt {/PATH/}magma_exp_gtex_v8_ts_general_avg_log2TPM.gsa.out \
-id {/PATH/}indexfile.txt
```

- To use the cache CADD and VEP, use the following command: 
```
python {/PATH/}/FLAMES/FLAMES.py annotate \
-o {/PATH/} \
-a {/PATH/}/Annotation_data \
-p {/PATH/}twinning_Mbareketal2024.preds \
-m {/PATH/}magma.genes.out \
-mt {/PATH/}magma_exp_gtex_v8_ts_general_avg_log2TPM.gsa.out \
-id {/PATH/}indexfile.txt \            
-t tabix \ #tabix executable \
-cf {/PATH/}/whole_genome_SNVs.tsv.gz" \
-cv vep \ #vep command line
-vc /.vep" #path to vep cache
```

### Step 5: Run FLAMES
```
python {PATH}/FLAMES/FLAMES.py FLAMES \
-id indexfile.txt \
-o {/PATH/}
```


