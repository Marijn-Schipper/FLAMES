FLAMES version 1.0.0.
Thank you for your interest in using FLAMES for GWAS gene prioritization.
The Python version of FLAMES is still being optimized. 
If you have any problem using/installing FLAMES please open an issue.

# Installation
1. Download the FLAMES from this GitHub
2. Download the required annotation data from [Zenodo](https://zenodo.org/records/10417198) 
3. Create virtual enviroment  with required packages (recommended) or install needed packages.
You can do this either by using conda by running:  
       1: ```conda env create --file environment.yml```  
       2: ```conda activate FLAMES```   
      When using Anaconda, you can select import and navigate to the environment.yml file contained in the downloaded FLAMES folder in git and run python within that environment.  
Or by running:  
       ```pip install -r requirements.txt```  
   We highly recommend installing miniconda and creating an environment within conda as this provides a stable environment for FLAMES to run.  
   

This is everything you need to run the tutorial
If you want to run FLAMES on your own GWAS results you will also need the following:
MAGMA
PoPS
Finemapping results from statistical fine-mapping software e.g. FINEMAP or SusieR

# Running FLAMES
## Running FLAMES on example data:
To run FLAMES on the provided example data navigate to the example_data folder in the downloaded FLAMES folder from this GitHub.
Make sure you are in an environment that contains the dependencies needed for FLAMES and that you have downloaded the reference data from Zenodo as described in the Installation section.

### Running annotation on example data
To run FLAMES on the provided example data run the following commands:
Change to the example data directory:  
```cd example_data```  
Run flames annotate. This should take around 1 minute.
```
python {PATH_TO_FLAMES_FOLDER}/FLAMES.py annotate \
-a {PATH_TO_DOWNLOADED_ANNOTATION_DATA}/Annotation_data/ \
-p PoPS.preds \
-m magma.genes.out \
-mt magma_exp_gtex_v8_ts_avg_log2TPM.txt.gsa.out \
-id indexfile.txt  \
-pc prob1 \
-sc cred1 \
-g genes.txt \
-c95 False
```

This will run FLAMES annotate on the four loci of twinning as described in the FLAMES paper.

### FLAMES scoring on created annotated loci
python {PATH_TO_FLAMES_FOLDER}/FLAMES.py FLAMES -id indexfile.txt -o ./

This will score the generated annotated loci and produce the results reported for twinning in the FLAMES paper.

## Running FLAMES on your data from scratch
Step 1 and 2 can be performed by uploading your summary statistics to FUMA and running MAGMA there, and downloading the final results.
### 1. Run MAGMA on your summary statistics to obtain gene-level Z-scores. 
You can find information on how to do this on the [MAGMA website](https://ctg.cncr.nl/software/magma).
in general your command to run MAGMA will look like:
```
./magma \
 --bfile {PATH_TO_REFERENCE_PANEL_PLINK} \
 --gene-annot {PATH_TO_MAGMA_ANNOT}.genes.annot \
 --pval {PATH_TO_SUMSTATS}.txt ncol=N \
 --gene-model snp-wise=mean \
 --out {DESIRED_ZSCORE_FILENAME}
```
   
### 2. Run MAGMA tissue type analysis using your MAGMA Z-scores on the preformatted GTEx tissue expression file. 
The command uses the previously generated MAGMA gene Z-scores and GTEx expression file which can be found on [Zenodo](https://zenodo.org/records/10409723)
```
./magma \
--gene-results {DESIRED_ZSCORE_FILENAME}.genes.raw \
--gene-covar {PATH_TO_DOWNLOADED_GTEx_FILE}/gtex_v8_ts_avg_log2TPM.txt \
--out {DESIRED_TISSUE_RELEVANCE_FILENAME}
```

### 3. Run PoPS on the generated MAGMA z-scores. 
The features used in the FLAMES manuscript can be downloaded from [Zenodo](https://zenodo.org/records/10409723). You can find the github for PoPS [here](https://github.com/FinucaneLab/pops).
```
python pops.py \
--gene_annot_path {PATH_TO_DOWNLOADED_FEATURES}\pops_features_pathway_naive/gene_annot.txt \
--feature_mat_prefix {PATH_TO_DOWNLOADED_FEATURES}\pops_features_pathway_naive/munged_features/pops_features \
--num_feature_chunks 99 \
--magma_prefix {PATH_TO_GENERATED_MAGMA_Z_SCORES}\{DESIRED_ZSCORE_FILENAME} \
--control_features {PATH_TO_DOWNLOADED_FEATURES}\pops_features_pathway_naive/control.features \
--out_prefix {DESIRED_POPS_OUTPUT_PREFIX)
```
   
### 4. Format your credible sets. 
The required format is two tab separated columns. 
The first column should contain the SNPs in the credible set in the format CHR:BP:A1:A2 (e.g. 2:2345123:A:T).
The second column should contain the fine-mapping PIP for each SNP.

### 5. Run FLAMES annotate on the credible set. 
This generates a file with all SNP-to-gene evidence from the credible set for genes in the locus, plus MAGMA-Z and PoPS scores. An example command could be:
```
python FLAMES.py annotate \
-o {DESIRED_OUTPUT_DIRECTORY} \
-a {PATH_TO_THE_DOWNLOADED_ANNOTATION_DATA_DIRECTORY} \
-p {DESIRED_POPS_OUTPUT_PREFIX}.preds \
-m {DESIRED_ZSCORE_FILENAME}.genes.out \
-mt {DESIRED_TISSUE_RELEVANCE_FILENAME}.gsa.out \
-id {PATH_TO_INDEXFILE} 
```

The INDEXFILE should contain the following column including the header:
- Filename : The path to the formatted credible set file

To predefine the output names of all the generated credible sets, also create the column Annotfiles.
- Annotfiles : path and outputname for credible set. (e.g. /home/annotated_loci/annotated_locus_1.txt)

To run with predefined locus definitions add the -l flag with a the GENOMIC_LOCI_FILE. This file should contain the folowing tab separated columns including headers:
- GenomicLocus : a unique identifier of a locus
- chr : the chromosome of the locus
- start : start of the locus location in bp
- end : the end of the locus location in bp
The GenomicLocus columns should now also be added to the INDEXFILE so that the credset matches the correct locus.
   The INDEXFILE should then contain:
   Filename : The path to the formatted credible set file
   GenomicLocus : the unique identifier that matches the GENOMIC_LOCI_FILE
   
### 6. Run FLAMES scoring on the previously generated annotation file. 
For running FLAMES scoring you should run FLAMES.py FLAMES with the following options:
-id or -i where :id points to a tab delimited txt file containing the column Annotfiles with the output from FLAMES annotate that you want to score. i will point to a txt file containing all the desired input files. Each file should be on their own line.
-o: the desired output directory or filename.
The command will look something like:
```
python FLAMES.py FLAMES \
-id {INDEX_FILE_NCLUDING_COLUMN Annotfiles} \
-o {DESIRED_OUTPUT_DIRECTORY}
```


Running FLAMES annotate faster:
Default FLAMES will query the VEP API for the variants within your credible sets for the VEP features of interest.
This is the biggest bottleneck for annotation speed. You can significantly speed up this process by running a command line version of VEP.
Please use the --cmd-vep and --vep-cache and --vep_docker flags to use the cmd line version of VEP for FLAMES.


