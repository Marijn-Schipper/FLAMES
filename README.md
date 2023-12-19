Instructions for installing FLAMES:
1. Download the github
2. Download the required annotation data from link
3. Download MAGMA
4. Download PoPS 
5. Create virtual enviroment  with required packages (recommended) or install needed packages.

How to run FLAMES:
1. Run MAGMA on your summary statistics to obtain gene-level Z-scores.
2. Run MAGMA tissue type analysis using your MAGMA Z-scores on the preformatted GTEx tissue expression file.
3. Run PoPS on the generated MAGMA z-scores. The features used in the FLAMES manuscript can be downloaded from link.
4. Format your credible sets. The required format is two tab separated columns. 
The first column should contain the SNPs in the credible set in the format CHR:BP:A1:A2 (e.g. 2:2345123:A:T).
The second column should contain the fine-mapping PIP for each SNP.
5. Run FLAMES annotate on the credible set. This generates a file with all SNP-to-gene evidence from the credible set for genes in the locus, plus MAGMA-Z and PoPS scores.
6. Run FLAMES score on the previously generated annotation file.

Running FLAMES annotate faster:
Default FLAMES will query the VEP API for the variants within your credible sets for the VEP features of interest.
This is the biggest bottleneck for annotation speed. You can significantly speed up this process by running a command line version of VEP.
Please use the --cmd-vep and --vep-cache and --vep_docker flags to use the cmd line version of VEP for FLAMES.


