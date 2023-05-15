#import modules
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
import os
import multiprocessing as mp
from pyliftover import LiftOver
import numpy as np
from tqdm import tqdm
from sys import argv
import sys
import time
import Query_api


###define functions helper functions
"""Function to time a process"""
def substep_timer(timestamp):
    start_time = timestamp
    end_time = time.time()
    print(f"Execution time: {end_time - start_time} seconds")
    return end_time

"""Function to merge dictionaries"""
def merge_dicts(*dicts):
    d = {}
    for dict in dicts:
        for key in dict.keys():
            val = dict[key][0] if isinstance(dict[key], list) else dict[key]
            try:
                d[key].append(val)
            except KeyError:
                d[key] = [val]
    return d

def validate_build(build):
    if build.upper() == 'HG19' or build.upper() == 'GRCH37':
        build = 'GRCH37'
    elif build.upper() == 'HG38' or build.upper() == 'GRCH38':
        build = 'GRCH38'
    else:
        print('Invalid build. Please use build hg19 or hg38')
        sys.exit()
    return build

"""Function for checking if a genomic location is within a range"""
def check_range(start, stop, loc):
    return start <= loc <= stop

"""Function that will print the runtime of a function"""
def runtime(func):
    def wrapper(*args, **kwargs):
        start = time.time()
        func(*args, **kwargs)
        end = time.time()
        print(f'Runtime: {end - start}')
    return wrapper

"""Annotate a credible set based on variant position within annotation range"""	
def create_positional_mapping_matrix(df1, df2, colnames=['start', 'end'], chr=False):
    if not chr == False:
        if chr == 'X':
            chr = 23
        df1 = df1[df1['chr'] == chr]
    if 'pos_37' in df2.columns:
        poscol = 'pos_37'
    else:
        poscol = 'pos'
    #reshape and perform matrix factorization to check whether a SNP is within a genomic range
    starts = df1[colnames[0]].values.reshape(-1, 1)
    ends = df1[colnames[1]].values.reshape(-1, 1)
    positions = df2[poscol].values.reshape(1, -1)
    matches = np.logical_and(starts <= positions, ends >= positions)
    return matches

"""Liftover genomic coordinates column from one genome build to another into a new column"""
def liftover_df(lo, df, incols, outcol):
    df['new'] = df.apply(lambda row: lo.convert_coordinate(f"chr{row[incols[0]]}", int(row[incols[1]]))[0][1], axis=1)
    return df

"""For any given locus, retrieve the genes present in that locus  """
def get_genes_in_locus(locus_no, path_to_genes_file, GenomicRiskLoci):
    genes = pd.read_csv(path_to_genes_file, sep='\t')
    genes['GenomicLocus'] = genes.GenomicLocus.str.split(":")
    genes = genes.explode('GenomicLocus')
    relevant_genes = pd.DataFrame()
    locusname = str(GenomicRiskLoci.iloc[int(locus_no)-1][0])
    for locus in locusname.split(':'):
        genes_in_loc = genes[genes['GenomicLocus'] == locus]
        relevant_genes = pd.concat([relevant_genes, genes_in_loc])
    return relevant_genes

"""For any given locus, retrieve the genes present in that locus  """
def Genes_in_locus_bp_based(locus_no, GenomicRiskLoci, baseline_genes, lo=None):
    start = GenomicRiskLoci.iloc[int(locus_no)-1]['start']
    end = GenomicRiskLoci.iloc[int(locus_no)-1]['end']
    chrom = GenomicRiskLoci.iloc[int(locus_no)-1]['chr']
    if chrom == 'X':
        chrom = '23'
    if not lo == None:
        start = lo.convert_coordinate(f"chr{chrom}", int(start))[0][1]
        end = lo.convert_coordinate(f"chr{chrom}", int(end))[0][1]
    genes = pd.read_csv(baseline_genes, sep='\t', engine='pyarrow')
    genes['chromosome_name'] = genes['chromosome_name'].replace(['X', 'x'], '23')

    genes = genes[(genes['end_position'] >= start) & (genes['start_position'] <= end) & (genes['chromosome_name'].astype(int) == int(chrom))]
    genes = genes.rename(columns={'chromosome_name': 'chr', 'start_position': 'start', 'end_position': 'end', 'ensembl_gene_id': 'ensg', 'external_gene_name': 'symbol'})
    return genes

"""Check if genes in locus are in a list of true postive genes for creating a TP set"""
def check_TP_genes(genes_in_locus, TP_genes):
    with open(TP_genes) as file:
        TP_genes = [line.rstrip() for line in file]
    ## check if genes in locus symbol or ensgid are in TP_genes and add column to df
    genes_in_locus = genes_in_locus[genes_in_locus['ensg'].isin(TP_genes) | genes_in_locus['symbol'].isin(TP_genes)]
    return genes_in_locus


"""For the genes in a given finemapped locus, determine the distance between the TSS and gene body 
    for the SNP with highest PiP and for PiP weighted centroid """
def get_TSS_distances(creds, genes_in_locus, build):
    distances = []
    weighted_distances = []
    TSS_distances =[]
    weighted_TSS_distances = []
    poscol = 'pos_37' if 'pos_37' in creds.columns else 'pos'
    centroid = sum(creds[poscol] * creds['prob1']) / sum(creds['prob1'])
    most_likely_snp = creds.loc[creds['prob1'].idxmax()][poscol]
    print('Calculate distance metrics for each gene')
    with tqdm(total=genes_in_locus.shape[0]) as pbar:
        for index, row in genes_in_locus.iterrows():
            if row['strand'] == '+':
                TSS_distances.append(abs(row['start'] - most_likely_snp))
                weighted_TSS_distances.append(abs(row['start'] - centroid))
            else:
                TSS_distances.append(abs(most_likely_snp - row['end']))
                weighted_TSS_distances.append(abs(centroid - row['end']))
            if row['start'] < centroid < row['end']:
                distances.append(0)
                weighted_distances.append(0)
            else:
                distances.append(min(abs(row['start'] - most_likely_snp), abs(row['end'] - most_likely_snp)))
                weighted_distances.append(min(abs(row['start'] - centroid), abs(row['end'] - centroid)))
            pbar.update(1)
    genes_in_locus['TSS_distance'] = TSS_distances
    genes_in_locus['weighted_TSS_distance'] = weighted_TSS_distances
    genes_in_locus['distance'] = distances
    genes_in_locus['weighted_distance'] = weighted_distances
    return(genes_in_locus)

"""Create tissue type scaling for eQTL based on MAGMA relevance"""
def create_relevance_dict(magma_scores):
    magma_scores['score'] = magma_scores['P'].rank(ascending=False) / len(magma_scores['P'])
    threshold = 0.05/len(magma_scores['score'])
    magma_scores['score'] = magma_scores.apply(lambda row: row['score'] * 0.5 if row['P'] > 0.05 else row['score'] * 2 if row['P'] <= threshold else row['score'], axis=1)
    magma_scores_dict = dict(zip(magma_scores['VARIABLE'], magma_scores['score']))
    return magma_scores_dict

""""For each variant in credible set, get the finemapped GTEx eqtls and score each gene in locus"""
def get_GTEX_eQTLs(creds, genes, finemapped_gtex, build, magma_scores=False):
    #prep data for input build
    if build.upper() == 'GRCH38':
        creds['eqtl_id'] = creds.apply(lambda row: f"chr{row['chr']}_{row['pos']}_{row['a1']}_{row['a2']}_b38", axis=1)
        creds['eqtl_id2'] = creds.apply(lambda row: f"chr{row['chr']}_{row['pos']}_{row['a2']}_{row['a1']}_b38", axis=1)
        header = 'eQTL'
    else:
        creds['eqtl_id'] = creds.apply(lambda row: f"{row['chr']}_{row['pos']}_{row['a1']}_{row['a2']}_b37", axis=1)
        creds['eqtl_id2'] = creds.apply(lambda row: f"{row['chr']}_{row['pos']}_{row['a2']}_{row['a1']}_b37", axis=1)
        header = 'eQTL_hg37'
    #GTEX get CLPP
    CLPP_sum = []
    CLPP_max = []
    tissue_specific_CLPP_sum = []
    tissue_specific_CLPP_max = []
    print('Get GTEx CLPP for each gene')
    with tqdm(total=genes.shape[0]) as pbar:
        for index, row in genes.iterrows():
            pbar.update(1)
            maxval, sumval = 0,0
            genepath = os.path.join(finemapped_gtex, row['ensg'].split('.')[0] + '.txt')
            if not os.path.isfile(genepath):
                CLPP_sum.append(sumval)
                CLPP_max.append(maxval)
                tissue_specific_CLPP_sum.append(sumval)
                tissue_specific_CLPP_max.append(maxval)
                continue
            eqtls = pd.read_csv(os.path.join(genepath), sep='\t',engine='pyarrow')
            overlap = pd.merge(creds, eqtls, left_on='eqtl_id', right_on=header, how='inner')
            overlap = pd.concat([overlap, pd.merge(creds, eqtls, left_on='eqtl_id2', right_on=header, how='inner')])
            if overlap.shape[0] > 0:
                sumval = sum(list(overlap['Probability'] * overlap['prob1']))
                maxval = max(list(overlap['Probability'] * overlap['prob1']))
            CLPP_sum.append(sumval)
            CLPP_max.append(maxval)
            if type(magma_scores) == dict and overlap.shape[0] > 0:
                overlap['scores'] = overlap['TISSUE'].map(magma_scores)
                overlap['Probability'] = overlap['Probability'] * overlap['scores']
                sumval = sum(list(overlap['Probability'] * overlap['prob1']))
                maxval = max(list(overlap['Probability'] * overlap['prob1']))
                tissue_specific_CLPP_sum.append(sumval)
                tissue_specific_CLPP_max.append(maxval)
            else:
                tissue_specific_CLPP_sum.append(sumval)
                tissue_specific_CLPP_max.append(maxval)
    genes['CLPP_GTEX_eQTL_sum'] = CLPP_sum
    genes['CLPP_GTEX_eQTL_max'] = CLPP_max
    if type(magma_scores) == dict:
        genes['CLPP_GTEX_tissue_weighted_eQTL_sum'] = tissue_specific_CLPP_sum
        genes['CLPP_GTEX_tissue_weighted_eQTL_max'] = tissue_specific_CLPP_max
    return genes

""""For each variant in credible set, get the finemapped eQTLgen eqtls and score each gene in locus"""
def get_eqtlgen_eQTLs(creds, genes, eQTL_dir , build):
    chrom = creds.iloc[0]['chr']
    if build.upper == 'GRCH38':
        lo = LiftOver('hg38', 'hg19')
        creds = liftover_df(lo, creds, ['chr', 'pos'], ['pos'])
    CLPP_max = []
    CLPP_sum = []
    #for each gene in locus:
    print('Get eQTLGen CLPP for each gene')
    with tqdm(total=genes.shape[0]) as pbar:    
        for index, row in genes.iterrows():
            maxval, sumval = 0, 0
            pbar.update(1)
            fname = os.path.join(f"{eQTL_dir}", f"chr{chrom}", "post" ,f"{row['ensg'].upper()}_new.out_post")
            # if that gene has finemapped eQTLgen variants
            if os.path.isfile(fname):
                posteriors = pd.read_csv(fname, sep="\t", engine='pyarrow')
                posteriors['SNP_ID'] = posteriors['SNP_ID'].astype(np.int32)
                #merge finemapped eQTLgen variants with locus finemapping variants
                gene_CLPP = pd.merge(creds, posteriors, left_on=['pos'], right_on=['SNP_ID'])
                #if there are no finemapped eqtlgen variants which are also finemapped variants in locus then CLPP = 0
                if gene_CLPP.shape[0] > 0:
                    maxval = max(list(gene_CLPP['Causal_Post._Prob.'] * creds['prob1']))
                    sumval = sum(list(gene_CLPP['Causal_Post._Prob.'] * creds['prob1']))
            CLPP_max.append(maxval)
            CLPP_sum.append(sumval)
    genes['CLPP_eQTLgen_sum'] = CLPP_sum
    genes['CLPP_eQTLgen_max'] = CLPP_max
    return genes

""""For each variant in credible set, get the VEP annotation score each gene in locus"""
def get_VEP(creds, genes, build):
    print('Query and assign VEP scores for each variant')
    with tqdm(total=creds.shape[0]) as pbar:
        VEP_dict = {'HIGH':1, 'MODERATE':0.66, 'LOW':0.33, 'MODIFIER':0}
        VEPs = creds.apply(lambda row: Query_api.query_VEP(row['chr'], row['pos'], row['a1'], row['a2'], build), axis=1)
        scores ={}
        for i in range(len(VEPs)):
            variant_consequences ={}
            variants = VEPs[i][0]
            PiP = creds['prob1'][i]
            for variant in variants['transcript_consequences']:
                if 'gene_id'in variant.keys():
                    max_consequence = max(variant_consequences.get(variant['gene_id'], float(0)), float(PiP) * float(VEP_dict[variant['impact']]))
                    variant_consequences[variant['gene_id']] = max_consequence
            scores = merge_dicts(scores, variant_consequences)
            pbar.update(1)
    genes[['VEP_max', 'VEP_sum']] = genes['ensg'].apply(lambda row: pd.Series([max(scores.get(row, [0])), sum(scores.get(row, [0]))]))
    return genes

""""For each variant in credible set, get the CADD annotation score each gene in locus"""
def get_CADD(creds, genes, build):
    print('Query and assign CADD scores for each variant')
    with tqdm(total=1) as pbar:
        pbar.update(1)
        creds['CADD'] = creds.apply(lambda row: Query_api.query_CADD(row['chr'], row['pos'], row['a1'], row['a2'], build), axis=1)
        matrix = create_positional_mapping_matrix(genes, creds)
        creds['CADD_weighted'] = creds['CADD'] * creds['prob1']
        matrix2 = pd.DataFrame(np.outer(np.ones(len(genes)), creds['CADD_weighted']))
        genes['CADD_sum'] = np.dot(matrix, creds['CADD'])
        genes['CADD_max'] = np.max(matrix2, axis=1)
    return genes

""""For each gene annotate association scores from Jung PCHiC with credible set variants"""
def get_jung_HiC(jung, genes, creds):
    sums = []
    maxs = []
    print('Get Jung HiC interactions for each gene')
    with tqdm(total=genes.shape[0]) as pbar:
        for gene in genes['ensg']:
            pbar.update(1)
            if not os.path.exists(os.path.join(jung, gene + '.txt')):
                sums.append(0)
                maxs.append(0)
                continue
            jungs = pd.read_csv(os.path.join(jung, gene + '.txt'), sep=" ", engine='pyarrow')
            ## split jungs 'raw' on '.' and split into three columns called chr, start and end
            jungs[['tmp', 'start', 'end']] = jungs['target'].str.split('.', expand=True)
            jungs['p'] = jungs['dist_pvalue']
            jungs['start'] = jungs['start'].astype(int)
            jungs['end'] = jungs['end'].astype(int)
            ## create a matrix that maps which variant is in which captured region
            matrix = create_positional_mapping_matrix(jungs, creds)
            matrix2 = pd.DataFrame(np.outer(jungs['p'], creds['prob1']))
            matrix2 = matrix2.where(matrix, other=0)
            sums.append(matrix2.values.sum())
            maxs.append(matrix2.values.max())
    genes['Jung_HiC_sum'] = sums
    genes['Jung_HiC_max'] = maxs      
    return genes

""""For each gene annotate association scores from Javierre PCHiC with credible set variants"""
def get_javierre_HiC(javierre, genes, creds):
    sums = []
    maxs = []
    print('Get Javierre HiC interactions for each gene')
    with tqdm(total=genes.shape[0]) as pbar:    
        for gene in genes['ensg']:
            pbar.update(1)
            if not os.path.exists(os.path.join(javierre, gene + '.txt')):
                sums.append(0.0)
                maxs.append(0.0)
                continue
            javierres = pd.read_csv(os.path.join(javierre, gene + '.txt'), sep="\t", engine='pyarrow')
            celltypes = [ 'Mon', 'Mac0', 'Mac1','Mac2', 'Neu', 'MK', 'EP', 'Ery', 'FoeT', 'nCD4', 'tCD4', 'aCD4','naCD4', 'nCD8', 'tCD8', 'nB', 'tB']
            ## create a max and sum column in javierres based on the celltypes values
            javierres['max'] = javierres[celltypes].max(axis=1)
            javierres['sum'] = javierres[celltypes].sum(axis=1)
            # melt the dataframe to have each celltype measurement into individual rows
            matrix = create_positional_mapping_matrix(javierres, creds, ['oeStart', 'oeEnd'])
            matrix2 = pd.DataFrame(np.outer(javierres['sum'], creds['prob1'])).where(matrix, other=0)
            matrix3 = pd.DataFrame(np.outer(javierres['max'], creds['prob1'])).where(matrix, other=0)
            sums.append(matrix2.values.sum())
            maxs.append(matrix3.values.max())
        genes['Jung_HiC_sum'] = sums
        genes['Jung_HiC_max'] = maxs      
    return genes

def get_ABC_EP(ABC_EP, genes, creds):
    print('Get ABC Enhancer Promotor interactions for each gene')
    sums = []
    maxs = []
    with tqdm(total=genes.shape[0]) as pbar:
        for gene in genes['ensg']:
            pbar.update(1)
            if not os.path.exists(os.path.join(ABC_EP, gene + '.txt')):
                sums.append(0.0)
                maxs.append(0.0)
                continue
            ABCs = pd.read_csv(os.path.join(ABC_EP, gene + '.txt'), sep="\t", engine='pyarrow')
            matrix = create_positional_mapping_matrix(ABCs, creds)
            matrix2 = pd.DataFrame(np.outer(ABCs['ABC.Score'], creds['prob1'])).where(matrix, other=0)
            sums.append(matrix2.values.sum())
            maxs.append(matrix2.values.max())
    genes['ABC_EP_sum'] = sums
    genes['ABC_EP_max'] = maxs
    return genes

def get_ABC_CRISPR(ABC_CRISPR, genes, creds):
    print('Get ABC CRISPR interactions for each gene')
    sums = []
    maxs = []
    with tqdm(total=genes.shape[0]) as pbar:
        for gene in genes['ensg']:
            pbar.update(1)
            if not os.path.exists(os.path.join(ABC_CRISPR, gene + '.txt')):
                sums.append(0.0)
                maxs.append(0.0)
                continue
            ABCs = pd.read_csv(os.path.join(ABC_CRISPR, gene + '.txt'), sep="\t", engine='pyarrow')
            matrix = create_positional_mapping_matrix(ABCs, creds)
            matrix2 = pd.DataFrame(np.outer(ABCs['ABC_max'], creds['prob1'])).where(matrix, other=0)
            sums.append(matrix2.values.sum())
            maxs.append(matrix2.values.max())
    genes['ABC_CRISPR_sum'] = sums
    genes['ABC_CRISPR_max'] = maxs
    return genes

def get_F5(F5_ep, genes, creds):
    print('Get Fanthom5 promoter enhancer correlation for each finemapped variant in enhancer to promoter of gene in locus')
    sums = []
    maxs = []
    with tqdm(total=genes.shape[0]) as pbar:
        for gene in genes['ensg']:
            pbar.update(1)
            if not os.path.exists(os.path.join(F5_ep, gene + '.txt')):
                sums.append(0.0)
                maxs.append(0.0)
                continue
            F5s = pd.read_csv(os.path.join(F5_ep, gene + '.txt'), sep=" ", engine='pyarrow')
            F5s[['chr', 'start', 'end']] = F5s['enhancer'].str.split(':|-', expand=True)
            F5s['start'] = F5s['start'].astype(int)
            F5s['end'] = F5s['end'].astype(int)
            matrix = create_positional_mapping_matrix(F5s, creds)
            matrix2 = pd.DataFrame(np.outer(F5s['correlation'], creds['prob1'])).where(matrix, other=0)
            sums.append(matrix2.values.sum())
            maxs.append(matrix2.values.max())
    genes['F5_ep_sum'] = sums
    genes['F5_ep_max'] = maxs
    return genes

def get_Genehancer(Epi, genes, creds):
    print('Get Genehancer interactions for each gene')
    sums = []
    maxs = []
    with tqdm(total=genes.shape[0]) as pbar:
        for gene in genes['ensg']:
            pbar.update(1)
            if not os.path.exists(os.path.join(Epi, gene + '.txt')):
                sums.append(0.0)
                maxs.append(0.0)
                continue
            EpiScores = pd.read_csv(os.path.join(Epi, gene + '.txt'), sep="\t", engine='pyarrow')
            matrix = create_positional_mapping_matrix(EpiScores, creds)
            matrix2 = pd.DataFrame(np.outer(EpiScores['score'], creds['prob1'])).where(matrix, other=0)
            sums.append(matrix2.values.sum())
            maxs.append(matrix2.values.max())
    genes['GeneHancer_max'] = maxs
    genes['GeneHancer_sum'] = sums
    return genes

def get_RoadmapEpi(Roadmap, genes, creds):
    print('Get Roadmap Epigenomics interactions for each gene')
    sums = []
    maxs = []
    with tqdm(total=genes.shape[0]) as pbar:
        for gene in genes['ensg']:
            pbar.update(1)
            if not os.path.exists(os.path.join(Roadmap, gene + '.txt')):
                sums.append(0.0)
                maxs.append(0.0)
                continue
            Epi = pd.read_csv(os.path.join(Roadmap, gene + '.txt'), sep="\t", engine='pyarrow')
            matrix = create_positional_mapping_matrix(Epi, creds)
            matrix2 = pd.DataFrame(np.outer(Epi['max'], creds['prob1'])).where(matrix, other=0)
            sums.append(matrix2.values.sum())
            maxs.append(matrix2.values.max())
    genes['RoadmapEpi_sum'] = sums
    genes['RoadmapEpi_max'] = maxs
    return genes

def get_EpiMap(EpiMap, genes, creds):
    print('Get Roadmap Epigenomics interactions for each gene')
    sums = []
    maxs = []
    with tqdm(total=genes.shape[0]) as pbar:
        for gene in genes['ensg']:
            pbar.update(1)
            if not os.path.exists(os.path.join(EpiMap, gene + '.txt')):
                sums.append(0.0)
                maxs.append(0.0)
                continue
            Epi = pd.read_csv(os.path.join(EpiMap, gene + '.txt'), sep="\t", engine='pyarrow')
            matrix = create_positional_mapping_matrix(Epi, creds)
            matrix2 = pd.DataFrame(np.outer(Epi['corr'], creds['prob1'])).where(matrix, other=0)
            sums.append(matrix2.values.sum())
            maxs.append(matrix2.values.max())
    genes['EpiMap_sums'] = sums
    genes['EpiMap_max'] = maxs
    return genes

def annotate_credset(genes,creds, build, Ann_path, magma=False, trainset=False):
    genes = get_TSS_distances(creds, genes, build)
    if trainset:
        check = check_dist_to_causal(genes)
        if not check:
            exit()
    genes = get_VEP(creds, genes, build)
    genes = get_CADD(creds,genes, build)
    genes = get_eqtlgen_eQTLs(creds, genes, os.path.join(Ann_path, 'eQTLGen', 'caviar_farhad'), build)
    genes = get_GTEX_eQTLs(creds, genes, os.path.join(Ann_path, 'GTEx_v8_finemapping_CAVIAR'), build, magma)
    genes = get_jung_HiC(os.path.join(Ann_path, 'Jung_PCHiC'), genes, creds)
    genes = get_javierre_HiC(os.path.join(Ann_path, 'Javierre_PCHiC'), genes, creds)
    genes = get_ABC_EP(os.path.join(Ann_path, 'ABC_EP'), genes, creds)
    genes = get_ABC_CRISPR(os.path.join(Ann_path, 'ABC_CRISPR'), genes, creds)
    genes = get_F5(os.path.join(Ann_path, 'EP_correlation_F5'), genes, creds)
    genes = get_Genehancer(os.path.join(Ann_path, 'GeneHancer'), genes, creds)
    genes = get_RoadmapEpi(os.path.join(Ann_path, 'RoadMapEpi'), genes, creds)
    genes = get_EpiMap(os.path.join(Ann_path, 'EpiMap'), genes, creds)
    return genes


def POPS_MAGMA_annotation(genes, magma, PoPS):
    print('Get MAGMA gene Z-scores')
    magma = magma[['GENE', 'ZSTAT']]
    genes = pd.merge(genes, magma, left_on='ensg', right_on='GENE', how='left')
    genes = genes.drop(columns=['GENE'])
    genes = genes.rename(columns={'ZSTAT': 'MAGMA_Z'}).fillna(0)
    print('Get PoPS scores')
    PoPS = PoPS[['ENSGID', 'PoPS_Score']]
    genes = pd.merge(genes, PoPS, left_on='ensg', right_on='ENSGID', how='left')
    genes = genes.drop(columns=['ENSGID'])
    return genes

def full_annotation_of_credset(path_to_credset, build,  GenomicRiskLoci, Ann_path, outdir, lo, magma, PoPS, genes, magma_scores=False, trainset=False):
    #load credset
    locno = path_to_credset[1]
    path_to_credset = path_to_credset[0]
    creds = pd.read_csv(path_to_credset, sep=' ', comment='#')
    creds[['chr', 'pos', 'alleles']] = creds['cred1'].str.split(':', expand=True)
    creds[['chr', 'pos']] = creds[['chr', 'pos']].astype(np.int32)
    creds[['a1', 'a2']] = creds['alleles'].str.split('_', expand=True)
    
    # liftover variants if build is GRCH38
    
    if build == 'GRCH38':
        creds = liftover_df(lo, creds, ['chr', 'pos'], ['pos_37'])

    # get genes in locus
    genes = get_genes_in_locus(locno, '/home/schipper/ML/annotation_test/genes.txt', GenomicRiskLoci)
    genes = Genes_in_locus_bp_based(locno, GenomicRiskLoci, Ann_path + '/ENSG/ENSG.v102.genes.txt')
    genes = annotate_credset(genes, creds, build, Ann_path, magma=magma_scores)
    genes = POPS_MAGMA_annotation(genes, magma, PoPS)
    if trainset:
        causal = (check_TP_genes(genes, trainset))

    genes.to_csv(os.path.join(outdir, 'annotated_credset_' + str(locno) + '.txt'), sep='\t', index=False)
    return

def check_dist_to_causal(genes, trainset, dist=500000):
    causal = check_TP_genes(genes, trainset)
    if causal['dist'].min() > dist:
        print('WARNING: Credible set is too far away from causal gene')
        return False
    return True



@runtime
def main(credset_dir, Annotation_dir, build, pops_out, MAGMA_outdir, outdir, loci, genes,  TP = False, filter=False , processors=False):
    #Path for annotation files
    Ann_path = os.path.normpath(Annotation_dir)

    #Load GenomicRiskLoci
    GenomicRiskLoci = pd.read_csv(os.path.normpath(loci), sep='\t', engine='pyarrow')

    #Load PoPS & MAGMA scores
    PoPS = pd.read_csv(pops_out, sep='\t', engine='pyarrow')
    magma = pd.read_csv(MAGMA_outdir + 'magma.genes.out', sep='\t', engine='pyarrow')
    magma_tissues = pd.read_csv(MAGMA_outdir + 'magma_exp_gtex_v8_ts_avg_log2TPM.txt.gsa.out', delim_whitespace=True, comment='#')
    magma_scores = create_relevance_dict(magma_tissues)

    #Check build
    build = validate_build(build)
    if build == 'GRCH38':
        lo = LiftOver('hg38', 'hg19')
    else:
        lo = False
    
    #Load credset file locations and create arguments for multiprocessing
    infiles = []
    if os.path.isdir(credset_dir):
        for file in os.listdir(credset_dir):
            if file.endswith('.cred'):
                splitfile = file.split('.')
                infiles.append(((os.path.join(credset_dir, file), splitfile[1]), build, GenomicRiskLoci, Ann_path, outdir, lo, magma, PoPS, genes, magma_scores, TP))
    else:
        with open(credset_dir, 'r') as f:
            for l in f:
                splitline = l.split('\t')
                infiles.append(((splitline[0], splitline[1].strip()),build, GenomicRiskLoci, Ann_path, outdir, lo, magma, PoPS, genes, magma_scores, TP))
                print(infiles)
    if processors:
        pool = mp.Pool(processes=processors)
    else:
        pool = mp.Pool()

    # apply the function to each file in parallel, passing additional arguments
    pool.starmap(full_annotation_of_credset, infiles)

    # wait for the worker processes to complete
    pool.close()
    pool.join()
    return

if __name__ == "__main__":
    main(argv[1:])