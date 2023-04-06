#import modules
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
import os
import glob
from pyliftover import LiftOver
import numpy as np
from tqdm import tqdm
from sys import argv
import sys
import time

###define functions

"""Parses a finemap_log output file to extract the posterior probabilities for k-causal SNPS"""
def parse_finemap_logs(log_file_path):
    with open(log_file_path) as f:
        posteriors = False
        posteriors_list = []
        for line in f:
            if "- Post-Pr(# of causal SNPs is k)   :" in line:
                posteriors = True
            if posteriors and "->" in line:
                if not "(" in line:
                    posteriors_list.append(line.split('->')[1].strip())
        return(posteriors_list)


"""Given a the FINEMAP locus log file and locus cred file in the same location:
    parses the cred file to return the credible sets of possible causal variants.
    The returned cred set is a df with variant ids and posterior probs  """
def get_credible_sets(log_file_path):
    posteriors = parse_finemap_logs(log_file_path)
    credsets_no = posteriors.index(max(posteriors)) + 1
    credfile = log_file_path.split('.')[0]
    credfile = glob.glob(credfile+ ".cred" + str(credsets_no))[0]
    credsets = pd.read_csv(credfile, comment='#', sep=' ')
    return credsets, credsets_no 

"""For any given locus, retrieve the genes present in that locus  """
def get_genes_in_locus(locus_no, path, GenomicRiskLoci):
    genes = pd.read_csv(path, sep='\t')
    genes['GenomicLocus'] = genes.GenomicLocus.str.split(":")
    genes = genes.explode('GenomicLocus')
    relevant_genes = pd.DataFrame()
    locusname = str(GenomicRiskLoci.iloc[int(locus_no)-1][0])
    for locus in locusname.split(':'):
        genes_in_loc = genes[genes['GenomicLocus'] == locus]
        relevant_genes = pd.concat([relevant_genes, genes_in_loc])
    return relevant_genes

""" Liftover row values from hg19 to hg38 """
def liftover_rows(row):
    lo = LiftOver('hg19', 'hg38')
    hg38_coordinates = lo.convert_coordinate(f"chr{row['chr']}", int(row['pos']))[0]
    return hg38_coordinates[0], hg38_coordinates[1]

"""For the genes in a given finemapped locus, determine the distance
    between the TSS and each finemapped variant, sum the distances weighted by posterior probabilities  """
def get_TSS_distances(credset, genes_in_locus):
    distances = []
    leadvar_distance =[]
    credset = credset.dropna()
    lo = LiftOver('hg19', 'hg38')
    for index, gene in genes_in_locus.iterrows():
        distance = 0
        for index, row in credset.iterrows():
            hg38_coordinates = lo.convert_coordinate(f"chr{row['chr']}", int(row['pos']))
            if distance == 0:
                if int(gene['strand']) == -1:
                    leadvar_distance.append(abs(int(gene['end']) - int(hg38_coordinates[0][1])))
                elif int(gene['strand']) == 1:
                    leadvar_distance.append(abs(int(gene['start']) - int(hg38_coordinates[0][1])))
            if int(gene['strand']) == -1:
                distance += abs(int(gene['end']) - int(hg38_coordinates[0][1])) * float(row[2])
            elif int(gene['strand']) == 1:
                distance += abs(int(gene['start']) - int(hg38_coordinates[0][1])) * float(row[2])
        distances.append(distance)
    genes_in_locus['weighted_distance_TSS_credset'] = distances
    genes_in_locus['distance_TSS_leadvar'] = leadvar_distance
    return(genes_in_locus)


def get_GTEX_eQTLs(cred_w_prob, genes, finemapped_gtex):
    chr_38=[]
    pos_38=[]
    lo = LiftOver('hg19', 'hg38')
    for index, row in cred_w_prob.iterrows():
        hg38_coordinates = lo.convert_coordinate(f"chr{row['chr']}", int(row['pos']))[0]
        chr_38.append(int(hg38_coordinates[0].strip('chr')))
        pos_38.append(int(hg38_coordinates[1]))
    cred_w_prob["chr_38"] = chr_38
    cred_w_prob["pos_38"] = pos_38
    cred_w_prob['pos'] =cred_w_prob['pos'].astype(int)

    #GTEX
    overlap = pd.merge(cred_w_prob, finemapped_gtex, left_on=['chr_38','pos_38'], right_on=['CHROM', 'POS'])
    CLPPs = []
    ENSGs = []
    print('Get eQTL CLPP for each gene with overlapping finemapped variants GWAS & GTEx ')
    with tqdm(total=overlap.shape[0]) as pbar:    
        for index, row in overlap.iterrows():
            pbar.update(1)
            CLPPs.append(row['Probability'] * row['prob'])
            ENSGs.append(row['GENE'].split('.')[0])
    CLPP = pd.DataFrame({'ensg':ENSGs, 'CLPP_GTEX':CLPPs})
    CLPP = CLPP.groupby('ensg').agg({'CLPP_GTEX': ['sum', 'max']}).droplevel(axis=1, level=0).reset_index().rename(columns = {'sum':'CLPP_GTEX_eQTL_sum', 'max':'CLPP_GTEX_eQTL_max'})
    genes = pd.merge(genes,CLPP, on='ensg', how='left')
    genes[['CLPP_GTEX_eQTL_sum','CLPP_GTEX_eQTL_max']] = genes[['CLPP_GTEX_eQTL_sum','CLPP_GTEX_eQTL_max']].fillna(0)
    return genes

def get_eqtlgen_eQTLs(cred_w_prob, genes, eqtlgen_dir):
    #eQTLgen
    #this is tested, this works
    eQTL_dir = eqtlgen_dir
    chrom = round(sum(genes["chr"])/len(genes["chr"]))
    CLPP_max = []
    CLPP_sum = []
    #for each gene in locus:
    print('Get eQTLGen CLPP for each gene')
    with tqdm(total=genes.shape[0]) as pbar:    
        for index, row in genes.iterrows():
            pbar.update(1)
            ensg = row['ensg'].upper()
            fname = os.path.join(f"{eQTL_dir}", f"chr{chrom}", "post" ,f"{ensg}_new.out_post")
            # if that gene has finemapped eQTLgen variants
            if os.path.isfile(fname):
                posteriors = pd.read_csv(fname, sep="\t", engine='pyarrow')
                posteriors['SNP_ID'] = posteriors['SNP_ID'].astype(np.int32)
                #merge finemapped eQTLgen variants with locus finemapping variants
                gene_CLPP = pd.merge(cred_w_prob, posteriors, left_on=['pos'], right_on=['SNP_ID'])
                #if there are no finemapped eqtlgen variants which are also finemapped variants in locus then CLPP = 0
                if gene_CLPP.shape[0] == 0:
                    CLPP_max.append(0.0)
                    CLPP_sum.append(0.0)
                #if there are finemapped eqtlgen variants which are also finemapped variants in locus then CLPP = sum posterior probabilities
                else:
                    CLPPs = []
                    for index, nrow in gene_CLPP.iterrows():
                        CLPPs.append(nrow['Causal_Post._Prob.'] * nrow['prob'])
                    CLPP_max.append(max(CLPPs))
                    CLPP_sum.append(sum(CLPPs))
            else:
                CLPP_max.append(0.0)
                CLPP_sum.append(0.0)
    genes['CLPP_eQTLgen_sum'] = CLPP_sum
    genes['CLPP_eQTLgen_max'] = CLPP_max
    return genes
    
def get_CADD_scores(snps_txt, genes, cred_w_prob):
    CADD = snps_txt
    CADD['chr'] = CADD['chr'].astype(str)
    CADD = CADD.merge(cred_w_prob[['chr', 'pos']], on=['chr', 'pos'])
    CADDS = pd.DataFrame()
    print('Get CADD for each gene:')
    with tqdm(total=genes.shape[0]) as pbar:    
        for index, row in genes.iterrows():
            pbar.update(1)
            variants_in_gene = CADD[CADD['pos'] >= row['start']-1000]
            variants_in_gene = variants_in_gene[variants_in_gene['pos'] <= row['end']+1000]
            variants_in_gene['ensg'] = row['ensg']
            CADDS = pd.concat([CADDS, variants_in_gene])
    CADDS = CADDS.groupby('ensg').agg({'CADD': ['sum', 'max']}).droplevel(axis=1, level=0).reset_index().rename(columns = {'sum':'CADD_sum', 'max':'CADD_max'})
    genes = pd.merge(genes,CADDS, on='ensg', how='left')
    genes[['CADD_sum','CADD_max']] = genes[['CADD_sum','CADD_max']].fillna(0)
    return genes

def prep_HiC_data(jung, javierre):
    #jung     
    jung = pd.read_csv(jung, sep='\t', skiprows=1, names=['cell', 'chr','gene', 'target', 'raw', 'dist', 'all_capture_res', 'dist_res', 'dist_pvalue'], engine='pyarrow')
    jung['gene'] = jung['gene'].apply(lambda x: x.strip(','))
    jung['gene'] = jung['gene'].apply(lambda x: x.split(','))
    jung = jung.explode('gene')
    #javierre
    javierre = pd.read_csv(javierre, sep='\t', engine='pyarrow')
    javierre.dropna(subset=['baitName'])
    javierre['baitName'] = javierre['baitName'].astype(str)
    javierre['baitName'] = javierre['baitName'].apply(lambda x: x.split(';'))
    javierre = javierre.explode('baitName')
    return jung,javierre

def get_jung_HiC(jung, genes, cred_w_prob):
    jung_max = []
    jung_sum = []
    print('Get jung HiC interactions for each gene')
    with tqdm(total=genes.shape[0]) as pbar:    
        for index, row in genes.iterrows():
            pbar.update(1)
            overlap_jung = jung[jung['gene'] == row['symbol']]
            if overlap_jung.shape[0] == 0:
                jung_max.append(0.0)
                jung_sum.append(0.0)
            else:
                totals = pd.DataFrame()
                overlap_jung[['chr','start', 'end']] = overlap_jung.target.str.split(".",expand=True,)
                overlap_jung['start'] = overlap_jung['start'].astype(int)
                overlap_jung['end'] = overlap_jung['end'].astype(int)
                for index, row in cred_w_prob.iterrows():
                    df = overlap_jung[overlap_jung['start'] <= row['pos']]
                    df = df[df['end'] >= row['pos']] 
                    totals = pd.concat([totals, df])
                totals= totals.drop_duplicates()
                if totals.shape[0] == 0:
                    jung_max.append(0.0)
                    jung_sum.append(0.0)
                else:
                    jung_max.append(max(totals['dist_pvalue']))
                    jung_sum.append(sum(totals['dist_pvalue']))
    genes['Jung_HiC_sum'] = jung_sum
    genes['Jung_HiC_max'] = jung_max
    return genes

def get_javierre_HiC(javierre, genes, cred_w_prob):
    javierre_max = []
    javierre_sum = []
    print('Get javierre HiC interactions for each gene')
    with tqdm(total=genes.shape[0]) as pbar:    
        for index, row in genes.iterrows():
            pbar.update(1)
            Chicagos = []
            overlap = javierre[javierre['baitName'] == row['symbol']]
            if overlap.shape[0] == 0:
                javierre_max.append(0.0)
                javierre_sum.append(0.0)
            else:
                totals = pd.DataFrame()
                overlap['oeStart'] = overlap['oeStart'].astype(int)
                overlap['oeEnd'] = overlap['oeEnd'].astype(int)
                for index, row in cred_w_prob.iterrows():
                    df = overlap[overlap['oeStart'] <= int(row['pos'])]
                    df['oeChr'] = df['oeChr'].astype(str)
                    df = df[df['oeEnd'] >= int(row['pos'])] 
                    df = df[df['oeChr'] == str(row['chr'])]
                    totals = pd.concat([totals, df])
                totals= totals.drop_duplicates()
                if totals.shape[0] == 0:
                    javierre_max.append(0.0)
                    javierre_sum.append(0.0)
                else:
                    celltypes = [ 'Mon', 'Mac0', 'Mac1','Mac2', 'Neu', 'MK', 'EP', 'Ery', 'FoeT', 'nCD4', 'tCD4', 'aCD4','naCD4', 'nCD8', 'tCD8', 'nB', 'tB']
                    for ct in celltypes:
                        Chicagos.extend(list(totals[ct]))
                    javierre_max.append(max(Chicagos))
                    javierre_sum.append(sum(Chicagos))
    genes['Javierre_HiC_sum'] = javierre_sum
    genes['Javierre_HiC_max'] = javierre_max
    return(genes)

def prep_Fanthom5_data(F5_organ, F5_celltype,chrom):
    F5_ct = pd.read_csv(F5_celltype, sep='\t', engine='pyarrow')
    F5_organ = pd.read_csv(F5_organ, sep='\t', engine='pyarrow')

    #load and transform celltype data
    F5_ct[['p_chr', 'p_coords']] = F5_ct.promoter.str.split(":",expand=True)
    F5_ct = F5_ct[F5_ct['p_chr'] == 'chr' + str(chrom)]
    F5_ct[['p_start','p_stop']] = F5_ct.p_coords.str.split("\.\.", regex=True,expand=True,)
    F5_ct['p_stop'] = F5_ct.p_stop.str.split(",").str[0]
    F5_ct[['e_chr','e_coords']] = F5_ct.enhancer.str.split(":",expand=True,)
    F5_ct = F5_ct[F5_ct['e_chr'] == 'chr' + str(chrom)]
    F5_ct[['e_start','e_stop']] = F5_ct.e_coords.str.split("-",expand=True,)

    #load and transform organ data
    F5_organ[['p_chr', 'p_coords']] = F5_organ.promoter.str.split(":",expand=True)
    F5_organ = F5_organ[F5_organ['p_chr'] == 'chr' + str(chrom)]
    F5_organ[['p_start','p_stop']] = F5_organ.p_coords.str.split("\.\.", regex=True,expand=True,)
    F5_organ['p_stop'] = F5_organ.p_stop.str.split(",").str[0]
    F5_organ[['e_chr','e_coords']] = F5_organ.enhancer.str.split(":",expand=True,)
    F5_organ = F5_organ[F5_organ['e_chr'] == 'chr' + str(chrom)]
    F5_organ[['e_start','e_stop']] = F5_organ.e_coords.str.split("-",expand=True,)

    #merge and tidy
    F5_ep = pd.concat([F5_organ, F5_ct])
    F5_ep = F5_ep[['correlation', 'e_chr', 'e_start', 'e_stop', 'p_chr', 'p_start', 'p_stop']]
    return F5_ep

def get_F5_ep(F5_ep, genes, cred_w_prob):
    print('Get Fanthom5 promoter enhancer correlation for each finemapped variant in enhancer to promoter of gene in locus')
    with tqdm(total=genes.shape[0]) as pbar:  
        EP_cor_sum = []
        EP_cor_max = []
        for index, row in genes.iterrows():
            pbar.update(1)
            if int(row['strand']) == 1:
                EP_overlap = F5_ep[F5_ep['p_start'].astype(int) >= int(row['start']) - 500 ]
                EP_overlap = EP_overlap[EP_overlap['p_start'].astype(int) < int(row['start']) ]
            else:
                EP_overlap = F5_ep[F5_ep['p_stop'].astype(int) <= int(row['end']) + 500]
                EP_overlap = EP_overlap[EP_overlap['p_stop'].astype(int) > int(row['end'])]
            if EP_overlap.shape[0] == 0:
                EP_cor_sum.append(0.0)
                EP_cor_max.append(0.0)
            else:
                totals = []
                for index, row in cred_w_prob.iterrows():
                    df= EP_overlap
                    #df = EP_overlap[EP_overlap['p_chr'] == 'chr'+ str(row['chr'])]
                    df = df[df['e_start'].astype(int) <= int(row['pos'])]
                    df = df[df['e_stop'].astype(int) >= int(row['pos'])] 
                    if df.shape[0] > 0:
                        totals.extend(list(df['correlation']))
                if len(totals) > 0:
                    EP_cor_max.append(max(totals))
                    EP_cor_sum.append(sum(totals))
                else:
                    EP_cor_sum.append(0.0)
                    EP_cor_max.append(0.0)
                    
    genes['Fanthom_EP_cor_sum'] = EP_cor_sum
    genes['Fanthom_EP_cor_max'] = EP_cor_max
    return genes

def get_activity_by_contact(abc, genes, cred_w_prob, chrom):
    ABC_max = []
    ABC_sum = []
    if chrom == 23:
        chrom = 'X'
    abc = abc[abc['chr'].astype(str) == 'chr' + str(chrom)]
    print('For gene in locus, get activity-by-contact scores to finemapped variants')
    with tqdm(total=genes.shape[0]) as pbar:    
        for index, row in genes.iterrows():
            pbar.update(1)
            overlap = abc[abc['TargetGene'] == row['symbol']]
            if overlap.shape[0] == 0:
                ABC_max.append(0.0)
                ABC_sum.append(0.0)
            else:
                totals = pd.DataFrame()
                for index, row in cred_w_prob.iterrows():
                    df = overlap[overlap['start'] <= row['pos']]
                    df = df[df['end'] >= row['pos']] 
                    totals = pd.concat([totals, df])
                totals= totals.drop_duplicates()
                if totals.shape[0] == 0:
                    ABC_max.append(0.0)
                    ABC_sum.append(0.0)
                else:
                    ABC_max.append(max(totals['ABC.Score']))
                    ABC_sum.append(sum(totals['ABC.Score']))
    genes['ABC_sum'] = ABC_sum
    genes['ABC_max'] = ABC_max
    return genes

def get_epi(Epi, genes, cred_w_prob, chrom):
    Epi_max = []
    Epi_sum = []
    print('For gene in locus, get epimap interactions to finemapped variants')
    Epi_chrom = Epi[Epi['chr'] == 'chr'+str(chrom)]
    min_end = min(list(cred_w_prob['pos'])) -1
    max_start = max(list(cred_w_prob['pos'])) +1
    Epi_chrom = Epi_chrom[Epi_chrom['end'] > min_end]
    Epi_chrom = Epi_chrom[Epi_chrom['start'] < max_start]
    with tqdm(total=genes.shape[0]) as pbar:    
        for index, row in genes.iterrows():
            pbar.update(1)
            overlap = Epi_chrom[Epi_chrom['gene'] == row['ensg']]
            if overlap.shape[0] == 0:
                Epi_sum.append(0.0)
                Epi_max.append(0.0)
            else:
                totals = pd.DataFrame()
                for index, row in cred_w_prob.iterrows():
                    df = overlap[overlap['start'] <= row['pos']]
                    df = df[df['end'] >= row['pos']] 
                    totals = pd.concat([totals, df])
                totals= totals.drop_duplicates()
                if totals.shape[0] == 0:
                    Epi_max.append(0.0)
                    Epi_sum.append(0.0)
                else:
                    Epi_max.append(max(totals['score']))
                    Epi_sum.append(sum(totals['score']))
                    
    genes['Epi_max'] = Epi_max
    genes['Epi_sum'] = Epi_sum
    return genes

def genes_to_ML_features(GenomicRiskLoci, loci_gene_index, genes, locus_no, basedir, cred_no, extension, filtname):
    loci_gene_index = loci_gene_index.astype(str)
    TP_gene_in_locus = loci_gene_index.to_dict()['symbol']
    locusname = str(GenomicRiskLoci.iloc[int(locus_no)-1][0])
    try:
        true_pos=TP_gene_in_locus[locusname]
    except:
        print(f'{locusname} not found in loci_gene_index in {basedir}')
        print(TP_gene_in_locus)
        exit
    true_pos=TP_gene_in_locus[locusname]
    genes['TruePositive'] = np.where(genes['symbol'] in true_pos, 1, 0)
    genes['rel_distance_TSS_leadvar'] = min(genes['distance_TSS_leadvar'])/genes['distance_TSS_leadvar']
    genes['rel_weighted_distance_TSS_credset'] = min(genes['weighted_distance_TSS_credset'])/genes['weighted_distance_TSS_credset']
    for feature in ['CLPP_GTEX_eQTL_sum', 'CLPP_GTEX_eQTL_max',
        'CLPP_eQTLgen_sum', 'CLPP_eQTLgen_max', 'CADD_sum', 'CADD_max',
        'Jung_HiC_sum', 'Jung_HiC_max', 'Javierre_HiC_sum', 'Javierre_HiC_max',
        'Fanthom_EP_cor_sum', 'Fanthom_EP_cor_max', 'ABC_sum', 'ABC_max',
        'Epi_max', 'Epi_sum']:
        if max(genes[feature].astype(float)) == 0.0:
            genes[f'rel_{feature}'] = genes[feature]
        else:
            genes[f'rel_{feature}'] = genes[feature].astype(float)/max(genes[feature].astype(float))
    print(os.path.join(basedir, 'finemap_output', f'{extension}.locus_{locus_no}.{cred_no}_ML_features{filtname}.tsv'))
    genes.to_csv(os.path.join(basedir, 'finemap_output', f'{extension}.locus_{locus_no}.{cred_no}_ML_features{filtname}.tsv'), sep='\t')
    return

def timer(start_time, process_name):
    timestamp = time.perf_counter() - start_time
    start_time = time.perf_counter()
    m, s = divmod(timestamp, 60)
    h, m = divmod(m, 60)
    timestamp_f = f'{int(h):d}:{int(m):02d}:{int(s):02d}'
    print(f'time elapsed for {process_name}: {timestamp_f}')
    return start_time

def main(input_log_path):
    #inputs:
    #annotations basedir
    #path to log file

    start_time = time.perf_counter()

    input_log_path= input_log_path

    basedir, locus_no = input_log_path.split('locus_')
    basedir = basedir.split('finemap/')[0]
    extension = os.path.normpath(basedir.split('finemapped_fuma/')[1])
    locus_no = locus_no.split('.')[0]
    GenomicRiskLoci = pd.read_csv(os.path.join(basedir,'finemap','GenomicRiskLociFinemap.txt'), sep='\t', engine='pyarrow')
    if os.path.exists(os.path.join(basedir,'loci_gene_index')):
        loci_gene_index = pd.read_csv(os.path.join(basedir,'loci_gene_index'),index_col=[0], header=None, sep = ' ', engine='pyarrow')
    else: 
        genes = pd.read_csv(os.path.join(basedir, 'genes.txt'), sep='\t', engine='pyarrow')
        fgenes = []
        with open(os.path.join(basedir,'geneslist.txt')) as f:
            for l in f:
                fgenes.append(l.strip())
        loci_gene_index = genes[(genes['symbol'].isin(fgenes)) | (genes['ensg'].isin(fgenes))].reset_index()
        print(loci_gene_index)
        loci_gene_index['GenomicLocus'] = loci_gene_index.GenomicLocus.str.split(":")
        print(loci_gene_index)
        loci_gene_index = loci_gene_index[['symbol', 'GenomicLocus']]
        loci_gene_index = loci_gene_index.explode('GenomicLocus')
        loci_gene_index = loci_gene_index.set_index('GenomicLocus')
        print(loci_gene_index)
    #Read in credible sets and genes in locus
    cred_sets, credsets_no = get_credible_sets(input_log_path)
    base_genes = get_genes_in_locus(locus_no, os.path.join(basedir, 'genes.txt'),GenomicRiskLoci)
    if base_genes.shape[0] > 0:
        chrom = round(sum(base_genes["chr"])/len(base_genes["chr"]))
    else:
        print(f'no genes in {input_log_path}')
        return

    #define paths to annotation data, load in pandas to save repeated loading in where possible
    annot_dir = '/home/schipper/ML/Annotation_data'
    
    #timer
    start_time = timer(start_time, 'reading locus data')

    finemapped_gtex = pd.read_csv(os.path.join(annot_dir, 'GTEx_v8_finemapping_CAVIAR', 'CAVIAR_Results_v8_GTEx_LD_ALL_NOCUTOFF_with_Allele.txt'), sep = "\t", engine="pyarrow")
    #timer
    start_time = timer(start_time, 'reading GTEx')

    eqtlgen_dir = os.path.join(annot_dir, 'eQTLGen', "caviar_farhad")
    snps_txt = pd.read_csv(os.path.join(basedir, "snps.txt"), sep='\t', engine='pyarrow')
    path_to_jung_PCHiC = os.path.join(annot_dir,'Jung_PCHiC', 'GSE86189_all_interaction.po.txt')
    path_to_javierre_PCHiC = os.path.join(annot_dir, 'Javierre_PCHiC','PCHiC_peak_matrix_cutoff5.tsv')
    jung,javierre = prep_HiC_data(path_to_jung_PCHiC, path_to_javierre_PCHiC)
    F5_organ = os.path.join(annot_dir, 'EP_correlation_F5', 'hg19_enhancer_promoter_correlations_distances_organ.txt')
    F5_celltype = os.path.join(annot_dir, 'EP_correlation_F5', 'hg19_enhancer_promoter_correlations_distances_cell_type.txt')
    Fanthom5_EP_cor = prep_Fanthom5_data(F5_organ, F5_celltype, chrom)
    activity_by_contact = pd.read_csv(os.path.join(annot_dir, 'ABC_EP', 'ABC_EP_reduced_columns.txt'), sep='\t', engine='pyarrow')

    #timer
    start_time = timer(start_time, 'reading EpiMap')
    Epi = pd.read_csv(os.path.join(annot_dir,'EpiMap','GeneHancer_all_tissues.txt'), sep='\t', engine='pyarrow')

    #timer
    start_time = timer(start_time, 'reading annotation data')

    #annotate genes per credible set, normally, now credset 1, as we only use these for true positives.
    #for i in range(credsets_no):
    for i in [0]:
        i = i + 1
        cred_w_prob = cred_sets[[f'cred{i}', f'prob{i}']]
        cred_w_prob = cred_w_prob.rename(columns={f'cred{i}': "cred", f'prob{i}': "prob"})
        #if cred_w_prob.iloc[0]['prob'] > 0.05:
        #    cred_w_prob = cred_w_prob[cred_w_prob['prob']> 0.05]
        if cred_w_prob.shape[0] > 0:
            filtname = 'nofilter'
            cred_w_prob[['chr','pos', 'allelle']] = cred_w_prob.cred.str.split(":",expand=True,)
            cred_w_prob = cred_w_prob.dropna()
            genes = get_TSS_distances(cred_w_prob[['chr', 'pos', 'prob']], base_genes)
            genes = get_GTEX_eQTLs(cred_w_prob, genes, finemapped_gtex)
            genes = get_eqtlgen_eQTLs(cred_w_prob, genes, eqtlgen_dir)
            genes = get_CADD_scores(snps_txt, genes, cred_w_prob)
            genes = get_jung_HiC(jung, genes, cred_w_prob)
            genes = get_javierre_HiC(javierre, genes, cred_w_prob)
            genes = get_F5_ep(Fanthom5_EP_cor, genes, cred_w_prob)
            genes = get_activity_by_contact(activity_by_contact, genes, cred_w_prob, chrom)
            genes = get_epi(Epi, genes, cred_w_prob, chrom)
            genes_to_ML_features(GenomicRiskLoci, loci_gene_index, genes, locus_no, basedir, i, extension, filtname)

    for i in [0]:
        i = i + 1
        cred_w_prob = cred_sets[[f'cred{i}', f'prob{i}']]
        cred_w_prob = cred_w_prob.rename(columns={f'cred{i}': "cred", f'prob{i}': "prob"})
        cred_w_prob = cred_w_prob[cred_w_prob['prob']> 0.05]
        if cred_w_prob.shape[0] > 0:
            filtname = 'pip05_filter'
            cred_w_prob[['chr','pos', 'allelle']] = cred_w_prob.cred.str.split(":",expand=True,)
            cred_w_prob = cred_w_prob.dropna()
            genes = get_TSS_distances(cred_w_prob[['chr', 'pos', 'prob']], base_genes)
            genes = get_GTEX_eQTLs(cred_w_prob, genes, finemapped_gtex)
            genes = get_eqtlgen_eQTLs(cred_w_prob, genes, eqtlgen_dir)
            genes = get_CADD_scores(snps_txt, genes, cred_w_prob)
            genes = get_jung_HiC(jung, genes, cred_w_prob)
            genes = get_javierre_HiC(javierre, genes, cred_w_prob)
            genes = get_F5_ep(Fanthom5_EP_cor, genes, cred_w_prob)
            genes = get_activity_by_contact(activity_by_contact, genes, cred_w_prob, chrom)
            genes = get_epi(Epi, genes, cred_w_prob, chrom)
            genes_to_ML_features(GenomicRiskLoci, loci_gene_index, genes, locus_no, basedir, i, extension,filtname)

    for i in [0]:
        i = i + 1
        cred_w_prob = cred_sets[[f'cred{i}', f'prob{i}']]
        cred_w_prob = cred_w_prob.rename(columns={f'cred{i}': "cred", f'prob{i}': "prob"})
        if cred_w_prob.shape[0] > 25:
            filtname = 'top25_filter'
            cred_w_prob = cred_w_prob.head(25)
        cred_w_prob[['chr','pos', 'allelle']] = cred_w_prob.cred.str.split(":",expand=True,)
        cred_w_prob = cred_w_prob.dropna()
        genes = get_TSS_distances(cred_w_prob[['chr', 'pos', 'prob']], base_genes)
        genes = get_GTEX_eQTLs(cred_w_prob, genes, finemapped_gtex)
        genes = get_eqtlgen_eQTLs(cred_w_prob, genes, eqtlgen_dir)
        genes = get_CADD_scores(snps_txt, genes, cred_w_prob)
        genes = get_jung_HiC(jung, genes, cred_w_prob)
        genes = get_javierre_HiC(javierre, genes, cred_w_prob)
        genes = get_F5_ep(Fanthom5_EP_cor, genes, cred_w_prob)
        genes = get_activity_by_contact(activity_by_contact, genes, cred_w_prob, chrom)
        genes = get_epi(Epi, genes, cred_w_prob, chrom)
        genes_to_ML_features(GenomicRiskLoci, loci_gene_index, genes, locus_no, basedir, i, extension, filtname)



if __name__ == "__main__":
    main(argv[1])
    