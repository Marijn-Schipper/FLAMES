import pandas as pd
import numpy as np
from sys import argv
import pickle as pk
import os
import random
import xgboost as xgb
from sklearn.preprocessing import MinMaxScaler

def add_into_single_df(list_of_dfs):
    df = pd.concat(list_of_dfs)
    df = df.reset_index(drop=True)
    df["locus"] = pd.factorize(df["filename"])[0]
    return df

def create_XY_train(df, filter=None, inclusion="all"):
    train_features = [
        "genes_in_locus",
        "credset_size",
        "genes_within_50kb",
        "genes_within_100kb",
        "genes_within_250kb",
        "Highest_PIP",
        "rel_MAGMA_Z",
        "MAGMA_Z",
    ]
    best_predictors = [
        "weighted_distance",
        "TSS_distance",
        "weighted_TSS_distance",
        "VEP_max",
        "CADD_sum",
        "CLPP_GTEX_tissue_weighted_eQTL_sum",
        "RoadmapEpi_sum",
        "HACER_PRO-seq_GRO-seq_sum",
        "Promoter_sum",
        "F5_ep_sum",
        "Jung_HiC_sum",
        "GeneHancer_sum",
        "EpiMap_sum",
        "mean_CLPP_eQTL_eQTLCatalog_max",
        "ABC_CRISPR_sum",
        "ABC_EP_sum",
        "max_CLPP_rQTL_eQTLCatalog_max",
        "Javierre_HiC_sum",
        "HACER_CAGE_max",
        "CLPP_eQTLgen_sum",
        "cicero_sum",
    ]
    if filter != None:
        df = df[df["distance"] <= int(filter)]
    else:
        df = df[df["weighted_distance"] <= 750000]
        df = df.groupby("filename").apply(adjust_PoPS_scores)
    if "gene_biotype" in df.columns:
        df = df[df["gene_biotype"] == "protein_coding"]
    if "type" in df.columns:
        df = df[df["type"] == "protein_coding"]
    df = df.groupby("filename").apply(ann.create_relative_annotation)

def Flames_scoring(df, pops_combi):
    df = df.drop_duplicates(keep="first")
    df = df.drop_duplicates(subset=["ensg", "symbol"], keep="first")
    max_in_loc = float(max(df["XGB_score"]))
    df["rel_XGB_score"] = df["XGB_score"] / max_in_loc
    max_in_loc = float(max(df["XGB_score"]))
    df["PoPS_prediction"] = df["rel_PoPS_Score"].apply(
        lambda x: 1 if float(x) == float(1) else 0
    )
    df["Magma_pred"] = df["rel_MAGMA_Z"].apply(
        lambda x: 1 if float(x) == float(1) else 0
    )
    df["dist_pred"] = df["rel_distance"].apply(lambda x: 1 if x == 1 else 0)
    scaler = MinMaxScaler(feature_range=(0.292, 1))
    scaled_values = scaler.fit_transform(df[['PoPS_Score']])
    df['rel_PoPS_Score'] = scaled_values
    df['rel_PoPS_Score'] = df['rel_PoPS_Score'].round(10)
    df['rel_PoPS_Score'] = df['rel_PoPS_Score']/(max(df['rel_PoPS_Score']))
    if len(df) ==1 :
        df['FLAMES_scaled'] = 1
        df['FLAMES_raw'] = 1.583 * df['XGB_score']
    elif max(df['rel_PoPS_Score']) == 0:
        df['FLAMES_raw'] = df['XGB_score']
        df['FLAMES_scaled'] = ((df['FLAMES_raw']/sum(df['FLAMES_raw'])))
    else:
        df['FLAMES_raw'] =  df['XGB_score'] * df['rel_PoPS_Score'] +  0.583 * df['XGB_score'] * df['rel_PoPS_Score'].astype(int)
        df['FLAMES_scaled'] = ((df['FLAMES_raw']/sum(df['FLAMES_raw']))) 
    max_in_loc = max(df['FLAMES_scaled'])
    df["FLAMES_highest"] = df["FLAMES_scaled"]/max_in_loc
    df['FLAMES_highest'] = df['FLAMES_highest'].astype(int)
    df['FLAMES_causal'] = ((df['FLAMES_raw'] > 0.134) & (df['FLAMES_scaled'] > 0.248)).astype(int) * df['FLAMES_highest']
    #Based on calibration as performed in paper (simple polynomial fit)
    df['estimated_cumulative_precision'] = 1.5776340061236926 * df['FLAMES_scaled']**3 -2.761585223042816 * df['FLAMES_scaled']**2 + 1.8136337684824269 * df['FLAMES_scaled'] + 0.44600372410005074 
    df['estimated_cumulative_precision'] = df['estimated_cumulative_precision'].apply(lambda x: 1 if x > 1 else x)
    return df


def optimize_params(
    y_preds,
    df,
    iterations,
):
    best_params = None
    b2 = {}
    df["XGB_score"] = y_preds
    results = df.copy()
    f1_best = 0
    fscores = {}
    for i in range(iterations):
        PoPS_combination = np.linspace(0, 1, iterations)[i]
        df = results.groupby("locus").apply(
            Flames_scoring,
            PoPS_combination,
        )
        predictions = sum(df["FLAMES_causal"])
        if predictions == 0:
            continue
        causal = sum(df["TP"])
        correct = sum(df["TP"] & df["FLAMES_causal"])
        df["correct"] = df["TP"] * df["FLAMES_causal"]
        df["closest_gene_correct"] = df["TP"] & df["dist_pred"] == 1
        precision_closest_gene = sum(df["closest_gene_correct"]) / sum(df["dist_pred"])
        precision = correct / predictions
        recall = correct / causal
        if recall > 0:
            f1 = 2 * (precision * recall) / (precision + recall)
            if f1 > f1_best:
                f1_best = f1
                best_params = PoPS_combination
        else:
            f1 = 0
        # correct for simple closest gene prediction model as minimal baseline performance so that splits with only closest genes don't get overweighted
        if precision_closest_gene == 0:
            f1_corrected = 0
        else:
            f1_corrected = f1_best / precision_closest_gene
    return best_params, f1_corrected


def main(model, trainfiles, outdir, filter):
    # Load training data and trained model for hyperparams of model
    XGB_model = pk.load(open(model, "rb"))
    data_file = os.path.join(trainfiles)
    data_list = []
    # Randomize order of loci
    with open(data_file, "r") as file:
        for line in file:
            data_list.append(line.strip())
    rng = random.Random(42)
    rng.shuffle(data_list)
    # Split the list into 50 parts for cross-validation
    num_splits = 50
    split_data = np.array_split(data_list, num_splits)
    data_list = []
    for s in split_data:
        s_list = []
        for f in s:
            df_s = pd.read_csv(f, sep="\t")
            df_s["filename"] = f
            s_list.append(df_s)
        data_list.append(s_list)

    # Create a for loop to retrain each separate model on split data
    best_params = []
    for i, data_split in enumerate(data_list):
        train = []
        for j, x in enumerate(data_list):
            if j == i:
                continue
            train.extend(data_list[j])
        val = data_split
        train = add_into_single_df(train)
        train_X, train_y, df = create_XY_train(train, filter)
        # Load features used in model
        if os.path.exists(outdir + "/features.txt"):
            with open(outdir + "/features.txt") as f:
                features = [l.strip() for l in f]
            train_X = train_X[features]
        # train model
        model = XGB_model.fit(train_X, train_y)
        val = add_into_single_df(val)
        val_X, val_y, df = create_XY_train(val, filter)
        if os.path.exists(outdir + "/features.txt"):
            with open(outdir + "/features.txt") as f:
                features = [l.strip() for l in f]
            val_X = val_X[features]
            y_preds = model.predict_proba(val_X)[:, 1]
        print(f"optimize split {i}")
        # Find the best PoPS + XGB linear combination
        iterations = 41

        print(f"total amount of tps = {sum(train_y)}")
        best_param, f1_corrected = optimize_params(
            y_preds,
            df,
            iterations,
        )
        best_params.append(best_param)
    best_params = np.median(best_params)
    with open(outdir + "/best_params.txt", "w") as f:
        f.write(
            f"{f1_corrected} :XGBoost weight {str(best_params)}, PoPS weight {1-best_params}"
        )
    return


if __name__ == "__main__":
    main(argv[1:])
