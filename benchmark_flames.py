import pandas as pd
import numpy as np
from sys import argv
import pickle as pk
import os
import Train_FLAMES as tf
import matplotlib.pyplot as plt
import optimize_FLAMES as optf
from sklearn.metrics import auc, precision_recall_curve, roc_curve, roc_auc_score


def main(
    model_dir, input_dir, flames_thres, pops_combi, final_thresh, diff_thresh, filt
):
    print(input_dir)
    for file in os.listdir(model_dir):
        if ".sav" in file:
            model = file
    FLAMES = pk.load(open(model_dir + "/" + model, "rb"))
    infiles = []
    if os.path.isdir(input_dir):
        for file in os.listdir(input_dir):
            if "annotated_" in file:
                infiles.append(input_dir + "/" + file)
    elif input_dir == "val":
        input_dir = os.path.join(model_dir, "val.txt")
    if os.path.isfile(input_dir):
        print("aight")
        with open(input_dir) as f:
            for l in f:
                infiles.append(l.strip())
    dfs = [pd.read_csv(f, sep="\t") for f in infiles]
    for i in range(len(dfs)):
        dfs[i]["filename"] = infiles[i]
    df = tf.add_into_single_df(dfs)
    if not "TP" in list(df.columns):
        df["TP"] = 1
    if filt:
        X_val, y_val, df = tf.create_XY_train(df, filt)
    else:
        X_val, y_val, df = tf.create_XY_train(df)

    if "logreg" in model_dir or not os.path.exists(model_dir + "/features.txt"):
        X_val = tf.create_feature_for_logreg(X_val)
        features = X_val.columns
    elif os.path.exists(model_dir + "/features.txt"):
        with open(model_dir + "/features.txt") as f:
            features = [l.strip() for l in f]
        X_val = X_val[features]
    select_features = [
        feature
        for feature in features
        if "rel" in feature
        and not "PoPS" in feature
        and not "MAGMA" in feature
        and not "dist" in feature
    ]
    tmpdf = df[df["TP"] == 1]
    tmpdf = tmpdf[select_features]
    tmpdf["sum"] = tmpdf.sum(axis=1)
    for feature in X_val.columns:
        X_val[feature] = X_val[feature].astype(float)
    y_pred = FLAMES.predict_proba(X_val)[:, 1]
    print(y_pred)
    df["FLAMES_score"] = y_pred
    df = df.groupby("locus").apply(
        optf.Flames_scoring, flames_thres, pops_combi, final_thresh, float(diff_thresh)
    )
    if not os.path.exists(model_dir + "/feature_importances.png"):
        feature_names = X_val.columns
        importances = FLAMES.feature_importances_
        indices = np.argsort(importances)[::-1]
        plt.figure(figsize=(10, 6))

        # Plot the feature importances
        plt.bar(range(len(indices)), importances[indices], color="b", align="center")

        # Add feature names as x-axis labels
        plt.xticks(range(len(indices)), feature_names[indices], rotation="vertical")

        # Set labels and title
        plt.xlabel("Features")
        plt.ylabel("Importance")
        plt.title("Feature Importances")
        plt.savefig(model_dir + "/feature_importances.png", bbox_inches="tight")
    print(
        df[
            [
                "symbol",
                "FLAMES_score",
                "rel_PoPS_Score",
                "rel_FLAMES_score",
                "PoPS_combined",
                "diff",
                "MAGMA_Z",
                "PoPS_Score",
            ]
        ]
    )
    print(f'length of TP: {len(df[df["TP"] == 1])}')
    TP_genes = set(df[df["TP"] == 1]["ensg"])
    print(len(TP_genes))
    print(f"TP genes: {len(TP_genes)}")
    predictions = sum(df["FLAMES_causal"])
    df["Flames_causal_correct"] = df["FLAMES_causal"].astype(int) & df["TP"]
    print(set(df[df["FLAMES_causal"] == 1]["symbol"]))
    Flames_max_preds = list(df[df["Flames_causal_correct"] == 1]["symbol"])
    # print(Flames_max_preds)
    # print(list(df[df['FLAMES_causal'] == 1]['symbol']))
    precision = sum(df["Flames_causal_correct"]) / sum(df["FLAMES_causal"])
    recall = len(set(Flames_max_preds)) / len(TP_genes)

    # Theoretical_min_recall_all_predictions_and_correct = causal / predictions
    print(f"total predictions: {predictions}")
    print(f"total causal: {len(TP_genes)}")
    print(f"correct predictions: {sum(df['Flames_causal_correct'])}")
    print(f"precision: {precision}")
    print(
        f"total correct genes: {len(set(df[df['Flames_causal_correct'] == 1]['symbol']))}"
    )
    print(f"total causal genes: {len(TP_genes)}")
    print(f"recall: {recall}")

    df["closest_gene_correct"] = df["TP"] & df["dist_pred"] == 1
    closest_gene_preds = list(df[df["closest_gene_correct"] == 1]["symbol"])
    cg_precision = sum(df["closest_gene_correct"]) / sum(df["dist_pred"])

    cg_recall = len(set(closest_gene_preds)) / len(TP_genes)
    print(f"precision closest gene: {cg_precision}")
    print(f"recall closest gene: {cg_recall}")
    print(f"total closest gene predictions: {sum(df['dist_pred'])}")
    print(f"total closest gene correct: {sum(df['closest_gene_correct'])}")
    df["locusname"] = df["filename"].str.split(".").str[0]
    popsdf = df.drop_duplicates(subset=["locusname", "symbol"])
    # popsdf = df
    popsdf["PoPS_correct"] = (
        popsdf["rel_PoPS_Score"].fillna(0).astype(int) * popsdf["TP"]
    )
    PoPS_preds = list(popsdf[popsdf["PoPS_correct"] == 1]["symbol"])
    pops_precision = sum(popsdf["PoPS_correct"]) / len(
        popsdf[popsdf["rel_PoPS_Score"] == 1]
    )
    pops_recall = len(set(PoPS_preds)) / len(TP_genes)
    print(f"precision PoPS: {pops_precision}")
    print(f"recall PoPS: {pops_recall}")
    print(f"total PoPS predictions: {sum(popsdf['PoPS_prediction'])}")
    print(f"total PoPS correct: {sum(popsdf['PoPS_correct'])}")

    df["Magma_correct"] = df["Magma_pred"] & df["TP"]
    PoPS_preds = list(df[df["Magma_correct"] == 1]["symbol"])
    mgm_precision = sum(df["Magma_correct"]) / sum(df["Magma_pred"])
    mgm_recall = len(set(PoPS_preds)) / len(TP_genes)
    print(f"precision MAGMA: {mgm_precision}")
    print(f"recall MAGMA: {mgm_recall}")
    print(f"total MAGMA predictions: {sum(df['Magma_pred'])}")
    print(f"total MAGMA correct: {sum(df['Magma_correct'])}")
    df["Flames_max_correct"] = df["FLAMES_max"].astype(int) & df["TP"]
    Flames_max_preds = list(df[df["Flames_max_correct"] == 1]["symbol"])
    precision = sum(df["Flames_max_correct"]) / sum(df["FLAMES_max"])
    recall = len(set(Flames_max_preds)) / len(TP_genes)
    print(f"precision Flames_max: {precision}")
    print(f"recall Flames_max: {recall}")
    print(f"total Flames_max predictions: {sum(df['FLAMES_max'])}")
    print(f"total Flames_max correct: {sum(df['Flames_max_correct'])}")
    df["summed_max"] = df["summed_max"].fillna(0)
    df["sum_correct"] = df["TP"] & df["summed_max"].astype(int) == 1
    summed_correct = list(df[df["sum_correct"] == 1]["symbol"])
    precision = sum(df["sum_correct"]) / sum(df["summed_max"])
    recall = len(set(summed_correct)) / len(set(TP_genes))
    print(f"precision summed_score gene: {precision}")
    print(f"recall summed_score gene: {recall}")
    print(f"total summed_score predictions: {sum(df['summed_max'])}")
    print(f"total summed_score correct: {sum(df['sum_correct'])}")
    FLAMES_precisions = []
    FLAMES_recalls = []
    combi_score_full_precision = []
    combi_score_full_recall = []
    PoPS_precisions = []
    PoPS_recalls = []
    print("A")
    print(df[df["FLAMES_causal"].astype(int) == 1][["symbol"]])
    print(
        f'TP: {len(df[(df["FLAMES_causal"].astype(int) == 1 )& ( df["TP"] == 1)])/len(df[df["FLAMES_causal"].astype(int) == 1])}'
    )
    print("_-_-_-_-_-_-_")
    print(
        df[
            df["filename"]
            == "/home/schipper/FLAMES_SCZ/SCZ/finemap/annotated_credset_203.txt"
        ][["symbol", "rel_FLAMES_score", "diff", "rel_PoPS_Score"]]
    )

    ndf = df[(df["dist_pred"] < 1) & (df["FLAMES_causal"] == 1)]
    ndf["ender"] = ndf["filename"].str.split("/").str[-1]
    print(
        ndf[
            [
                "ender",
                "symbol",
                "diff",
                "TP",
                "FLAMES_causal",
                "FLAMES_score",
                "rel_FLAMES_score",
                "rel_PoPS_Score",
            ]
        ]
    )
    print(sum(ndf["TP"]) / len(ndf))
    diffs = []
    for i in np.linspace(5, -5, 3001):
        ndf = df[df["PoPS_combined"] >= i]
        try:
            precision1 = len(ndf[ndf["TP"] == 1]) / len(ndf)
        except ZeroDivisionError:
            precision1 = 0
        try:
            recall1 = len(ndf[ndf["TP"] == 1]) / len(df[df["TP"] == 1])
        except ZeroDivisionError:
            recall1 = 0
        i = round(i, 5)
        diffs.append(i)
        ndf = df[df["diff"] >= i]
        try:
            precision = len(ndf[ndf["TP"] == 1]) / len(ndf)
        except ZeroDivisionError:
            precision = 0
        try:
            recall = len(ndf[ndf["TP"] == 1]) / len(df[df["TP"] == 1])
        except ZeroDivisionError:
            recall = 0
        ndf = popsdf[popsdf["pdiff"] >= i]
        try:
            PoPS_precision = len(ndf[ndf["TP"] == 1]) / len(ndf)
        except ZeroDivisionError:
            PoPS_precision = 0
        try:
            PoPS_recall = len(ndf[ndf["TP"] == 1]) / len(df[df["TP"] == 1])
        except ZeroDivisionError:
            PoPS_recall = 0
        FLAMES_precisions.append(precision)
        FLAMES_recalls.append(recall)
        combi_score_full_precision.append(precision1)
        combi_score_full_recall.append(recall1)
        PoPS_precisions.append(PoPS_precision)
        PoPS_recalls.append(PoPS_recall)
    outdf = pd.DataFrame()
    outdf["FLAMES_precision"] = FLAMES_precisions
    outdf["FLAMES_recall"] = FLAMES_recalls
    outdf["combi_score_precision"] = combi_score_full_precision
    outdf["combi_score_recall"] = combi_score_full_recall
    outdf["PoPS_precision"] = PoPS_precisions
    outdf["PoPS_recall"] = PoPS_recalls
    outdf["diff"] = diffs
    outdf["closest_gene_precision"] = cg_precision
    outdf["closest_gene_recall"] = cg_recall
    outdf["pops_max_precision"] = pops_precision
    outdf["pops_max_recall"] = pops_recall
    outdf["magma_max_precision"] = mgm_precision
    outdf["magma_max_recall"] = mgm_recall
    name = argv[2]
    if "/" in argv[2]:
        name = argv[2].split("/")[-1]
    if "." in name:
        name = name.split(".")[0]
    if "_" in name:
        name = name.split("_")[-1]
    randoms = []
    for fname in df["filename"]:
        randoms.append(1 / len(df[df["filename"] == fname]))
    random_pred = np.mean(randoms)
    outdf["random"] = random_pred
    print("writing to file")
    outdf.to_csv(f"{model_dir}/FLAMES_precision_recall_{name}.csv")
    if not "PoPS_rank" in df.columns:
        df["PoPS_rank"] = "NA"
    traindf = df[
        [
            "ensg",
            "filename",
            "PoPS_rank",
            "rel_PoPS_Score",
            "FLAMES_score",
            "rel_FLAMES_score",
            "diff_first",
            "PoPS_Score",
            "TP",
            "diff",
            "rel_weighted_distance",
            "symbol",
            "rel_MAGMA_Z",
            "has_causal",
        ]
    ]
    traindf.to_csv(f"{model_dir}/retrain_{name}.csv")
    return


if __name__ == "__main__":
    try:
        filt = argv[7]
    except IndexError:
        filt = None
    main(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], filt)
