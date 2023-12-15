import pandas as pd
import numpy as np
from sys import argv
import pickle as pk
import os
import train_FLAMES as tf
import random
from sklearn.utils.parallel import delayed
import xgboost as xgb


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
    df["Combined_score"] = (1 - float(pops_combi)) * (df["rel_PoPS_Score"]) + float(
        pops_combi
    ) * df["rel_XGB_score"]
    max_in_loc = float(max(df["Combined_score"]))
    df["raw_FLAMES_score"] = df["Combined_score"]
    if len(df) == 1:
        df["FLAMES_score"] = df["XGB_score"]
    else:
        df["diff"] = df["Combined_score"] - df["Combined_score"].nlargest(2).values[-1]
        df["diff2"] = df["Combined_score"].nlargest(1).values[0] - df["Combined_score"]
        df["FLAMES_score"] = df["diff"] - df["diff2"]
    max_in_loc = float(max(df["FLAMES_score"]))
    df["FLAMES_causal"] = df["FLAMES_score"].apply(
        lambda x: 1 if x == max_in_loc and max_in_loc > 0 else 0
    )
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
        f1_corrected = f1_best / precision_closest_gene
    return best_params, f1_corrected


def main(model, trainfiles, outdir, filter):
    # Load training data and trained model for hyperparams of model
    XGB_model = pk.load(model)
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
    names = data_list
    data_list = [pd.read_csv(f, sep="\t") for f in data_list]
    for i in range(len(data_list)):
        data_list[i]["filename"] = names[i]
    split_data = np.array_split(data_list, num_splits)

    # Create a for loop to retrain each separate model on split data
    best_params = []
    for i, data_split in enumerate(split_data):
        train = []
        for j, x in enumerate(split_data):
            if j == i:
                continue
            train.extend(split_data[j])
        val = data_split
        print("read data")
        train = tf.add_into_single_df(train)
        train_X, train_y, df = tf.create_XY_train(train, filter)
        # Load features used in model
        if os.path.exists(outdir + "/features.txt"):
            with open(outdir + "/features.txt") as f:
                features = [l.strip() for l in f]
            train_X = train_X[features]
        # train model
        print("train_model")
        model = XGB_model.fit(train_X, train_y)
        val = tf.add_into_single_df(val)
        val_X, val_y, df = tf.create_XY_train(val, filter)
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
