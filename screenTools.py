from ast import Str
from typing import Tuple
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import argparse
from datetime import datetime

# for a 96 well plate
ROWS_96 = ["A", "B", "C", "D", "E", "F", "G", "H"]
ROWS_384 = [
    "A",
    "A",
    "B",
    "B",
    "C",
    "C",
    "D",
    "D",
    "E",
    "E",
    "F",
    "F",
    "G",
    "G",
    "H",
    "H",
]
# list of numbers 1 to 12 repeating twice
COLUMNS_384 = []
for c in range(2):
    for i in range(1, 13):
        COLUMNS_384.append(str(i))


def parsePlate(well_data):
    """
    Codes column and row data of plate into a list format with each datapoint as a row.
    """
    well_data = well_data.reset_index(drop=True)
    values = []
    # wells = []
    rows = []
    columns = []
    for i, row in well_data.iterrows():
        for n, item in enumerate(row):
            values.append(item)
            # wells.append(ROWS_96[i]+str(n+1))
            rows.append(ROWS_96[i])
            columns.append(str(n + 1))
    # df = pd.DataFrame(wells, columns=['well'])
    df = pd.DataFrame(rows, columns=["row"])
    df["column"] = columns
    df["value"] = values

    return df


def expandPlate(plate_df, plate_number, substrate_number):
    """
    Reshape the df so that each cells is in a list.
    """
    l = []
    for i, row in plate_df.iterrows():
        for n, item in enumerate(row):
            l.append(item)
    df = pd.DataFrame(l, columns=["value"])
    df["plate_number"] = plate_number
    df["peptide"] = substrate_number

    return df


def generate384list():
    """
    Generates a list of 384 well locations in a dataframe. Columns are numbered and rows are letters.
    """
    rows = ROWS_384
    columns = COLUMNS_384
    columns = [str(i) for i in columns]
    wells = []
    for i, row in enumerate(rows):
        for n, item in enumerate(columns):
            wells.append(row + item)
    l = np.array(wells)

    # now reshape the array to have 24 columns and 16 rows
    l = l.reshape(16, 24)

    # return it as a dataframe
    df = pd.DataFrame(l, columns=columns, index=rows)

    return df


def importPlates(
    xls_path: Str,
    plate_list: list,
    peptide_list: list,
    positive_controls: list = [
        ["A1", "B1", "C1"],
        ["A1", "B1", "C1"],
    ],  # these are the positive controls for each plate
    negative_controls: list = [
        ["D1", "E1", "F1"],
        ["D1", "E1", "F1"],
    ],  # these are the negative controls for each plate
):
    """
    Imports a 384 well plate from an excel file. The first 12 columns are designated as plate 1 and
    the last 12 columns are designated as plate 2. Alternating rows are designated
    as substrate 1 and 2 respectively. The plate is then parsed into a list of 384 well plates.
    Parameters
    ----------
    xls_path : Str
        Path to the excel file containing the plate data.
    plate_list : list
        List of plate numbers. The first plate number corresponds to the first 96 wells of the plate.
    peptide_list : list
        List of peptides. The first peptide corresponds to the first 96 wells of the plate.
    positive_controls : list, optional
        List of positive control wells. The first list corresponds to the first plate and the second list corresponds to the second plate. The default is [["A1", "B1", "C1"], ["A1", "B1", "C1"]].
    negative_controls : list, optional
        List of negative control wells. The first list corresponds to the first plate and the second list corresponds to the second plate. The default is [["D1", "E1", "F1"], ["D1", "E1", "F1"]].
    Returns
    -------
    plates : list
        List of 384 well plates.
    """
    print("importing: ", xls_path)
    # upper_left_location = (50, 2) # these should be right if the plate reader saved correctly
    upper_left_location = (
        52,
        2,
    )  # these should be right if the plate reader saved correctly
    lower_right_location = (67, 25)

    df = pd.read_excel(xls_path)
    plate_df = df.iloc[
        upper_left_location[0] - 1 : lower_right_location[0],
        upper_left_location[1] : lower_right_location[1] + 1,
    ]  # this is the dataframe of the plate with the shape (16, 24)

    plates = parsePlate384(
        plate_df,
        plate_list,
        peptide_list,
        positive_controls=positive_controls,
        negative_controls=negative_controls,
    )

    return plates


def parsePlate384(
    plate_df: pd.DataFrame,
    plate_list: list,
    peptide_list: list,
    positive_controls: list = [["A1", "B1", "C1"], ["A1", "B1", "C1"]],
    negative_controls: list = [["D1", "E1", "F1"], ["D1", "E1", "F1"]],
):
    """
    Parses a 384 well plate. The first 12 columns are designated as plate 1 and
    the last 12 columns are designated as plate 2. Alternating rows are designated
    as substrate 1 and 2 respectively.
    """
    # get the first 12 columns of plate_df. These contain plate 1.
    plate_1 = plate_df.iloc[:, :12]
    # get the last 12 columns of plate_df. These contain plate 2.
    plate_2 = plate_df.iloc[:, 12:]
    # get the first set of alternating rows of plate_1. These contain substrate 1.
    plate_1_sub1 = plate_1.iloc[::2, :]
    # get the second set of alternating rows of plate_1. These contain substrate 2.
    plate_1_sub2 = plate_1.iloc[1::2, :]
    # get the first set of alternating rows of plate_2. These contain substrate 1.
    plate_2_sub1 = plate_2.iloc[::2, :]
    # get the second set of alternating rows of plate_2. These contain substrate 2.
    plate_2_sub2 = plate_2.iloc[1::2, :]

    # parse each plate as a 96 well plate
    plate_1_sub1 = parsePlate(plate_1_sub1)
    plate_1_sub2 = parsePlate(plate_1_sub2)
    plate_2_sub1 = parsePlate(plate_2_sub1)
    plate_2_sub2 = parsePlate(plate_2_sub2)

    # assign control wells to each plate
    plate_1_sub1 = assignControls(
        plate_1_sub1, positive=positive_controls[0], negative=negative_controls[0]
    )
    plate_1_sub2 = assignControls(
        plate_1_sub2, positive=positive_controls[0], negative=negative_controls[0]
    )
    plate_2_sub1 = assignControls(
        plate_2_sub1, positive=positive_controls[1], negative=negative_controls[1]
    )
    plate_2_sub2 = assignControls(
        plate_2_sub2, positive=positive_controls[1], negative=negative_controls[1]
    )

    # add the plate number and peptide to each plate
    plate_1_sub1["plate_number"] = plate_list[0]
    plate_1_sub1["peptide"] = peptide_list[0]
    plate_1_sub2["plate_number"] = plate_list[0]
    plate_1_sub2["peptide"] = peptide_list[1]
    plate_2_sub1["plate_number"] = plate_list[1]
    plate_2_sub1["peptide"] = peptide_list[0]
    plate_2_sub2["plate_number"] = plate_list[1]
    plate_2_sub2["peptide"] = peptide_list[1]

    # # now expand the plates into a list of wells
    # plate_1_sub1 = expandPlate(plate_1_sub1, plate_list[0], peptide_list[0])
    # plate_1_sub2 = expandPlate(plate_1_sub2, plate_list[0], peptide_list[1])
    # plate_2_sub1 = expandPlate(plate_2_sub1, plate_list[1], peptide_list[0])
    # plate_2_sub2 = expandPlate(plate_2_sub2, plate_list[1], peptide_list[1])

    # concatenate the plates into a single dataframe
    plate_df = pd.concat([plate_1_sub1, plate_1_sub2, plate_2_sub1, plate_2_sub2])

    return plate_df


def assignControls(plate_df, positive=None, negative=None):
    """
    Assigns positive and negative conditions to wells. Params are given as a list of well locations: ['A1', 'A2', 'A3']
    """
    if positive == None:
        positive = []
    if negative == None:
        negative = []

    conditions = []

    for i, row in plate_df.iterrows():
        if row["row"] + row["column"] in positive:
            conditions.append("positive")
        elif row["row"] + row["column"] in negative:
            conditions.append("negative")
        else:
            conditions.append("experimental")

    plate_df["condition"] = conditions

    return plate_df


def plotSpread(plate_df):
    n = list(range(plate_df.shape[0]))
    plate_df["n"] = n

    return sns.scatterplot(data=plate_df, x="n", y="value", hue="condition")


def pivotPlates(plate_df):
    """
    Pivot the data so that the matching wells are in the same row so that we can compare peptides to eachother.

    Does not operate in place.
    """
    compare = plate_df.pivot(
        index=["plate_number", "row", "column", "condition"],
        values=["value"],
        columns=["peptide"],
    )
    compare = compare.reset_index()

    return compare


def computeRatios(pivoted_data):
    """
    computes two ratios between respective peptides and adds two columns.
    """
    pep0 = pivoted_data["value"].iloc[:, 0]
    pep1 = pivoted_data["value"].iloc[:, 1]

    ratio0 = pep0 / pep1
    ratio1 = pep1 / pep0

    log0 = np.log10(ratio0)
    log1 = np.log10(ratio1)

    ratios_data = pivoted_data.copy()

    ratios_data[str(pep0.name) + "/" + str(pep1.name)] = log0
    ratios_data[str(pep1.name) + "/" + str(pep0.name)] = log1

    return ratios_data


def find_hits(
    data,
    stdev_threshold_brightness,
    stdev_threshold_selectivity_0,
    stdev_threshold_selectivity_1,
):
    """
    Finds hits and labels data based on the mean and stdev of the positive controls.
    """

    pep0 = data["value"].iloc[:, 0].name
    pep1 = data["value"].iloc[:, 1].name
    ratiostr0 = str(pep0) + "/" + str(pep1)
    ratiostr1 = str(pep1) + "/" + str(pep0)

    pos = data[data["condition"] == "positive"]

    mean0 = np.mean(pos["value"][pep0])
    stdev0 = np.std(pos["value"][pep0])
    mean_ratio0 = np.mean(pos[ratiostr0])
    stdev_ratio0 = np.std(pos[ratiostr0])

    mean1 = np.mean(pos["value"][pep1])
    stdev1 = np.std(pos["value"][pep1])
    mean_ratio1 = np.mean(pos[ratiostr1])
    stdev_ratio1 = np.std(pos[ratiostr1])

    hits = []

    for i, row in data.iterrows():
        # check if control
        if row["condition"].values != "experimental":
            hits.append(row["condition"].values[0])
        # check for pep0 success
        elif (row["value"][pep0] > stdev_threshold_brightness * stdev0 + mean0) and (
            row[ratiostr0].values
            > stdev_threshold_selectivity_0 * stdev_ratio0 + mean_ratio0
        ):
            hits.append(pep0)
        elif (row["value"][pep1] > stdev_threshold_brightness * stdev1 + mean1) and (
            row[ratiostr1].values
            > stdev_threshold_selectivity_1 * stdev_ratio1 + mean_ratio1
        ):
            hits.append(pep1)
        else:
            hits.append("not significant")

    return data.assign(to_pick=hits)


def find_hits_by_plate(
    data,
    stdev_threshold_brightness,
    stdev_threshold_selectivity_0,
    stdev_threshold_selectivity_1,
    pep0=None,
    pep1=None,
):
    """
    Finds hits using only the WT wells on the corresponding plate where the hit resides.
    """
    hitlist = []
    for platenum in set(data["plate_number"]):
        plate = data[data["plate_number"] == platenum]
        plate_hits = find_hits(
            plate,
            stdev_threshold_brightness,
            stdev_threshold_selectivity_0,
            stdev_threshold_selectivity_1,
            # pep0=pep0,
            # pep1=pep1,
        )
        hitlist.append(plate_hits)

    return pd.concat(hitlist)

def top_n_hits(plate_df, pep0n, pep1n, n=10, stdev_threshold=1):
    """
    Returns the top n hits for each plate in a dataframe.
    """
    # make the ratio strings
    ratiostr0 = str(pep0n) + "/" + str(pep1n)
    ratiostr1 = str(pep1n) + "/" + str(pep0n)

    # calculate the thresholds for each peptide by calculating the mean and adding a multiple of the standard deviation
    pep0_threshold = plate_df['value'][pep0n].mean() + stdev_threshold * plate_df['value'][pep0n].std()
    pep1_threshold = plate_df['value'][pep1n].mean() + stdev_threshold * plate_df['value'][pep1n].std()
    print("Thresholding plate", plate_df['plate_number'].unique()[0], "at", pep0_threshold, "for", pep0n,"and", pep1_threshold, "for", pep1n)

    # drop control wells
    plate_df = plate_df[plate_df['condition'] == 'experimental']

    # sort by selectivity for pep0
    pep0 = plate_df.sort_values(ratiostr0, ascending=False)
    # drop any values below the threshold and take the top n
    pep0 = pep0[pep0['value'][pep0n] > pep0_threshold].head(n)
    pep0['to_pick'] = pep0n
    
    # sort by selectivity for pep1
    pep1 = plate_df.sort_values(ratiostr1, ascending=False)
    # drop any values below the threshold and take the top n
    pep1 = pep1[pep1['value'][pep1n] > pep1_threshold].head(n)
    pep1['to_pick'] = pep1n

    # concatenate the two dataframes
    hits = pd.concat([pep0, pep1])

    # print plate number and number of hits if less than n*2
    if len(hits) < n*2:
        print("Only found", len(hits), "hits on plate", plate_df['plate_number'].unique()[0])

    return hits

def plot_performance(data, peptide0, peptide1):

    fig = plt.figure(figsize=(10, 10))
    peptide0 = str(peptide0)
    peptide1 = str(peptide1)
    ratiostr = peptide0 + "/" + peptide1

    sns.scatterplot(x=("value", peptide0), y=ratiostr, data=data, hue="to_pick")


def export_to_pick(
    data: pd.DataFrame, peptide0: Str, peptide1: Str, fname_prefix: Str = "./output/"
) -> pd.DataFrame:
    """
    To a folder named "output" in the current path:
    CSV containing the selected mutants.
    PNGs of each performance plot.
    """
    stamp = datetime.now().strftime("%Y.%m.%d-%H.%M")
    fname = fname_prefix + str(peptide0) + "_" + str(peptide1) + "_" + stamp

    plot_performance(data, peptide0, peptide1)
    # plt.savefig(fname + "_" + str(peptide0) + ".png")
    plt.savefig(fname + "_" + str(peptide0) + ".pdf")

    plot_performance(data, peptide1, peptide0)
    # plt.savefig(fname + "_" + str(peptide1) + ".png")
    plt.savefig(fname + "_" + str(peptide1) + ".pdf")

    trimmed_data = data[(data["to_pick"] == peptide0) | (data["to_pick"] == peptide1)]
    sorted_data = trimmed_data.sort_values(
        ["to_pick", "plate_number", "row", "column"]
    )[["to_pick", "plate_number", "row", "column"]]

    sorted_data.to_csv(fname + ".csv")

    return sorted_data

# prep and output hits as a csv output 160 at a time
def to_pick_csv(hits, file_prefix='to_pick'):
    # cast the plate_number column as an int
    hits['plate_number'] = hits['plate_number'].astype(int)
    # sort by plate_number
    hits = hits.sort_values('plate_number')
    # for each plate_number in the hits dataframe, add a new column with the deck position, starting at 1 and repeating after 8
    d = {hits['plate_number'].unique()[i]: i%8+1 for i in range(len(hits['plate_number'].unique()))}
    hits['deck_position'] = hits['plate_number'].map(d)

    # rename the plate_number column to plate_label
    hits = hits.rename(columns={'plate_number': 'plate_label'})
    # keep only the columns we want: to_pick, deck_position, row, column, and plate_label in that order
    hits = hits[['to_pick', 'deck_position', 'row', 'column', 'plate_label']]

    # output the hits as a csv file, splitting into a new file every eight unique plate_labels
    labels = hits['plate_label'].unique()
    for i in range(0, len(labels), 8):
        hits[hits['plate_label'].isin(labels[i:i+8])].to_csv(file_prefix + "_" + str(i) + ".csv", index=True)

    return hits

def main():

    parser = argparse.ArgumentParser(
        description="""
        Processes luminescence plates from a single excel file and outputs a list of "winners". It's important to have the proper order of the sheets in your Excel file. 
        The order should proceed as plate 1, peptide 1 --> plate 1, peptide 2 --> plate 2, peptdide 1 etc.
        """
    )

    parser.add_argument(
        "path", type=str, help="Path to the excel file with sheets for each plate."
    )
    parser.add_argument(
        "number of plates",
        type=int,
        help="The total number of deepwell plates that were tested.",
    )
    parser.add_argument(
        "peptide1", type=str, help="Name of the first peptide.", default="peptide1"
    )
    parser.add_argument(
        "selectivity threshold 1",
        type=float,
        help="Standard Deviation selectivity threshold for the first peptide.",
    )
    parser.add_argument(
        "peptide2", type=str, help="Name of the second peptide.", default="peptide2"
    )
    parser.add_argument(
        "selectivity threshold 2",
        type=float,
        help="Standard Deviation selectivity threshold for the second peptide.",
    )
    parser.add_argument(
        "--controls",
        metavar=str,
        nargs="+",
        type=str,
        help="List of positive control wells.",
        default=["A1", "A2", "A3", "H10", "H11", "H12"],
    )
    parser.add_argument(
        "--byplate",
        type=bool,
        default=True,
        help="Whether to find hits using all positive wells, or only based on the positives that reside on the same plate.",
    )

    args = vars(parser.parse_args())

    plate_order = []
    for plate in range(args["number of plates"]):
        plate_order.append(plate + 1)
        plate_order.append(plate + 1)

    print("Expected plate order:", plate_order)

    peptide_order = [args["peptide1"], args["peptide2"]] * args["number of plates"]

    print("Expected peptide order:", peptide_order)

    data = importPlates(args["path"], plate_order, peptide_order)
    # positive_wells = ["A1", "A2", "A3", "H10", "H11", "H12"]
    print("Got", data.shape[0], "datapoints.")
    assignControls(data, args["controls"])
    print("Control wells assigned: ", args["controls"])
    pp = pivotPlates(data)
    print(pp)
    ratios = computeRatios(pp)

    if args["byplate"]:
        a = find_hits_by_plate(
            ratios, -1, args["selectivity threshold 1"], args["selectivity threshold 2"]
        )
    else:
        a = find_hits(
            ratios, -1, args["selectivity threshold 1"], args["selectivity threshold 2"]
        )

    export_to_pick(a, args["peptide1"], args["peptide2"], args["path"])


if __name__ == "__main__":
    exit(main())
