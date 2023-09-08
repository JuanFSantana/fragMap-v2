import json
import subprocess
import os
import argparse
import concurrent.futures
import sys
import uuid
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd


def check_regions_bed(regions_bed_path: str) -> Tuple[int, int]:
    """
    Checks if the bed file follows certain conditions.

    The conditions checked are:
        1. The bed file should not have more than 6 columns.
        2. For each row, the second column's value should be greater than the first column's value.
        3. All regions should be of the same size and even.
        4. Column 3 should have unique labels.

    Returns
    -------
    int start of region relative to center.
    int end of region relative to center.
    """

    df = pd.read_csv(regions_bed_path, sep="\t", header=None)
    if df.shape[1] > 6:
        return sys.exit("Your bed file has more than 6 columns. Exiting")
    elif any(df[2] - df[1] < 0):
        return sys.exit(
            "Coordinate in column 3 should be greater than coordinate in column 2 (column indexing starts at 1, not 0). Exiting"
        )
    elif len(set(df[2] - df[1])) > 1:
        return sys.exit("All regions should be of the same size. Exiting")
    elif list(set(df[2] - df[1]))[0] % 2 != 0:
        return sys.exit("All regions should be of even length. Exiting")
    elif len(set(df[3])) != len(df[3]):
        return sys.exit(
            "All regions should have a unique identifier in column 4 (column indexing starts at 1, not 0). Exiting"
        )

    center = int((df[1].iloc[0] + df[2].iloc[0]) / 2)
    start = int(df[1].iloc[0] - center)
    end = int(df[2].iloc[0] - center)

    return start, end


def run_featureCount(args: Tuple[str, Dict]) -> str:
    """
    :param args: tuple of key and value from the dictionary
    :return: path to the temporary bed file with feature counts

    The function runs bedtools intersect and featureCount
    """
    # key is in the format of (path_to_regions_bed, path_to_reads_bed, spikein, overlap_type)
    # value is in the format of {name: [[position_left, position_right, size_left, size_right]]}
    key, value = args
    path_to_regions_bed, path_to_reads_bed, spikein, overlap_type = key
    json_data = json.dumps(value)

    # create a random name for the temporary bed file
    random_name_1, random_name_2 = str(uuid.uuid4()), str(uuid.uuid4())
    temp_data_bedtools, temp_data = Path(Path.cwd(), random_name_1 + ".bed"), Path(
        Path.cwd(), random_name_2 + ".bed"
    )

    # run bedtools
    cmd = f"bedtools intersect -a {path_to_regions_bed} -b {path_to_reads_bed} -wa -wb > {temp_data_bedtools}"
    completed_process = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)

    if completed_process.returncode != 0:
        error_message = completed_process.stderr.decode("utf-8")
        print(f"{error_message}")
        sys.exit(1)

    # call featureCount
    if os.name == "posix":  # Linux or macOS
        featureCount = "./featureCount"
    elif os.name == "nt":  # Windows
        featureCount = "./featureCount.exe"

    subprocess.run(
        [
            featureCount,
            temp_data_bedtools,
            json_data,
            temp_data,
            overlap_type,
            str(spikein),
        ]
    )

    # remove the temporary bed file
    os.remove(temp_data_bedtools)

    return temp_data


def consolidate_files(results: List[str], output_directory: str) -> None:
    """
    :param results: list of paths to the temporary bed files with counts
    :param output_directory: path to the output directory
    :return: None

    The function consolidates the temporary bed files with counts into one bed file
    """

    files = [open(file, "r") for file in results]
    with open(output_directory + "feature_counts.txt", "w") as output_file:
        for lines in zip(*files):
            for line in lines:
                joined = [(line.split(), idx) for idx, line in enumerate(lines)]
                final_line = ""
                for each_line, idx in joined:
                    if idx == 0:
                        for each in each_line:
                            final_line += each + "\t"
                    if idx > 0:
                        for each in each_line[1:]:
                            final_line += each + "\t"
            output_file.write(final_line + "\n")

    for f in files:
        f.close()

    # remove the temporary bed files
    for path in results:
        os.remove(path)


def parse_args():
    """
    Get arguments and make multiple checks
    """

    parser = argparse.ArgumentParser(
        prog="dff_feature_count.py",
        description="The program tallies the number of DFF-Seq reads with specific length ranges that overlap within a designated genomic interval",
    )
    parser.add_argument(
        "regions",
        type=str,
        help="Bed file of genomic regions of chosen length. The regions should be of even length and the MaxTSS should be in the middle of the region.",
    )
    parser.add_argument(
        "-f",
        dest="fragments_range_size_position",
        metavar="\b",
        nargs="*",
        required=True,
        action="append",
        help="Singular fragment dataset, followed by position range, followed by fragment range\
                        Examples: -f /home/reads.bed 20 1000 400 800 \
                                 -f /home/reads.bed 20 1000 400 800 20 1000 300 600",
    )

    parser.add_argument(
        "-n",
        dest="names",
        metavar="\b",
        type=str,
        required=True,
        nargs="*",
        help="Provide a name for each feature (space sperated). The names should be in the same order as the features provided with -f.",
    )
    parser.add_argument(
        "-t",
        dest="overlap_type",
        choices=["centers", "full", "partial"],
        required=True,
        help="Type of overlap: center, full, or partial",
    )
    parser.add_argument(
        "-o",
        dest="output_dir",
        metavar="\b",
        type=str,
        required=True,
        nargs=1,
        help="Path to output, for example -o /home/user/dir",
    )
    parser.add_argument(
        "-s",
        dest="spikein",
        metavar="\b",
        type=float,
        required=False,
        nargs="*",
        default=[1.0],  # Setting the default value to 1.0
        help="Correction factors - must be 1 per dataset (-f) space separated. The correction factors should be in the same order as the datasets provided with -f.",
    )

    args = parser.parse_args()
    regions = args.regions
    reads_sizes_positions = args.fragments_range_size_position
    output_directory = args.output_dir[0]
    names = args.names
    spikein = args.spikein
    overlap_type = args.overlap_type

    # checks
    # check that names are unique
    if len(names) != len(set(names)):
        sys.exit("The names must be unique")
    # check if the regions file exists
    if not os.path.isfile(regions):
        sys.exit("The regions file does not exist")
    if not os.path.isfile(reads_sizes_positions[0][0]):
        sys.exit("The reads file does not exist")
    # check if the output directory exists
    if not os.path.isdir(output_directory):
        sys.exit("The output directory does not exist")
    if not output_directory.endswith("/"):
        output_directory = output_directory + "/"
    # spike-in values must be 1 per dataset
    if len(spikein) != len(reads_sizes_positions):
        sys.exit(
            "The number of spike-in values must be equal to the number of datasets"
        )

    to_analyze_dict = {}
    count = 0
    for idx, dataset in enumerate(reads_sizes_positions):
        reads_bed_path = dataset[0]
        features = dataset[1:]

        if len(features) % 4 != 0:
            sys.exit("Each length and positions must have two values")
        for feature in range(0, len(features), 4):
            position_left, position_right, size_left, size_right = (
                int(features[feature]),
                int(features[feature + 1]),
                int(features[feature + 2]),
                int(features[feature + 3]),
            )
            # add to dictionary
            key = (regions, reads_bed_path, spikein[idx], overlap_type)
            sub_key = names[count]
            if key not in to_analyze_dict:
                to_analyze_dict[key] = {sub_key: []}

            if sub_key not in to_analyze_dict[key]:
                to_analyze_dict[key][sub_key] = []

            to_analyze_dict[key][sub_key].extend(
                [[position_left, position_right], [size_left, size_right]]
            )
            count += 1

    total_features = int(sum([len(x[1:]) / 4 for x in reads_sizes_positions]))
    if total_features != len(names):
        sys.exit("The number of names must be equal to the number of features")

    return (
        regions,
        to_analyze_dict,
        output_directory,
    )


if __name__ == "__main__":
    (
        regions,
        to_analyze_dict,
        output_directory,
    ) = parse_args()

    # check bed file, and get start and end of regions relative to center
    region_start, region_end = check_regions_bed(regions)

    # iterate over the dictionary, dictionry will be in the form:
    # {(path_to_bed, spikein): {name: [[position_left, position_right, size_left, size_right]]}}
    final_dict = {}
    data = [(key, value) for key, value in to_analyze_dict.items()]

    # run featureCount
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = list(executor.map(run_featureCount, data))

    # merge files
    consolidate_files(results, output_directory)
