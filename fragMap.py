import argparse
import concurrent.futures
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
from PIL import Image

Image.MAX_IMAGE_PIXELS = 1000000000
import numpy as np
import pandas as pd


def check_regions_bed(regions_bed_path) -> tuple:
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


def fragMap_matrix(arg) -> str:
    """
    :param args: tuple
    :return: (str,str) identifier and path to the temporary bed file with matrix

    The function runs bedtools intersect and fragMapMatrix
    """
    (
        reads_path,
        corr_factor,
        name,
        overlap_type,
        region_start,
        region_end,
        size_left,
        size_right,
        path_to_regions_bed,
    ) = arg
    # create temporary path files
    temp_data_bedtools = Path(Path.cwd(), name + ".bed")

    # run bedtools
    cmd = f"bedtools intersect -a {path_to_regions_bed} -b {reads_path} -wa -wb > {temp_data_bedtools}"
    completed_process = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)

    if completed_process.returncode != 0:
        error_message = completed_process.stderr.decode("utf-8")
        print(f"{error_message}")
        sys.exit(1)

    # call featureCount
    if os.name == "posix":  # Linux or macOS
        fragMapMatrix = "./fragMapMatrix"
    elif os.name == "nt":  # Windows
        fragMapMatrix = "./fragMapMatrix.exe"

    json_data = json.dumps({str(temp_data_bedtools): corr_factor})

    # run fragMapMatrix
    subprocess.run(
        [
            fragMapMatrix,
            json_data,
            overlap_type,
            name,
            str(region_start),
            str(region_end),
            str(size_left),
            str(size_right),
        ]
    )

    return name, temp_data_bedtools


def modifiy_matrix(args2) -> tuple:
    """
    :param args: tuple
    :return: (str, numpy.ndarray)

    Calculates the vertical and horizontal lines per base pair
    """
    name_and_path, height, width = args2
    table_name, path_to_matrix = name_and_path
    df = pd.read_csv(path_to_matrix, sep="\t", index_col=0)
    # Aspect of height and width
    if height >= 1 and width == 1:
        vertically_repeated = df.reindex(df.index.repeat(height))
        final_matrix = vertically_repeated

    elif height >= 1 and width < 1:
        # average first
        # pixels per base
        width_rolling_avg = int(1 / width)
        # rolling average and select rows containing the average window HORIZONTALLY
        df_matrix_width_avg = (
            df.rolling(width_rolling_avg, axis=1).mean().dropna(axis=1, how="any")
        )
        avg_matrix = df_matrix_width_avg[
            df_matrix_width_avg.columns[::width_rolling_avg]
        ]
        # repeat array vertically
        vertically_repeated = avg_matrix.reindex(avg_matrix.index.repeat(height))
        final_matrix = vertically_repeated

    return table_name, final_matrix.to_numpy()


def image(tentnuple):
    """
    :param tentnuple: tuple
    :return: None

    Creates the image
    """
    (
        label,
        matrix,
        black_val,
        size_left,
        size_right,
        output_directory,
        height_modified,
        width_modified,
        gamma,
    ) = tentnuple

    if black_val == "default":
        black_val = int(np.amax(matrix))
    else:
        black_val = int(black_val)

    plt.rcParams["font.size"] = "5"
    plt.rcParams["figure.facecolor"] = "white"

    fig, ax = plt.subplots(dpi=1200)
    im = ax.imshow(matrix, vmin=0, vmax=black_val, cmap="binary")
    ax.tick_params(direction="out", length=1.8, width=0.3)

    matrix_height, matrix_length = matrix.shape
    steps = int(matrix_length / 10)
    real_xcoor = [i for i in range(0, matrix_length + 1, int(steps))]

    def conversion_y_fragmap(num_list):
        true_yticks = [(i - int(size_left)) * int(height_modified) for i in num_list]
        return true_yticks

    def get_ylabels_fragmap():
        if int(size_right) - int(size_left) <= 500:
            ylabels = [
                i for i in range(int(size_left), int(size_right) + 1) if i % 50 == 0
            ]
        else:
            ylabels = [
                i for i in range(int(size_left), int(size_right) + 1) if i % 100 == 0
            ]
        return ylabels

    def x_conver(matrix_length, width):
        if width_modified <= 1:
            xlabels_ = [
                i
                for i in range(
                    int(-matrix_length / (2 * float(width_modified))),
                    int(matrix_length / (2 * float(width_modified))) + 1,
                    int(steps * (1 / float(width_modified))),
                )
            ]
        else:
            xlabels_ = [
                i
                for i in range(
                    int(-matrix_length / 2 / float(width_modified)),
                    int(matrix_length / 2 / float(width_modified) + 1),
                    int(steps / float(width_modified)),
                )
            ]
        return xlabels_

    x_axis_converted_nums = x_conver(matrix_length, width_modified)
    sorted_x_axis_converted_nums = sorted(x_axis_converted_nums)
    second_min = abs(sorted_x_axis_converted_nums[1])
    if second_min < 1000:
        xlabels = [i if i != 0 else 1 for i in x_axis_converted_nums]
    else:
        xlabels = [
            str(int(i / 1000)) + "k" if i != 0 else 1 for i in x_axis_converted_nums
        ]

    plt.xticks(real_xcoor, xlabels)
    ylabels = get_ylabels_fragmap()
    plt.yticks(conversion_y_fragmap(ylabels), ylabels)

    if gamma != 1:
        image_path = Path(
            output_directory,
            "-".join(
                [
                    str(label),
                    "Max",
                    str(black_val),
                    "X",
                    str(width_modified),
                    "Y",
                    str(height_modified),
                    "Gamma",
                    str(gamma),
                ]
            )
            + ".png",
        )
    else:
        image_path = Path(
            output_directory,
            "-".join(
                [
                    str(label),
                    "Max",
                    str(black_val),
                    "X",
                    str(width_modified),
                    "Y",
                    str(height_modified),
                ]
            )
            + ".png",
        )

    ax = fig.gca()

    for axis in ["top", "bottom", "left", "right"]:
        ax.spines[axis].set_linewidth(0.2)
    ax.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.tick_params(which="minor", length=0.8, width=0.3)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.05)
    cbar = plt.colorbar(im, cax=cax, orientation="vertical")
    cbar.ax.tick_params(size=1, width=0.3)
    cbar.outline.set_linewidth(0.1)

    plt.savefig(image_path, format="png", facecolor="w", bbox_inches="tight", dpi=1200)

    def gammma(x, r):
        """
        From: https://linuxtut.com/en/c3cd5475663c41d9f154/
        Gamma correction y=255*(x/255)
        x Input image
        r Gamma correction coefficient
        """
        y = x / 255
        y = y ** (1 / r)

        return np.uint8(255 * y)

    if gamma != 1.0:
        img = Image.open(image_path).convert("RGB")
        arr = np.asarray(img, dtype=float)
        img_gamma = gammma(arr, gamma)
        plt.imsave(image_path, img_gamma)

    plt.close()


def parse_args():
    """
    Gets arguments and makes multiple checks
    """

    parser = argparse.ArgumentParser(
        prog="fragmap.py",
        description="Generates fragMaps from specific range of fragment sizes over a chosen genomic interval. Multiple replicates can be compared in a single fragMap\
            A combined fragMap is automatically generated by summing replicates",
    )
    parser.add_argument(
        "regions", type=str, help="Bed file of genomic regions of chosen length"
    )
    parser.add_argument(
        "-f",
        dest="fragments",
        metavar="\b",
        type=str,
        required=True,
        nargs="*",
        help="Bed file of reads/fragments",
    )
    parser.add_argument(
        "-r",
        dest="range",
        metavar="\b",
        type=int,
        nargs=2,
        required=True,
        help="Range of fragment sizes, for exmaple -r 20 400",
    )
    parser.add_argument(
        "-s",
        dest="spikein",
        metavar="\b",
        type=float,
        nargs="*",
        required=True,
        help="Spike-in or correction factors",
    )
    parser.add_argument(
        "-b",
        dest="black",
        metavar="\b",
        default="default",
        help="Sets the chosen value as black, default is largest number in the matrix",
    )
    parser.add_argument(
        "-c",
        dest="centers",
        action="store_true",
        default=False,
        help="If argument is invoked, the output will be a fragMap of centers of fragments",
    )
    parser.add_argument(
        "-y",
        dest="y_axis",
        metavar="\b",
        type=int,
        default=1,
        help="Horizontal lines/bp for each fragment length | Can only be greater or equal than 1",
    )
    parser.add_argument(
        "-x",
        dest="x_axis",
        metavar="\b",
        type=float,
        default=1.0,
        help="Vertical lines/bp for each genomic interval displayed, for example -x 1 is one vertical line/bp; -x 0.1 is one vertical line/10 bp | Can't be greater than 1",
    )
    parser.add_argument(
        "-g",
        dest="gamma",
        metavar="\b",
        type=float,
        default=1.0,
        help="Gamma correction",
    )
    parser.add_argument(
        "-o",
        dest="output_dir",
        metavar="\b",
        type=str,
        required=True,
        nargs=1,
        help="Path to output",
    )
    parser.add_argument(
        "-n",
        dest="names",
        metavar="\b",
        type=str,
        required=True,
        nargs="*",
        help="Image output names (no spaces, no file extensions)",
    )

    args = parser.parse_args()

    read_file = args.fragments
    size_left, size_right = args.range
    max_val = args.black
    width = args.x_axis
    height = args.y_axis
    spikeins = args.spikein
    identifier = args.names

    if float(width) > 1:
        sys.exit("Missing -x argument. x must be int or float less than or equal to 1")

    if float(height) < 1:
        sys.exit(
            "Missing -y argument. y must be int or float greater than or equal to 1"
        )

    if len(read_file) != len(spikeins) and len(read_file) != len(identifier):
        sys.exit(
            "The number of bed files, spike-ins or image output names do not match"
        )

    if size_left > size_right:
        sys.exit("Fragment size range is incorrect")

    return args


def main(args):
    regions_to_analyze = args.regions
    read_file = args.fragments
    max_val = args.black
    height = args.y_axis
    width = args.x_axis
    gamma = args.gamma
    output_directory = args.output_dir[0]
    identifier = args.names
    size_left, size_right = args.range
    spikeins = args.spikein
    centers = args.centers

    # check bed file, and get start and end of regions relative to center
    region_start, region_end = check_regions_bed(regions_to_analyze)

    # check if output directory exists
    if not Path(output_directory).exists():
        sys.exit("Output directory does not exist")

    # center or full fragment
    if centers:
        overlap_type = "centers"
    else:
        overlap_type = "full"

    # create iterable for multiprocessing

    info_iterable = [
        (
            read,
            corr,
            name,
            overlap_type,
            str(region_start),
            str(region_end),
            str(size_left),
            str(size_right),
            regions_to_analyze,
        )
        for read, corr, name in zip(read_file, spikeins, identifier)
    ]

    # run fragMap-matrix.cpp
    with concurrent.futures.ProcessPoolExecutor() as executor:
        # make fragMap matrix
        # path_to_matrix is a list of tuples (name - without file extension, path_to_matrix)
        path_to_matrix = list(executor.map(fragMap_matrix, info_iterable))
        # repeat/average fragMap matrix
        args2 = [(matrix_path, height, width) for matrix_path in path_to_matrix]
        # modified matrix and names
        # modified_matrix is a list of tuples (name - without file extension, numpy 2d array)
        modified_matrix = list(executor.map(modifiy_matrix, args2))
        # if more than one rep, sum matrices
        if len(modified_matrix) > 1:
            # sum matrices
            final_matrix = np.sum([matrix[1] for matrix in modified_matrix], axis=0)
            # get name
            names = "_".join([name for name in identifier]) + "_combined"
            # add to list
            modified_matrix.append((names, final_matrix))
        # create image
        arguments = [
            (
                matrix,
                label,
                max_val,
                size_left,
                size_right,
                output_directory,
                height,
                width,
                gamma,
            )
            for matrix, label in modified_matrix
        ]
        results = list(executor.map(image, arguments))

    # remove temporary files
    for matrix in path_to_matrix:
        os.remove(matrix[1])


if __name__ == "__main__":
    args = parse_args()
    main(args)
