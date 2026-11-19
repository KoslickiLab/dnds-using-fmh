#!/usr/bin/env python3
"""Check that --mode reference reports exactly what --mode main reports.

    usage:  python tests/check_reference_mode.py --fasta_input_list genomes.csv --reference NAME \
                --directory OUT [--scaled 1] [--downsample_to 5,10] [--cores 8] [--threshold 0.05]

1. main mode sketches the genomes at --scaled and compares every pair (OUT/main).
2. reference mode compares every pair from the same sketches (OUT/reference/reference_scaled<S>).
3. For each scale factor in --downsample_to, main mode sketches the genomes again at that scale
   factor (OUT/main_scaled<S>), and reference mode downsamples the sketches of step 1 to it.
With --sketch_shards N, reference mode makes its own sketches instead, in N shards (its
--sketch_shards option), so that the sharded sketching is checked against main mode too.
With --reuse_main, main-mode runs whose fmh_omega_<ksize>.csv already exists are not repeated.
Every pair, and every column of fmh_omega_<ksize>.csv (DNA_max_Cfrac, AA_max_Cfrac, PdN, PdS,
PdN/PdS, dN/dS), must be identical. Exits with status 1 on any difference.
"""
import argparse
import os
import subprocess
import sys

import numpy as np
import pandas as pd

TOOL = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "fmh_dnds.py")
COLUMNS = ["DNA_max_Cfrac", "AA_max_Cfrac", "PdN", "PdS", "PdN/PdS", "dN/dS"]


def run(*arguments):
    """Run the tool from the folder of the input csv, so that relative paths in it resolve."""
    command = [sys.executable, TOOL, *[str(a) for a in arguments]]
    result = subprocess.run(command, capture_output=True, text=True, cwd=WORKING_FOLDER)
    if result.returncode != 0:
        sys.exit(f"failed: {' '.join(command)}\n{result.stderr[-3000:]}")


def keyed(table):
    table = table.copy()
    table["pair"] = [" | ".join(sorted((a, b))) for a, b in zip(table.query_name.astype(str),
                                                                  table.match_name.astype(str))]
    return table.set_index("pair")[COLUMNS].sort_index()


def compare(main_csv, reference_folder, label):
    main = keyed(pd.read_csv(main_csv, index_col=0))
    reference = keyed(pd.read_parquet(f"{reference_folder}/pairs"))
    problems = []
    if not main.index.equals(reference.index):
        problems.append(f"pairs differ: {len(main)} in main, {len(reference)} in reference, "
                        f"{len(main.index.symmetric_difference(reference.index))} not in both")
    else:
        for column in COLUMNS:
            a, b = main[column].to_numpy(), reference[column].to_numpy()
            same = (a == b) | (np.isnan(a) & np.isnan(b))
            if not same.all():
                worst = np.nanmax(np.abs(a - b)[~same]) if (~same).any() else 0
                problems.append(f"{column}: {int((~same).sum())} of {len(a)} values differ (largest {worst:.3g})")
    print(f"{label}: {len(main)} pairs, " + ("identical in every column" if not problems else "; ".join(problems)))
    return not problems


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta_input_list", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--directory", required=True)
    parser.add_argument("--scaled", type=int, default=1)
    parser.add_argument("--downsample_to", default="")
    parser.add_argument("--cores", type=int, default=8)
    parser.add_argument("--threshold", type=float, default=0.05)
    parser.add_argument("--ksize", type=int, default=7)
    parser.add_argument("--sketch_shards", type=int, default=1)
    parser.add_argument("--reuse_main", action="store_true")
    arguments = parser.parse_args()
    out = os.path.abspath(arguments.directory)
    genomes = os.path.abspath(arguments.fasta_input_list)
    WORKING_FOLDER = os.path.dirname(genomes)
    common = ["--ksize", arguments.ksize, "--cores", arguments.cores, "--threshold", arguments.threshold]
    scales = [int(s) for s in arguments.downsample_to.split(",") if s]

    def run_main(folder, scale):
        os.makedirs(folder, exist_ok=True)
        if not (arguments.reuse_main and os.path.exists(f"{folder}/fmh_omega_{arguments.ksize}.csv")):
            run("--fasta_input_list", genomes, "--scaled_input", scale, "--directory", folder, *common)

    run_main(f"{out}/main", arguments.scaled)
    reference_arguments = ["--mode", "reference", "--reference", arguments.reference,
                           "--scaled_input", arguments.scaled, "--directory", f"{out}/reference", *common]
    if arguments.sketch_shards > 1:
        reference_arguments += ["--fasta_input_list", genomes, "--sketch_shards", arguments.sketch_shards]
    else:
        reference_arguments += ["--sketches", f"{out}/main/data.zip"]
    if scales:
        reference_arguments += ["--downsample_to", ",".join(str(s) for s in scales)]
    os.makedirs(f"{out}/reference", exist_ok=True)
    run(*reference_arguments)

    ok = compare(f"{out}/main/fmh_omega_{arguments.ksize}.csv",
                 f"{out}/reference/reference_scaled{arguments.scaled}", f"scale factor {arguments.scaled}")
    for scale in scales:
        run_main(f"{out}/main_scaled{scale}", scale)
        ok &= compare(f"{out}/main_scaled{scale}/fmh_omega_{arguments.ksize}.csv",
                      f"{out}/reference/reference_scaled{scale}",
                      f"scale factor {scale} (main sketched at {scale}; reference downsampled from {arguments.scaled})")
    sys.exit(0 if ok else 1)
