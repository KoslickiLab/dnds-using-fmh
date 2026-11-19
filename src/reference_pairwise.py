"""All-pairs FracMinHash dN/dS for a collection of closely related genomes (--mode reference).

The estimate is exactly the one `main` mode reports. Sketches are made by the same
`sourmash scripts manysketch` call (sourmash_ext.run_manysketch); containment is max containment of
the nucleotide (k = 3 * ksize) and protein (k = ksize) sketches; and PdN, PdS, PdN/PdS and dN/dS
come from the functions in src/dnds.py, which this module imports and does not change. The only
difference is how the containments of all pairs are computed.

`main` mode hands the sketches to `sourmash scripts pairwise`, which intersects the two complete
sketches of every pair. For genomes that are nearly identical that repeats the same work for every
pair, because almost every hash is shared by almost every genome. This mode compares each sketch
once with the sketch of one chosen genome, the reference R, and keeps only

    lacks_i = hashes of R that sketch S_i does not contain
    adds_i  = hashes of S_i that R does not contain

Then, for any two sketches, exactly (lacks are hashes of R, adds are not):

    |S_i & S_j| = |R| - |lacks_i| - |lacks_j| + |lacks_i & lacks_j| + |adds_i & adds_j|
    |S_i|       = |R| - |lacks_i| + |adds_i|

The two "&" counts of every pair come from one sparse matrix product over the lacks and adds of all
sketches. The result is the same integer intersection as sourmash computes, so the containments,
and everything derived from them, are the same numbers, not approximations. The identity holds
whatever the genomes; it is only faster than `main` when they are close to the reference.

A sketch with scale factor S keeps the hashes at or below sourmash's max_hash for S, so a coarser
scale factor than the sketches were made with (--downsample_to) is obtained by dropping the lacks,
adds and reference hashes above that value; this gives the same sets as sourmash's own downsample().

Output, one folder per scale factor:
    <directory>/reference_scaled<S>/pairs/part-<block>.parquet   one row per pair, the columns of
        `main` mode's fmh_omega_<ksize>.csv: query_name, match_name, DNA_max_Cfrac, AA_max_Cfrac,
        ksize, PdN, PdS, PdN/PdS, dN/dS (read the folder with pandas.read_parquet)
    <directory>/reference_scaled<S>/run.json               sketches, pairs, timings
    <directory>/reference_scaled<S>/names.txt              the genomes compared, one per line
    <directory>/reference_scaled<S>/fmh_omega_<ksize>.csv  the same rows as one CSV (--write_csv)
"""
import csv
import json
import os
import time
import zipfile
from multiprocessing import Pool

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from loguru import logger
from scipy import sparse
from sourmash.minhash import _get_max_hash_for_scaled
from sourmash.signature import load_signatures_from_json

from src import dnds

PAIRS_PER_BLOCK = 5_000_000


# ------------------------------------------------------------------ sketching in shards

def sketch_in_shards(fasta_input_list, ksize, scaled, cores, directory, shards):
    """Sketch exactly as `main` does (sourmash_ext.run_manysketch), but with the input csv split into
    shards that are sketched side by side, each into its own zip. The sketches are the same; what
    changes is that writing one zip file, which manysketch does from a single thread, is no longer
    the bottleneck when there are tens of thousands of large sketches. Returns the zip paths."""
    from concurrent.futures import ThreadPoolExecutor
    from src import sourmash_ext
    lines = [line for line in open(fasta_input_list).read().splitlines() if line.strip()]
    header, genomes = lines[0], lines[1:]
    shards = min(shards, len(genomes))
    folders = []
    for shard in range(shards):
        folder = f"{directory}/sketch_shards/{shard:04d}"
        os.makedirs(folder, exist_ok=True)
        open(f"{folder}/genomes.csv", "w").write("\n".join([header] + genomes[shard::shards]) + "\n")
        folders.append(folder)
    running = min(shards, cores)
    threads = max(1, cores // running)

    def sketch(folder):
        sourmash_ext.run_manysketch(fasta_file_csv=f"{folder}/genomes.csv", ksize=ksize, scaled=scaled,
                                    cores=threads, working_dir=folder)

    with ThreadPoolExecutor(max_workers=running) as pool:
        list(pool.map(sketch, folders))
    return [f"{folder}/data.zip" for folder in folders]


# ------------------------------------------------------------------ reading sketches

def read_manifest(zip_path):
    """The rows of the sourmash manifest of a sketch zip, read with Python's zipfile. A zip of more
    than 4 GB is written in the ZIP64 format, which sourmash's own zip reader (ZipStorage) rejects
    ("Extra data field contains disk number") while Python's reads it."""
    with zipfile.ZipFile(zip_path) as archive:
        text = archive.read("SOURMASH-MANIFEST.csv").decode()
    rows = list(csv.DictReader(line for line in text.splitlines() if not line.startswith("#")))
    for row in rows:
        row["ksize"], row["scaled"] = int(row["ksize"]), int(row["scaled"])
    return rows


def sketch_locations(zip_paths, ksize, names=None):
    """{(name, "DNA" | "protein"): (zip path, location in the zip, scaled)} for every sketch at the
    two k-mer sizes, optionally only for the given names."""
    wanted = {"DNA": 3 * ksize, "protein": ksize}
    locations = {}
    for zip_path in zip_paths:
        for row in read_manifest(zip_path):
            if row["moltype"] not in wanted or row["ksize"] != wanted[row["moltype"]]:
                continue
            if names is not None and row["name"] not in names:
                continue
            key = (row["name"], row["moltype"])
            if key in locations:
                raise ValueError(f"more than one {row['moltype']} sketch named {row['name']!r}")
            locations[key] = (zip_path, row["internal_location"], row["scaled"])
    return locations


_open_zips = {}


def read_hashes(zip_path, location):
    """The hashes of one sketch, as a sorted uint64 array."""
    if zip_path not in _open_zips:
        _open_zips[zip_path] = zipfile.ZipFile(zip_path)
    signature = next(iter(load_signatures_from_json(_open_zips[zip_path].read(location))))
    minhash = signature.minhash
    return np.fromiter(minhash.hashes.keys(), dtype=np.uint64, count=len(minhash))


def _set_reference(reference):
    global REFERENCE_HASHES
    REFERENCE_HASHES = reference
    # a forked worker must not share the parent's open zip files (one file offset between processes)
    _open_zips.clear()


def _difference(task):
    index, zip_path, location = task
    hashes = read_hashes(zip_path, location)
    lacks = np.setdiff1d(REFERENCE_HASHES, hashes, assume_unique=True)
    adds = np.setdiff1d(hashes, REFERENCE_HASHES, assume_unique=True)
    return index, lacks, adds


def difference_matrix(names, locations, reference_name, moltype, cores):
    """Sparse 0/1 matrix (sketches x hashes) of the lacks and adds of every sketch, the sorted hashes
    the columns stand for, which of those columns are reference hashes, and the reference hashes."""
    zip_path, location, _ = locations[(reference_name, moltype)]
    reference = read_hashes(zip_path, location)
    tasks = [(i, *locations[(name, moltype)][:2]) for i, name in enumerate(names)]
    lacks, adds = [None] * len(names), [None] * len(names)
    with Pool(cores, initializer=_set_reference, initargs=(reference,)) as pool:
        for index, lacked, added in pool.imap_unordered(_difference, tasks, chunksize=4):
            lacks[index], adds[index] = lacked, added
    universe = np.unique(np.concatenate(lacks + adds))
    row_lengths = np.array([len(a) + len(b) for a, b in zip(lacks, adds)], dtype=np.int64)
    indptr = np.zeros(len(names) + 1, dtype=np.int64)
    np.cumsum(row_lengths, out=indptr[1:])
    columns = np.empty(int(indptr[-1]), dtype=np.int32)
    for row, (lacked, added) in enumerate(zip(lacks, adds)):
        columns[indptr[row]:indptr[row + 1]] = np.searchsorted(universe, np.concatenate([lacked, added]))
    matrix = sparse.csr_matrix((np.ones(len(columns), dtype=np.int32), columns, indptr),
                               shape=(len(names), len(universe)))
    matrix.sort_indices()
    position = np.minimum(np.searchsorted(reference, universe), len(reference) - 1)
    in_reference = (reference[position] == universe).astype(np.int32)
    return matrix, universe, in_reference, reference


# ------------------------------------------------------------------ all pairs at one scale factor

def at_scale(matrix, universe, in_reference, reference, max_hash):
    """The difference matrix, its transpose, |R|, |lacks| and |S| after keeping hashes <= max_hash.
    The columns are sorted by hash, so this keeps a leading block of columns."""
    kept = int(np.searchsorted(universe, np.uint64(max_hash), side="right"))
    if kept < matrix.shape[1]:
        matrix = matrix[:, :kept].tocsr()
    reference_size = int(np.searchsorted(reference, np.uint64(max_hash), side="right"))
    lacks = np.asarray(matrix @ in_reference[:kept]).ravel().astype(np.int64)
    adds = np.diff(matrix.indptr).astype(np.int64) - lacks
    return dict(matrix=matrix, transposed=matrix.T.tocsr(), reference_size=reference_size,
                lacks=lacks, size=reference_size - lacks + adds)


def _set_scale(state):
    global STATE
    STATE = state


def _block(start):
    """Every pair (i, j) with start <= i < start + rows and j > i, written as one parquet part."""
    names, ksize, threshold = STATE["names"], STATE["ksize"], STATE["threshold"]
    n = len(names)
    stop = min(start + STATE["rows"], n)
    rows = np.arange(start, stop)
    upper = rows[:, None] < np.arange(n)[None, :]
    first, second = np.nonzero(upper)
    first = rows[first]
    containment = {}
    for moltype in ("DNA", "protein"):
        part = STATE[moltype]
        shared = np.asarray((part["matrix"][start:stop] @ part["transposed"]).todense())[upper]
        intersection = (part["reference_size"] - part["lacks"][first] - part["lacks"][second] + shared)
        containment[moltype] = intersection / np.minimum(part["size"][first], part["size"][second])
    # `main` mode keeps a pair for a molecule when its max containment exceeds the threshold, and
    # joins the two molecules' tables, so a pair kept for one molecule only has no value for the other
    passes = {moltype: values > threshold for moltype, values in containment.items()}
    keep = passes["DNA"] | passes["protein"]
    nt = np.where(passes["DNA"], containment["DNA"], np.nan)[keep]
    aa = np.where(passes["protein"], containment["protein"], np.nan)[keep]
    first, second = first[keep], second[keep]
    with np.errstate(divide="ignore", invalid="ignore"):
        table = pd.DataFrame({"DNA_max_Cfrac": nt, "AA_max_Cfrac": aa})
        table["ksize"] = ksize
        table["PdN"] = dnds.calc_PdN(protein_containment=table["AA_max_Cfrac"], k=table["ksize"])
        table["PdS"] = dnds.calc_PdS(protein_containment=table["AA_max_Cfrac"],
                                     nt_containment=table["DNA_max_Cfrac"], k=table["ksize"])
        table["PdN/PdS"] = dnds.dNdS_ratio(nt_containment=table["DNA_max_Cfrac"],
                                           protein_containment=table["AA_max_Cfrac"], k=table["ksize"])
        table["dN/dS"] = dnds.dNdS_ratio_with_constant(nt_containment=table["DNA_max_Cfrac"],
                                                       protein_containment=table["AA_max_Cfrac"],
                                                       k=table["ksize"])
    dictionary = pa.array(names)
    columns = {"query_name": pa.DictionaryArray.from_arrays(pa.array(first.astype(np.int32)), dictionary),
               "match_name": pa.DictionaryArray.from_arrays(pa.array(second.astype(np.int32)), dictionary)}
    for column in table.columns:
        columns[column] = pa.array(table[column].to_numpy())
    pq.write_table(pa.table(columns), f"{STATE['out']}/pairs/part-{start:09d}.parquet")
    return len(table)


# ------------------------------------------------------------------ the mode

def run(zip_paths, reference_name, ksize, scaled, directory, cores, threshold,
        names_file=None, downsample_to=(), write_csv=False, sketch_seconds=0.0):
    """Compute every pair of the sketches in zip_paths (optionally only the names in names_file) at
    the scale factor `scaled` and at each scale factor in downsample_to."""
    names = None
    if names_file:
        names = {line.strip() for line in open(names_file) if line.strip()} | {reference_name}
    locations = sketch_locations(zip_paths, ksize, names)
    sketch_names = sorted({name for name, _ in locations})
    for name in sketch_names:
        for moltype in ("DNA", "protein"):
            if (name, moltype) not in locations:
                raise ValueError(f"no {moltype} sketch for {name!r}")
    if reference_name not in sketch_names:
        raise ValueError(f"the reference {reference_name!r} is not among the sketches")
    sketch_scaled = {scale for _, _, scale in locations.values()}
    if len(sketch_scaled) != 1:
        raise ValueError(f"the sketches were made with different scale factors: {sorted(sketch_scaled)}")
    sketch_scaled = sketch_scaled.pop()
    scales = sorted({scaled, *downsample_to})
    if scales[0] < sketch_scaled:
        raise ValueError(f"cannot compare at scale factor {scales[0]}: the sketches have {sketch_scaled}")
    logger.info(f"{len(sketch_names)} sketches, reference {reference_name!r}, scale factors {scales}")

    start_time = time.time()
    full = {}
    for moltype in ("DNA", "protein"):
        full[moltype] = difference_matrix(sketch_names, locations, reference_name, moltype, cores)
        logger.info(f"{moltype}: {full[moltype][0].nnz} differences over {len(full[moltype][1])} distinct hashes")
    load_seconds = time.time() - start_time

    outputs = {}
    for scale in scales:
        max_hash = _get_max_hash_for_scaled(scale)
        out = f"{directory}/reference_scaled{scale}"
        os.makedirs(f"{out}/pairs", exist_ok=True)
        for old in os.listdir(f"{out}/pairs"):
            os.remove(f"{out}/pairs/{old}")
        scale_start = time.time()
        state = {moltype: at_scale(*full[moltype], max_hash) for moltype in ("DNA", "protein")}
        state.update(names=sketch_names, ksize=ksize, threshold=threshold, out=out,
                     rows=max(1, PAIRS_PER_BLOCK // len(sketch_names)))
        blocks = list(range(0, len(sketch_names), state["rows"]))
        with Pool(cores, initializer=_set_scale, initargs=(state,)) as pool:
            pairs = sum(pool.imap_unordered(_block, blocks, chunksize=1))
        del state
        seconds = time.time() - scale_start
        report = dict(mode="reference", reference=reference_name, sketches=len(sketch_names),
                      pairs=int(pairs), ksize=ksize, sketch_scaled=sketch_scaled, scaled=scale,
                      threshold=threshold, sketch_seconds=round(sketch_seconds, 1),
                      load_seconds=round(load_seconds, 1),
                      pairs_seconds=round(seconds, 1), cores=cores, sketch_files=list(zip_paths))
        json.dump(report, open(f"{out}/run.json", "w"), indent=1)
        open(f"{out}/names.txt", "w").write("\n".join(sketch_names) + "\n")
        logger.success(f"scale factor {scale}: {pairs} pairs in {seconds:.0f} s, written to {out}")
        if write_csv:
            table = pd.read_parquet(f"{out}/pairs")
            table.to_csv(f"{out}/fmh_omega_{ksize}.csv")
        outputs[scale] = out
    return outputs
