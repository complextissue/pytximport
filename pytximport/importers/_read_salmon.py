import gzip
import json
import os
from pathlib import Path
from typing import Literal, Union

import numpy as np

from ..definitions import InferentialReplicates, TranscriptData
from ._read_tsv import read_tsv


def read_inferential_replicates_salmon(
    file_path: Union[str, Path],
    aux_dir_name: Literal["aux_info", "aux"] = "aux_info",
) -> InferentialReplicates:
    """Read inferential replicates from a salmon quantification file.

    Args:
        file_path (Union[str, Path]): The path to the quantification file.
        aux_dir_name (Literal["aux_info", "aux"], optional): The name of the aux directory. Defaults to "aux_info".

    Returns:
        InferentialReplicates: The inferential replicates.
    """
    if not isinstance(file_path, Path):
        file_path = Path(file_path)

    if not file_path.is_dir():
        file_path = file_path.parent

    cmd_info_path = file_path / "cmd_info.json"
    if not os.path.exists(cmd_info_path):
        raise ImportError("cmd_info.json not found.")

    with open(cmd_info_path, "r") as f:
        cmd_info = json.load(f)

    if "auxDir" in cmd_info:
        aux_dir_name = cmd_info["auxDir"]

    aux_dir = file_path / aux_dir_name

    if not os.path.exists(aux_dir):
        raise ImportError("Auxiliary directory not found.")

    meta_info_path = aux_dir / "meta_info.json"
    with open(meta_info_path) as f:
        meta_info = json.load(f)

    if "salmon_version" in meta_info:
        assert meta_info["salmon_version"] >= "0.8.0", "Salmon version must be >= 0.8.0 to use inferential replicates."
    if "sailfish_version" in meta_info:
        assert (
            meta_info["sailfish_version"] >= "0.9.0"
        ), "Sailfish version must be >= 0.9.0 to use inferential replicates."

    bootstrap_count = meta_info.get("num_bootstraps", 0)

    if bootstrap_count == 0:
        raise ImportError("No bootstraps found.")

    bootstrap_path = aux_dir / "bootstrap" / "bootstraps.gz"

    if not os.path.exists(bootstrap_path):
        raise ImportError("Bootstraps file not found.")

    if "num_valid_targets" in meta_info:
        meta_info["num_targets"] = meta_info["num_valid_targets"]

    target_count = meta_info.get("num_targets", 0)

    if target_count == 0:
        raise ImportError("No inferential replicate targets found.")

    expected_n = target_count * bootstrap_count

    try:
        # Try to read as floats
        with gzip.open(bootstrap_path, "rb") as f:
            bootstrap_data = np.frombuffer(f.read(), dtype=np.float64, count=expected_n)
        assert len(bootstrap_data) == expected_n
    except (AssertionError, ValueError):
        # Try to read as integers
        with gzip.open(bootstrap_path, "rb") as f:
            bootstrap_data = np.frombuffer(f.read(), dtype=np.int32, count=expected_n)

    bootstrap_data = bootstrap_data.reshape((bootstrap_count, target_count)).T

    return InferentialReplicates(
        variance=np.var(bootstrap_data, axis=1, ddof=1),
        replicates=bootstrap_data,
    )


def read_salmon(
    file_path: Union[str, Path],
    id_column: str = "Name",
    counts_column: str = "NumReads",
    length_column: str = "EffectiveLength",
    abundance_column: str = "TPM",
    aux_dir_name: Literal["aux_info", "aux"] = "aux_info",
    inferential_replicates: bool = False,
    recompute_counts: bool = False,
) -> TranscriptData:
    """Read a salmon quantification file.

    Args:
        file_path (Union[str, Path]): The path to the quantification file.
        id_column (str, optional): The column name for the transcript id. Defaults to "Name".
        counts_column (str, optional): The column name for the counts. Defaults to "NumReads".
        length_column (str, optional): The column name for the length. Defaults to "EffectiveLength".
        abundance_column (str, optional): The column name for the abundance. Defaults to "TPM".
        aux_dir_name (Literal["aux_info", "aux"], optional): The name of the aux directory. Defaults to "aux_info".
        inferential_replicates (bool, optional): Whether to read inferential replicates. Defaults to False.
        recompute_counts (bool, optional): Whether inferential replicates will be used to recompute counts and
            abundances. If true, the counts and abundances will not be read from the file. Defaults to False.

    Returns:
        TranscriptData: The transcript-level expression.
    """
    if not isinstance(file_path, Path):
        file_path = Path(file_path)

    if file_path.is_dir():
        # Add quant.sf to the file path
        file_path = file_path / "quant.sf"

    # Check that we are importing a .sf file
    if not file_path.suffix == ".sf" and not file_path.suffix == ".gz":
        raise ImportError("Only .sf and .gz files are supported.")

    transcript_data = read_tsv(
        file_path,
        id_column=id_column,
        counts_column=counts_column,
        length_column=length_column,
        abundance_column=abundance_column,
        recompute_counts=recompute_counts,
    )

    if inferential_replicates:
        transcript_data["inferential_replicates"] = read_inferential_replicates_salmon(
            file_path,
            aux_dir_name=aux_dir_name,
        )

    return transcript_data
