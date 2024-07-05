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
) -> Union[InferentialReplicates, None]:
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
        return None

    with open(cmd_info_path, "r") as f:
        cmd_info = json.load(f)

    if "auxDir" in cmd_info:
        aux_dir_name = cmd_info["auxDir"]

    aux_dir = file_path / aux_dir_name

    if not os.path.exists(aux_dir):
        return None

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
        return None

    bootstrap_path = aux_dir / "bootstrap" / "bootstraps.gz"

    if not os.path.exists(bootstrap_path):
        return None

    if "num_valid_targets" in meta_info:
        meta_info["num_targets"] = meta_info["num_valid_targets"]

    target_count = meta_info.get("num_targets", 0)

    if target_count == 0:
        return None

    expected_n = target_count * bootstrap_count

    try:
        # try to read as floats
        with gzip.open(bootstrap_path, "rb") as f:
            bootstrap_data = np.frombuffer(f.read(), dtype=np.float64, count=expected_n)
        assert len(bootstrap_data) == expected_n
    except (AssertionError, ValueError):
        # try to read as integers
        with gzip.open(bootstrap_path, "rb") as f:
            bootstrap_data = np.frombuffer(f.read(), dtype=np.int32, count=expected_n)

    bootstrap_data = bootstrap_data.reshape((bootstrap_count, target_count)).T
    variance = np.var(bootstrap_data, axis=1, ddof=1)

    return InferentialReplicates(
        variance=variance,
        replicates=bootstrap_data,
    )


def read_salmon(
    file_path: Union[str, Path],
    id_column: str = "Name",
    counts_column: str = "NumReads",
    length_column: str = "EffectiveLength",
    abundance_column: str = "TPM",
    aux_dir_name: Literal["aux_info", "aux"] = "aux_info",
) -> TranscriptData:
    """Read a salmon quantification file.

    Args:
        file_path (Union[str, Path]): The path to the quantification file.

    Returns:
        TranscriptData: The transcript-level expression.
    """
    if not isinstance(file_path, Path):
        file_path = Path(file_path)

    if file_path.is_dir():
        # add quant.sf to the file path
        file_path = file_path / "quant.sf"

    # check that we are importing a .sf file
    if not file_path.suffix == ".sf" and not file_path.suffix == ".gz":
        raise ImportError("Only .sf files are supported.")

    # unzip the file if it is compressed
    if file_path.suffix == ".gz":
        try:
            with gzip.open(file_path, "rt") as f:
                file_content = f.read()
            file_path = file_path.with_suffix(".sf")
            with open(file_path, "w") as f:
                f.write(file_content)
        except Exception as e:
            raise ImportError(f"Could not unzip the file: {file_path}") from e

    transcript_data = read_tsv(
        file_path,
        id_column=id_column,
        counts_column=counts_column,
        length_column=length_column,
        abundance_column=abundance_column,
    )

    transcript_data["inferential_replicates"] = read_inferential_replicates_salmon(
        file_path,
        aux_dir_name=aux_dir_name,
    )

    return transcript_data
