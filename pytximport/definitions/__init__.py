"""Type definitions for internal use in `pytximport`."""

from typing import Any, List, Optional, TypedDict

from numpy.typing import ArrayLike


class InferentialReplicates(TypedDict):
    """Inferential replicates for a set of samples."""

    variance: ArrayLike
    replicates: Any


class OmicData(TypedDict):
    """Omic-level expression data from multiple samples."""

    abundance: Optional[ArrayLike]
    counts: Optional[ArrayLike]
    length: ArrayLike

    inferential_replicates: Optional[InferentialReplicates]


class TranscriptData(OmicData):
    """Transcript-level expression data from multiple samples."""

    transcript_id: List[str]


class GeneData(OmicData):
    """Gene-level expression data from multiple samples."""

    gene_id: List[str]
