from __future__ import annotations

import re
from collections import defaultdict
from os import PathLike

from loguru import logger
from pyteomics.fasta import FASTA
from tqdm.auto import tqdm

FASTA_NAME_REGEX = re.compile(r"^.*\|(.*)\|.*$")
UNIPROT_ACC_REGEX = re.compile(r"^[A-Z0-9]{6}(-\d+)?$")


class ProteinNGram:
    """Implements an n-gram to fast lookup of proteins that match a peptide.

    Usage
    -----
    ```python
    ngram = ProteinNGram.from_fasta("path/to/fasta")
    ngram.search_ngram("AAC")
    ['Prot1']
    ```

    Examples
    --------
    >>> base_ngram = {'AA': {1,3}, 'AB': {2,3}, 'AC': {1, 4}, 'CA': {4}}
    >>> inv_index = {1: "Prot1", 2: "Prot2", 3: "Prot3", 4: "Prot4"}
    >>> inv_seq = {1: "AACAA", 2: "AABAA", 3: "ABDAA", 4: "CACAA"}
    >>> ngram = ProteinNGram(
    ...     ngram = base_ngram,
    ...     inv_alias = inv_index,
    ...     inv_seq = inv_seq)
    >>> ngram.search_ngram("AAC")
    ['Prot1']
    >>> ngram.search_ngram("CAC")
    ['Prot4']
    """

    __slots__ = ("ngram_size", "ngram", "inv_alias", "inv_seq")

    def __init__(
        self,
        ngram: dict[str, set[int]],
        inv_alias: dict[int, str],
        inv_seq: dict[int, str],
    ) -> None:
        """Initialized an ngram for fast lookup.

        For details check the main class docstring.
        """
        keys = list(ngram)
        if not all(len(keys[0]) == len(k) for k in ngram):
            raise ValueError("All ngram keys need to be the same length")
        self.ngram_size: int = len(keys[0])
        self.ngram = ngram
        self.inv_alias = inv_alias
        self.inv_seq = inv_seq

    def search_ngram(self, entry: str) -> list[str]:
        """Searches a sequence using the n-gram and returns the matches."""

        if len(entry) < self.ngram_size:
            raise ValueError(
                f"Entry {entry} is shorter than the n-gram size ({self.ngram_size})"
            )

        candidates = None
        for x in [
            entry[x : x + self.ngram_size]
            for x in range(1 + len(entry) - self.ngram_size)
        ]:
            if len(x) < self.ngram_size:
                raise

            if candidates is None:
                candidates = self.ngram.get(x, set())
            else:
                candidates = candidates.intersection(self.ngram[x])
                if len(candidates) <= 1:
                    break

        # This makes sure the whole sequence is matched.
        # For instance ... "BAAAAB" and "BAAAAAAAAAB" share all the same length 2
        # ngrams, but do not match the same sequence (if the protein is "PEPTIDEBAAAB",
        # only the first would be kept).
        out = [
            self.inv_alias[x] for x in candidates if entry in self.inv_seq[x]
        ]
        return out

    @staticmethod
    def from_fasta(
        fasta_file: PathLike | str,
        ngram_size: int = 4,
        progress: bool = True,
        proteins_keep: None | set[str] = None,
    ) -> ProteinNGram:
        """Builds a protein n-gram from a fasta file.

        Parameters
        ----------
        fasta_file:
            Path-like or string representing the fasta file to read in order
            to build the index.
        ngram_size:
            Size of the chunks that will be used to build the n-gram, should
            be smaller than the smallest peptide to be searched. Longer sequences
            should give a more unique aspect to it but a larger index is built.
        progress:
            Whether to show a progress bar while building the index.
        proteins_keep: set[str]:
            If not None, only keep the proteins in the set, by matching the
            uniprot ID. Example: {'Q92804', 'Q92804-2', 'P13639'}

        """
        ngram = defaultdict(set)
        inv_alias = {}
        inv_seq = {}
        skipped = 0
        kept = 0

        for i, entry in tqdm(
            enumerate(FASTA(str(fasta_file))),
            disable=not progress,
            desc="Building peptide n-gram index",
        ):
            entry_name = FASTA_NAME_REGEX.search(entry.description).group(1)
            if proteins_keep is not None and entry_name not in proteins_keep:
                skipped += 1
                continue
            else:
                kept += 1
            sequence = entry.sequence
            if len(sequence) < ngram_size:
                logger.warning(
                    f"Skipping {entry_name} because it is shorter than the n-gram size"
                )
                continue

            inv_alias[i] = entry_name
            inv_seq[i] = sequence
            for x in [
                sequence[x : x + ngram_size]
                for x in range(1 + len(sequence) - ngram_size)
            ]:
                ngram[x].add(i)

        if proteins_keep is not None:
            logger.info(
                f"Kept {kept} proteins and skipped {skipped} "
                "proteins when importing the fasta file",
            )

        return ProteinNGram(ngram=ngram, inv_alias=inv_alias, inv_seq=inv_seq)


def get_protein_accessions(protein_list: list[str]) -> set[str]:
    """Extracts all protein accessions from a list of protein groups.

    This is meant to be used on the protein groups column in a mokapot
    results file.

    Examples
    --------
        >>> tmp = [
        ...     "sp|Q92804|RBP56_HUMAN",
        ...     "sp|Q92804-2|RBP56_HUMAN",
        ...     "sp|P13639|EF2_HUMAN"]
        >>> got = get_protein_accessions(tmp)
        >>> exp = {'Q92804', 'Q92804-2', 'P13639'}
        >>> exp == got
        True
    """
    out = set()
    for protein in protein_list:
        for accession in protein.split(","):
            acc = accession.split("|")[1]
            out.add(acc)
    return out
