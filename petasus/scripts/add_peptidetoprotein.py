import click
import sqlite3
import os
import pandas as pd
from petasus.protein_ngram import ProteinNGram, get_protein_accessions


def _add_peptidetoprotein(
    mokapot_proteins: os.PathLike,
    fasta_file: os.PathLike,
    speclib: os.PathLike,
    ngram_size=4,
    q_threshold=0.1,
):
    """Add peptidetoprotein table to dlib/elib file."""
    mokapot_proteins = pd.read_csv(mokapot_proteins, sep="\t")
    mokapot_proteins = mokapot_proteins[
        mokapot_proteins["mokapot q-value"] < q_threshold
    ]

    accessions_keep = get_protein_accessions(
        mokapot_proteins["mokapot protein group"].tolist()
    )

    ngram = ProteinNGram.from_fasta(
        fasta_file,
        proteins_keep=accessions_keep,
        ngram_size=ngram_size,
    )

    conn = sqlite3.connect(speclib)
    df = pd.read_sql("SELECT PeptideSeq FROM entries", conn)
    conn.close()
    peps = df["PeptideSeq"].unique().tolist()
    del df

    matches = [ngram.search_ngram(p) for p in peps]
    num_marches = sum([len(m) for m in matches])

    peptide_col = [None] * num_marches
    protein_col = [None] * num_marches

    i = 0
    for p, m in zip(peps, matches):
        for n in m:
            peptide_col[i] = p
            protein_col[i] = n
            i += 1

    df = pd.DataFrame(
        {"PeptideSeq": peptide_col, "ProteinAccession": protein_col}
    )
    df["isDecoy"] = False
    conn = sqlite3.connect(speclib)
    conn.execute("DROP TABLE IF EXISTS peptidetoprotein")
    conn.execute(
        (
            "CREATE TABLE peptidetoprotein "
            "(PeptideSeq string not null, "
            "isDecoy boolean, "
            "ProteinAccession string not null)"
        )
    )
    df.to_sql("peptidetoprotein", conn, if_exists="append", index=False)
    conn.close()


@click.command()
@click.argument("mokapot_proteins", type=click.Path(exists=True))
@click.argument("fasta_file", type=click.Path(exists=True))
@click.argument("speclib", type=click.Path(exists=True))
@click.option(
    "--ngram_size",
    type=int,
    default=4,
    help=(
        "Size of the chunks that will be used to build the n-gram,"
        " should be smaller than the smallest peptide to be searched."
        " Longer sequences should give a more unique aspect to it"
        " but a larger index is built.",
    ),
)
@click.option(
    "--q_threshold",
    type=float,
    default=0.1,
    help="Threshold for mokapot q-value.",
)
def add_peptidetoprotein(
    mokapot_proteins,
    fasta_file,
    speclib,
    ngram_size=4,
    q_threshold=0.1,
):
    """Add peptidetoprotein table to dlib/elib file."""
    _add_peptidetoprotein(
        mokapot_proteins=mokapot_proteins,
        fasta_file=fasta_file,
        speclib=speclib,
        ngram_size=ngram_size,
        q_threshold=q_threshold,
    )


if __name__ == "__main__":
    add_peptidetoprotein()
