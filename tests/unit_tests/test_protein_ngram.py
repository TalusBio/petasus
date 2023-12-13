from petasus.protein_ngram import ProteinNGram, get_protein_accessions
import pytest


def test_ngram(shared_datadir):
    fasta_file = shared_datadir / "Q99536.fasta"

    ngram = ProteinNGram.from_fasta(str(fasta_file), 4)
    assert ngram.ngram_size == 4
    assert len(ngram.search_ngram("PEPTIDE")) == 0
    assert len(ngram.search_ngram("LNRSGMWQEEVTVP")) == 1
    assert len(ngram.search_ngram("QEEVTVPLNRSGMW")) == 0

    res = ngram.search_ngram("LNRSGMWQEEVTVP")
    assert res[0] == "Q99536"

    res = ngram.search_ngram("VEALYLVCGERG")
    assert res[0] == "P01308"

    with pytest.raises(ValueError):
        res = ngram.search_ngram("VEA")


def test_id_extractor():
    tmp = [
        "sp|Q92804|RBP56_HUMAN",
        "sp|Q92804-2|RBP56_HUMAN",
        "sp|P13639|EF2_HUMAN",
    ]
    got = get_protein_accessions(tmp)
    exp = {"Q92804", "Q92804-2", "P13639"}
    assert exp == got
