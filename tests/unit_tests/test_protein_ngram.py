from petasus.protein_ngram import ProteinNGram


def test_ngram(shared_datadir):
    fasta_file = shared_datadir / "Q99536.fasta"

    ngram = ProteinNGram.from_fasta(str(fasta_file), 4)
    assert ngram.ngram_size == 4
    assert len(ngram.search_ngram("PEPTIDE")) == 0
    assert len(ngram.search_ngram("LNRSGMWQEEVTVP")) == 1
    assert len(ngram.search_ngram("QEEVTVPLNRSGMW")) == 0
