from ref_builder.repo import Repo
from tests.fixtures.factories import OTUFactory


def test_commit(empty_repo: Repo, otu_factory: OTUFactory):
    otu = otu_factory.build()

    with empty_repo.use_transaction() as transaction:
        empty_repo.create_otu(
            otu.acronym,
            otu.legacy_id,
            molecule=otu.molecule,
            name=otu.name,
            plan=otu.plan,
            taxid=otu.taxid,
        )

        empty_repo.create_isolate(
            otu.id,
            "isolate_id",
            name=otu.isolates[0].name,
        )

        # Make sure the event hasn't been applied to the repo.
        assert empty_repo.last_id == 1

        assert len(list(transaction.iter_events())) == 1

    assert empty_repo.last_id == 2
    assert len(list(empty_repo.iter_otus())) == 0


def test_abort(empty_repo: Repo, otu_factory: OTUFactory):
    otu = otu_factory.build()

    with empty_repo.use_transaction() as transaction:
        empty_repo.create_otu(
            otu.acronym,
            otu.legacy_id,
            molecule=otu.molecule,
            name=otu.name,
            plan=otu.plan,
            taxid=otu.taxid,
        )

        assert len(list(transaction.iter_events())) == 1

        transaction.abort()

        # This shouldn't be raised because the transaction was aborted.
        raise ValueError("Transaction was not aborted.")

    assert empty_repo.last_id == 1
    assert len(list(empty_repo.iter_otus())) == 0
