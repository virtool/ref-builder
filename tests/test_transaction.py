from ref_builder.repo import Repo
from tests.fixtures.factories import OTUFactory


def test_commit(empty_repo: Repo, otu_factory: OTUFactory):
    fake_otu = otu_factory.build()

    with empty_repo.lock(), empty_repo.use_transaction():
        otu = empty_repo.create_otu(
            fake_otu.acronym,
            fake_otu.legacy_id,
            molecule=fake_otu.molecule,
            name=fake_otu.name,
            plan=fake_otu.plan,
            taxid=fake_otu.taxid,
        )

        empty_repo.create_isolate(
            otu.id,
            "isolate_id",
            name=fake_otu.isolates[0].name,
        )

    assert empty_repo.last_id == 3
    assert len(list(empty_repo.iter_otus())) == 1


def test_abort(empty_repo: Repo, otu_factory: OTUFactory):
    otu = otu_factory.build()

    with empty_repo.lock(), empty_repo.use_transaction() as transaction:
        empty_repo.create_otu(
            otu.acronym,
            otu.legacy_id,
            molecule=otu.molecule,
            name=otu.name,
            plan=otu.plan,
            taxid=otu.taxid,
        )

        assert empty_repo.last_id == 2
        assert len(list(empty_repo.iter_otus())) == 1

        transaction.abort()

        # This shouldn't be raised because the transaction was aborted.
        raise ValueError("Transaction was not aborted.")

    assert empty_repo.last_id == 1
    assert len(list(empty_repo.iter_otus())) == 0
