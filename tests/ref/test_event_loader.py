from ref_builder.events.loader import event_adapter


def test_load_event(scratch_repo):
    """Test loading all events from the scratch repo."""
    for path in (scratch_repo.path / "src").glob("*.json"):
        with path.open() as f:
            event = event_adapter.validate_json(f.read())

        assert type(event).__name__ == str(event.type)
