import pytest

from ref_builder.legacy.utils import ErrorHandledResult
from ref_builder.legacy.validate import OTUValidationResult, log_otu_validation_result


def test_log_ok(capfd, mocker):
    """Test that the console logs OK when validation passes."""
    result = OTUValidationResult(
        handler_results=[],
        repaired_otu={"name": "Test virus"},
    )

    log_otu_validation_result("Test virus", result, False)

    assert capfd.readouterr().out == (
        "Test virus "
        "─────────────────────────────────────────────────────────────────────\n\n"
        "  ‣ OK                                                                   "
        "       \n\n"
    )


def test_log_no_ok(capfd):
    """Test that the console logs nothing when the validation passes and ``no_ok`` is
    set.
    """
    result = OTUValidationResult(
        handler_results=[],
        repaired_otu={"name": "Test virus"},
    )

    log_otu_validation_result("Test virus", result, True)

    assert capfd.readouterr().out == ""


@pytest.mark.parametrize("no_ok", [False, True])
def test_log_error_single(no_ok: bool, capfd):
    """Test that the console logs OK correctly."""
    result = OTUValidationResult(
        handler_results=[
            ErrorHandledResult(
                message="Error message",
                fixed=False,
            ),
        ],
        repaired_otu={"name": "Test virus"},
    )

    log_otu_validation_result("Test virus", result, no_ok)

    assert capfd.readouterr().out == (
        "Test virus "
        "─────────────────────────────────────────────────────────────────────\n\n"
        "  ‣ [ERROR] Error message                                                "
        "       \n\n"
    )


@pytest.mark.parametrize("no_ok", [False, True])
def test_log_error_multi(no_ok: bool, capfd):
    """Test that the console logs OK correctly."""
    result = OTUValidationResult(
        handler_results=[
            ErrorHandledResult(
                message="There was a problem",
                fixed=False,
            ),
            ErrorHandledResult(
                message="There was an even bigger problem",
                fixed=False,
            ),
        ],
        repaired_otu={"name": "Test virus"},
    )

    log_otu_validation_result("Test virus", result, no_ok)

    assert capfd.readouterr().out == (
        "Test virus "
        "─────────────────────────────────────────────────────────────────────\n\n"
        "  ‣ [ERROR] There was a problem                                          "
        "       \n"
        "  ‣ [ERROR] There was an even bigger problem                             "
        "       \n\n"
    )


@pytest.mark.parametrize("no_ok", [False, True])
def test_log_error_fixed(no_ok: bool, capfd):
    """Test that the console logs OK correctly."""
    result = OTUValidationResult(
        handler_results=[
            ErrorHandledResult(
                message="A fixable problem occurred",
                fixed=True,
            ),
        ],
        repaired_otu={"name": "Test virus"},
    )

    log_otu_validation_result("Test virus", result, no_ok)

    assert capfd.readouterr().out == (
        "Test virus "
        "─────────────────────────────────────────────────────────────────────\n\n"
        "  ‣ [FIXED] A fixable problem occurred                                   "
        "       \n\n"
    )
