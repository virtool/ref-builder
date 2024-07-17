import logging
import os
import sys

import click
import structlog

error_message_style = click.style("ERROR: ", fg="red")

logging.basicConfig(
    format="%(message)s",
    stream=sys.stderr,
    level=logging.DEBUG,
)


def configure_logger(debug: bool = False) -> None:
    """Configure structlog-based logging.

    :param debug: Debugging flag. Defaults to False.
    """
    processors = [
        structlog.stdlib.filter_by_level,
        structlog.stdlib.add_logger_name,
        structlog.stdlib.add_log_level,
        structlog.stdlib.PositionalArgumentsFormatter(),
        structlog.processors.StackInfoRenderer(),
        structlog.processors.UnicodeDecoder(),
    ]

    if debug:
        processors.append(
            structlog.processors.CallsiteParameterAdder(
                {
                    structlog.processors.CallsiteParameter.FILENAME,
                    structlog.processors.CallsiteParameter.FUNC_NAME,
                    structlog.processors.CallsiteParameter.LINENO,
                },
            ),
        )

    if sys.stderr.isatty() and not os.environ.get("NO_COLOR"):
        processors.append(structlog.dev.ConsoleRenderer())
    else:
        processors += [
            structlog.processors.TimeStamper(fmt="iso"),
            structlog.processors.format_exc_info,
            structlog.processors.dict_tracebacks,
            structlog.processors.JSONRenderer(),
        ]

    structlog.configure(
        processors=processors,
        logger_factory=structlog.stdlib.LoggerFactory(),
        cache_logger_on_first_use=True,
    )

    logger_level = logging.DEBUG if debug else logging.INFO

    structlog.configure(
        wrapper_class=structlog.make_filtering_bound_logger(logger_level),
    )
