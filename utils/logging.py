"""
utils/logging.py – Structured stage-aware logger.
"""

import logging
import sys
import time
from contextlib import contextmanager
from typing import Generator


_FMT = "[%(levelname)s %(asctime)s] %(message)s"
_DATE_FMT = "%H:%M:%S"


def get_logger(name: str = "ncrp") -> logging.Logger:
    logger = logging.getLogger(name)
    if not logger.handlers:
        handler = logging.StreamHandler(sys.stderr)
        handler.setFormatter(logging.Formatter(_FMT, _DATE_FMT))
        logger.addHandler(handler)
    logger.setLevel(logging.INFO)
    return logger


@contextmanager
def stage_timer(logger: logging.Logger, stage_name: str) -> Generator[None, None, None]:
    """Context manager that logs stage start/end with elapsed time."""
    logger.info(f"{'='*10} {stage_name} {'='*10}")
    t0 = time.perf_counter()
    try:
        yield
    finally:
        elapsed = time.perf_counter() - t0
        logger.info(f"[{stage_name}] Done in {elapsed:.2f}s")
