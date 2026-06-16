#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Standard logging helper for the pipeline.

The numbered scripts historically use ``print`` for their human-facing progress
output; this module provides a consistent, timestamped logger for the
orchestrator (``run_pipeline.py``) and any new code, optionally tee-ing to a
file under ``logs/``. Keeping it here makes the logging format a single source
of truth.
"""

from __future__ import annotations

import logging
from pathlib import Path


def get_logger(name: str = "amr", logfile: Path | str | None = None,
               level: int = logging.INFO) -> logging.Logger:
    """
    Return a configured logger that writes to stderr and, optionally, a file.

    Args:
        name:    logger name (use one per component, e.g. "pipeline").
        logfile: if given, also append to this path (parents are created).
        level:   logging level (default INFO).
    """
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.handlers.clear()
    logger.propagate = False

    fmt = logging.Formatter("%(asctime)s %(levelname)-7s %(name)s | %(message)s",
                            datefmt="%Y-%m-%d %H:%M:%S")

    stream = logging.StreamHandler()
    stream.setFormatter(fmt)
    logger.addHandler(stream)

    if logfile is not None:
        logfile = Path(logfile)
        logfile.parent.mkdir(parents=True, exist_ok=True)
        fileh = logging.FileHandler(logfile, encoding="utf-8")
        fileh.setFormatter(fmt)
        logger.addHandler(fileh)

    return logger
