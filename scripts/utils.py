#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Backward-compatibility shim.

The canonical implementations now live in the scripts/lib/ package
(SCALE_MLOPS_PLAN.md §5):
    get_y_chunk  -> lib.chunking
    run_command  -> lib.io_utils

This module is kept so existing ``from utils import get_y_chunk, run_command``
imports continue to work. New code should import from lib directly.
"""

from lib.chunking import get_y_chunk
from lib.io_utils import run_command

__all__ = ["get_y_chunk", "run_command"]
