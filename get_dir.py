#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: get_dir.py
"""
Created on Fri Sep 11 11:06:14 2020

@author: Neo(liuniu@smail.nju.edu.cn)

Get the data directory
"""

import sys
from sys import platform as _platform


# -----------------------------  MAIN -----------------------------
def get_data_dir():
    """Check the OS and get the data directory.

    I have several computers with different OS.
    The path to the data file is different.

    Returns
    -------
    datadir: string
        directory which stores Gaia catalogs.

    """

    # Check the type of OS
    if _platform in ["linux", "linux2"]:
        # linux
        data_dir = "/home/neo/Astronomy/data/catalogs"
    elif _platform in ["darwin"]:
        # MAC OS X
        data_dir = "/Users/Neo/Astronomy/data/catalogs"
    elif _platform in ["win32", "win64"]:
        # Windows
        print("Not implemented yet")
        sys.exit()

    return data_dir

# --------------------------------- END --------------------------------
