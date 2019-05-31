#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: vec_mod.py
"""
Created on Wed May 15 09:06:57 2019

@author: Neo(liuniu@smail.nju.edu.cn)

Calculate the modulus of vector

"""

import numpy as np

__all__ = ["vec_mod_calc"]


# -----------------------------  FUNCTIONS -----------------------------
def vec_mod_calc(x, xerr=None):
    """
    """

    xmod = np.sqrt(np.sum(x**2))

    if xerr is None:
        return xmod
    else:
        xmoderr = np.sqrt(np.dot(x**2, xerr**2)) / xmod
        return xmod, xmoderr


if __name__ == "__main__":
    pass
# --------------------------------- END --------------------------------
