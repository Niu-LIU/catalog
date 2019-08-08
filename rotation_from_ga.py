#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: rotation_from_ga.py
"""
Created on Sat Jul 20 11:14:05 2019

@author: Neo(liuniu@smail.nju.edu.cn)

The rotation of the celestial frame induced by the Galactic aberration (GA).


"""

import numpy as np
# My modules
from .glide_calc import glide_gen, glide_field_gen, glide_calc
from .vsh_deg1_cor import vsh_deg01_fitting


# -----------------------------  FUNCTIONS -----------------------------
def rotation_from_ga(cat):
    """Rotation of the celestial frame induced by the Galactic aberration (GA).

    Parameter
    --------
    cat: astropy.table.Table object
        contain at least "ra" and "dec" column

    Returns
    --------
    rot: np.array of (3,)
        rotational angles around three axis
    sig: np.array of (3,)
        formal error of the rotation
    cov: np.array of (3,3)
        covariance matrix
    """
    radeg = np.array(cat["ra"])
    decdeg = np.array(cat["dec"])

    rarad = np.deg2rad(radeg)
    decrad = np.deg2rad(decdeg)

    # GA constant (Liu et al. 2012)
    raGC = 266.4
    decGC = -28.4
    # ICRF3 adopted value=
    gmod = 5.8

    # Generate a dipolar field
    gvec = glide_gen(gmod, raGC, decGC)
    pmra, pmdec = glide_field_gen(radeg, decdeg, gvec)

    # Fit the rotation
    rot, sig, cov, _, _, _ = vsh_deg01_fitting(pmra, pmdec, rarad, decrad,
                                               elim_flag=None,
                                               fit_type="rotation")

    return rot, sig, cov


if __name__ == '__main__':
    rotation_from_ga()

# --------------------------------- END --------------------------------
