#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: read_gaia.py
"""
Created on Thu Oct  4 16:00:02 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""

from astropy.table import Table, Column
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np
import sys
# My modules
from .pos_err import pos_err_calc


__all__ = ["read_gdr2_qso"]


# -----------------------------  FUNCTIONS -----------------------------
def read_gdr2_qso(gaiadr2_file="/Users/Neo/Astronomy/Data/catalogs/"
                     "Gaia_DR2/gaiadr2_iers.fits"):
    """Read the positional information of Gaia DR2 auxiliary IERS catalog.

    """

    # Read Gaia DR2 IERS quasar data
    gaiadr2 = Table.read(gaiadr2_file)

    # Only the positional information are kept.
    gaiadr2.keep_columns(["iers_name",
                          "source_id",
                          "ra",
                          "ra_error",
                          "dec",
                          "dec_error",
                          "parallax",
                          "parallax_error",
                          "pmra",
                          "pmra_error",
                          "pmdec",
                          "pmdec_error",
                          "ra_dec_corr",
                          "ra_parallax_corr",
                          "ra_pmra_corr",
                          "ra_pmdec_corr",
                          "dec_parallax_corr",
                          "dec_pmra_corr",
                          "dec_pmdec_corr",
                          "parallax_pmra_corr",
                          "parallax_pmdec_corr",
                          "pmra_pmdec_corr",
                          "phot_g_mean_mag",
                          "phot_bp_mean_mag",
                          "phot_rp_mean_mag"])

    # Calculate the semi-major axis of error ellipse
    pos_err = pos_err_calc(
        gaiadr2["ra_error"], gaiadr2["dec_error"], gaiadr2["ra_dec_corr"])

    # Add the semi-major axis of error ellipse to the table
    gaiadr2.add_column(pos_err, name="pos_error", index=6)
    gaiadr2["pos_error"].unit = u.mas

    return gaiadr2

def read_gdr2_allwise(gaiadr2_file="/Users/Neo/Astronomy/Data/catalogs/"
                     "Gaia_DR2/gaiadr2_iers.fits"):
    """Read the positional information of Gaia DR2 auxiliary IERS catalog.

    """

    # Read Gaia DR2 IERS quasar data
    gaiadr2 = Table.read(gaiadr2_file)

    # Only the positional information are kept.
    gaiadr2.keep_columns(["iers_name",
                          "source_id",
                          "ra",
                          "ra_error",
                          "dec",
                          "dec_error",
                          "parallax",
                          "parallax_error",
                          "pmra",
                          "pmra_error",
                          "pmdec",
                          "pmdec_error",
                          "ra_dec_corr",
                          "ra_parallax_corr",
                          "ra_pmra_corr",
                          "ra_pmdec_corr",
                          "dec_parallax_corr",
                          "dec_pmra_corr",
                          "dec_pmdec_corr",
                          "parallax_pmra_corr",
                          "parallax_pmdec_corr",
                          "pmra_pmdec_corr",
                          "phot_g_mean_mag",
                          "phot_bp_mean_mag",
                          "phot_rp_mean_mag"])

    # Calculate the semi-major axis of error ellipse
    pos_err = pos_err_calc(
        gaiadr2["ra_error"], gaiadr2["dec_error"], gaiadr2["ra_dec_corr"])

    # Add the semi-major axis of error ellipse to the table
    gaiadr2.add_column(pos_err, name="pos_error", index=6)
    gaiadr2["pos_error"].unit = u.mas

    return gaiadr2
# --------------------------------- END --------------------------------
