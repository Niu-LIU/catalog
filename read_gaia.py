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
from my_progs.catalog.pos_err import pos_err_calc


__all__ = ["read_dr1_qso", "read_dr2_qso", "read_dr2_allwise"]


# -----------------------------  FUNCTIONS -----------------------------
def read_dr1_qso(dr1_file="/Users/Neo/Astronomy/data/catalogs/"
                 "gaia/dr1/qso.dat"):
    """Read the position information of Gaia DR1 quasar auxiliary solution.
    """

    gdr1 = Table.read(dr1_file, format="ascii.fixed_width_no_header",
                      names=["solution_id", "source_id", "ref_epoch",
                             "ra", "ra_err", "dec", "dec_err",
                             "ra_dec_corr", "phot_g_mean_mag",
                             "astrometric_priors_used", "icrf2_match",
                             "rot_flag"])

    # Add unit information
    gdr1["ra_err"].unit = u.mas
    gdr1["dec_err"].unit = u.mas

    # Calculate the semi-major axis of error ellipse
    pos_err = pos_err_calc(
        gdr1["ra_err"], gdr1["dec_err"], gdr1["ra_dec_corr"])

    # Add the semi-major axis of error ellipse to the table
    gdr1.add_column(pos_err, name="pos_err", index=9)
    gdr1["pos_err"].unit = u.mas

    return gdr1


def read_dr2_qso(dr2_file="/Users/Neo/Astronomy/data/catalogs/"
                 "gaia/dr2/gaiadr2_iers.fits"):
    """Read the positional information of Gaia DR2 auxiliary IERS catalog.
    """

    # Read Gaia DR2 IERS quasar data
    gdr2 = Table.read(dr2_file)

    # Only the positional information are kept.
    gdr2.keep_columns(["iers_name",
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

    # Rename the column names
    gdr2.rename_column("ra_error", "ra_err")
    gdr2.rename_column("dec_error", "dec_err")
    gdr2.rename_column("parallax_error", "parallax_err")
    gdr2.rename_column("pmra_error", "pmra_err")
    gdr2.rename_column("pmdec_error", "pmdec_err")

    # Calculate the semi-major axis of error ellipse
    pos_err = pos_err_calc(
        gdr2["ra_err"], gdr2["dec_err"], gdr2["ra_dec_corr"])

    # Add the semi-major axis of error ellipse to the table
    gdr2.add_column(pos_err, name="pos_err", index=6)
    gdr2["pos_err"].unit = u.mas

    return gdr2


def read_dr2_allwise(dr2_file="/Users/Neo/Astronomy/Data/catalogs/"
                     "Gaia_DR2/gaiadr2_iers.fits"):
    """Read the positional information of Gaia DR2 auxiliary IERS catalog.

    """

    # Read Gaia DR2 IERS quasar data
    gdr2 = Table.read(dr2_file)

    # Only the positional information are kept.
    gdr2.keep_columns(["iers_name",
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

    # Rename the column names
    gdr2.rename_column("ra_error", "ra_err")
    gdr2.rename_column("dec_error", "dec_err")
    gdr2.rename_column("parallax_error", "parallax_err")
    gdr2.rename_column("pmra_error", "pmra_err")
    gdr2.rename_column("pmdec_error", "pmdec_err")

    # Calculate the semi-major axis of error ellipse
    pos_err = pos_err_calc(
        gdr2["ra_err"], gdr2["dec_err"], gdr2["ra_dec_corr"])

    # Add the semi-major axis of error ellipse to the table
    gdr2.add_column(pos_err, name="pos_err", index=6)
    gdr2["pos_err"].unit = u.mas

    return gdr2
# --------------------------------- END --------------------------------
