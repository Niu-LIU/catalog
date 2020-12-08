#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: read_gaia.py
"""
Created on Thu Oct  4 16:00:02 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""

from astropy.table import Table, Column
from astropy import units as u
# from astropy.coordinates import SkyCoord
import numpy as np

# My modules
from .pos_err import error_ellipse_array
from .get_dir import get_data_dir

__all__ = ["read_dr1_qso", "read_dr2_iers", "read_dr2_allwise",
           "read_edr3_crf", "read_edr3_agn"]


# -----------------------------  FUNCTIONS -----------------------------
def read_dr1_qso(dr1_qso_file=None):
    """Read the position information of Gaia DR1 quasar auxiliary solution.

    Parameter
    ---------
    dr1_qso_file : string
        file name and path of the Gaia DR1 quasar auxiliary solution

    Return
    ------
    gdr1 : an astropy.Table object
        data in the catalog
    """

    if dr1_qso_file is None:
        data_dir = get_data_dir()
        dr1_qso_file = "{}/gaia/dr1/qso.dat".format(data_dir)

    gdr1 = Table.read(dr1_qso_file, format="ascii.fixed_width_no_header",
                      names=["solution_id", "source_id", "ref_epoch",
                             "ra", "ra_err", "dec", "dec_err",
                             "ra_dec_corr", "phot_g_mean_mag",
                             "astrometric_priors_used", "icrf2_match",
                             "rot_flag"])

    # Add unit information
    gdr1["ra_err"].unit = u.mas
    gdr1["dec_err"].unit = u.mas

    # Calculate the semi-major axis of error ellipse
    pos_err, pos_err_min, pa = error_ellipse_array(
        gdr1["ra_err"], gdr1["dec_err"], gdr1["ra_dec_corr"])
    del pos_err_min

    # Add the semi-major axis of error ellipse to the table
    pos_err = Column(pos_err, name="pos_err", unit=u.mas)
    pa = Column(pa, name="eepa", unit=u.deg)
    gdr1.add_columns([pos_err, pa], indexes=[9, 9])

    return gdr1


def modify_dr2_iers(dr2_qso_file=None):
    """Correct the wrong name in the Gaia-CRF2 subset.

    Parameter
    ---------
    dr2_qso_file : string
        file name and path of the Gaia DR2 auxiliary IERS catalog

    Return
    ------
    gdr2 : an astropy.Table object
        data in the catalog
    """

    if dr2_qso_file is None:
        data_dir = get_data_dir()
        dr2_qso_file = "{}/gaia/dr2/gaiadr2_iers0.fits".format(data_dir)

    # Read Gaia DR2 IERS quasar data
    gdr2 = Table.read(dr2_qso_file)

    # There are two small errore in the colomn iers_name in this sample
    # 0548+37A --> 0548+377
    # 1954+188 --> 1954+187

    oldnames = ["0548+37A", "1954+188"]
    newnames = ["0548+377", "1954+187"]

    for oldname, newname in zip(oldnames, newnames):
        idx = np.where(gdr2["iers_name"] == oldname)[0][0]
        gdr2[idx]["iers_name"] = newname

    gdr2.write(dr2_qso_file, overwrite=True)


def read_dr2_iers(dr2_qso_file=None, errscaling=False):
    """Read the positional information of Gaia DR2 auxiliary IERS catalog.

    Parameter
    ---------
    dr2_qso_file : string
        file name and path of the Gaia DR2 auxiliary IERS catalog
    errscaling : boolean
        if true, scale gaia position error by chi2/ndf

    Return
    ------
    gdr2 : an astropy.Table object
        data in the catalog
    """

    if dr2_qso_file is None:
        data_dir = get_data_dir()
        dr2_qso_file = "{}/gaia/dr2/gaiadr2_iers.fits".format(data_dir)

    # Read Gaia DR2 IERS quasar data
    gdr2 = Table.read(dr2_qso_file)

    # There are two small errore in the colomn iers_name in this sample
    # 0548+37A --> 0548+377
    # 1954+188 --> 1954+187

    # oldnames = ["0548+37A", "1954+188"]
    # newnames = ["0548+377", "1954+187"]
    #
    # for oldname, newname in zip(oldnames, newnames):
    #     idx = np.where(gdr2["iers_name"] == oldname)[0][0]
    #     gdr2[idx]["iers_name"] = newname

    # Only the positional information are kept.
    gdr2.keep_columns(["iers_name", "source_id",
                       "ra", "ra_error", "dec", "dec_error",
                       "parallax", "parallax_error",
                       "pmra", "pmra_error", "pmdec", "pmdec_error",
                       "ra_dec_corr", "ra_parallax_corr", "ra_pmra_corr",
                       "ra_pmdec_corr", "dec_parallax_corr", "dec_pmra_corr",
                       "dec_pmdec_corr", "parallax_pmra_corr",
                       "parallax_pmdec_corr", "pmra_pmdec_corr",
                       "phot_g_mean_mag", "phot_bp_mean_mag", "phot_rp_mean_mag",
                       "bp_rp", "astrometric_n_obs_al", "astrometric_n_good_obs_al",
                       "astrometric_matched_observations", "astrometric_chi2_al",
                       "astrometric_params_solved"])

    # Rename the column names
    gdr2.rename_column("ra_error", "ra_err")
    gdr2.rename_column("dec_error", "dec_err")
    gdr2.rename_column("parallax_error", "parallax_err")
    gdr2.rename_column("pmra_error", "pmra_err")
    gdr2.rename_column("pmdec_error", "pmdec_err")

    # Determine whether to scale the position error
    if errscaling:
        # Degree of freedom = #observation - #parameter
        dof = gdr2["astrometric_n_good_obs_al"] - \
            gdr2["astrometric_params_solved"]
        scale = np.sqrt(gdr2["astrometric_chi2_al"] / dof)

        gdr2["ra_err"] = gdr2["ra_err"] * scale
        gdr2["dec_err"] = gdr2["dec_err"] * scale

    # Calculate the semi-major axis of error ellipse
    pos_err, pos_err_min, pa = error_ellipse_array(
        gdr2["ra_err"], gdr2["dec_err"], gdr2["ra_dec_corr"])
    del pos_err_min

    # Add the semi-major axis of error ellipse to the table
    pos_err = Column(pos_err, name="pos_err", unit=u.mas)
    pa = Column(pa, name="eepa", unit=u.deg)
    gdr2.add_columns([pos_err, pa], indexes=[9, 9])

    return gdr2


def read_dr2_allwise(dr2_qso_file=None):
    """Read the positional information of Gaia DR2 auxiliary AllWISE catalog.

    Parameter
    ---------
    dr2_qso_file : string
        file name and path of the Gaia DR2 auxiliary AllWISE catalog

    Return
    ------
    gdr2 : an astropy.Table object
        data in the catalog
    """

    if dr2_qso_file is None:
        data_dir = get_data_dir()
        dr2_qso_file = "{}/gaia/dr2/gaiadr2_qso_all.fits".format(data_dir)

    gdr2 = read_dr2_iers(dr2_qso_file)

    # Read Gaia DR2 IERS quasar data
    # gdr2 = Table.read(dr2_file)

    # # Only the positional information are kept.
    # gdr2.keep_columns(["iers_name",
    #                    "source_id",
    #                    "ra",
    #                    "ra_error",
    #                    "dec",
    #                    "dec_error",
    #                    "parallax",
    #                    "parallax_error",
    #                    "pmra",
    #                    "pmra_error",
    #                    "pmdec",
    #                    "pmdec_error",
    #                    "ra_dec_corr",
    #                    "ra_parallax_corr",
    #                    "ra_pmra_corr",
    #                    "ra_pmdec_corr",
    #                    "dec_parallax_corr",
    #                    "dec_pmra_corr",
    #                    "dec_pmdec_corr",
    #                    "parallax_pmra_corr",
    #                    "parallax_pmdec_corr",
    #                    "pmra_pmdec_corr",
    #                    "phot_g_mean_mag",
    #                    "phot_bp_mean_mag",
    #                    "phot_rp_mean_mag"])

    # # Rename the column names
    # gdr2.rename_column("ra_error", "ra_err")
    # gdr2.rename_column("dec_error", "dec_err")
    # gdr2.rename_column("parallax_error", "parallax_err")
    # gdr2.rename_column("pmra_error", "pmra_err")
    # gdr2.rename_column("pmdec_error", "pmdec_err")
    #
    # # Calculate the semi-major axis of error ellipse
    # pos_err = pos_err_calc(
    #     gdr2["ra_err"], gdr2["dec_err"], gdr2["ra_dec_corr"])
    #
    # # Add the semi-major axis of error ellipse to the table
    # gdr2.add_column(pos_err, name="pos_err", index=6)
    # gdr2["pos_err"].unit = u.mas

    return gdr2


def read_edr3_agn(data_file=None, astro_only=False):
    """Read the positional information of Gaia EDR3 AGN catalogs

    Parameter
    ---------
    data_file : string
        file name and path of the Gaia EDR3 AGN catalogs

    Return
    ------
    gedr3_tab : an astropy.Table object
        data in the catalog
    """

    if data_file is None:
        data_dir = get_data_dir()
        data_file = "{}/gaia/edr3/gedr3_agn.fits".format(data_dir)

    # Read data
    gedr3_tab = Table.read(data_file)

    # # Only the positional information are kept.
    # gedr3_tab.keep_columns(["iers_name",
    #                    "source_id",
    #                    "ra",
    #                    "ra_error",
    #                    "dec",
    #                    "dec_error",
    #                    "parallax",
    #                    "parallax_error",
    #                    "pmra",
    #                    "pmra_error",
    #                    "pmdec",
    #                    "pmdec_error",
    #                    "ra_dec_corr",
    #                    "ra_parallax_corr",
    #                    "ra_pmra_corr",
    #                    "ra_pmdec_corr",
    #                    "dec_parallax_corr",
    #                    "dec_pmra_corr",
    #                    "dec_pmdec_corr",
    #                    "parallax_pmra_corr",
    #                    "parallax_pmdec_corr",
    #                    "pmra_pmdec_corr",
    #                    "phot_g_mean_mag",
    #                    "phot_bp_mean_mag",
    #                    "phot_rp_mean_mag"])

    # Rename the column names
    gedr3_tab.rename_column("ra_error", "ra_err")
    gedr3_tab.rename_column("dec_error", "dec_err")
    gedr3_tab.rename_column("parallax_error", "parallax_err")
    gedr3_tab.rename_column("pmra_error", "pmra_err")
    gedr3_tab.rename_column("pmdec_error", "pmdec_err")

#     # Calculate the semi-major axis of error ellipse
#     pos_err = pos_err_calc(
#         gedr3_tab["ra_err"], gedr3_tab["dec_err"], gedr3_tab["ra_dec_corr"])
#
#     # Add the semi-major axis of error ellipse to the table
#     gedr3_tab.add_column(pos_err, name="pos_err", index=6)
#     gedr3_tab["pos_err"].unit = u.mas

    return gedr3_tab


def read_edr3_crf(data_file=None, astro_only=False, type="All"):
    """Read the positional information of Gaia EDR3 CRF source catalog

    Parameter
    ---------
    data_file : string
        file name and path of the Gaia-CRF3 catalogs

    Return
    ------
    gedr3_tab : an astropy.Table object
        data in the catalog
    """

    if data_file is None:
        data_dir = get_data_dir()
        data_file = "{}/gaia/edr3/gedr3_frame_rotator_source.fits".format(
            data_dir)

    # Read data
    gedr3_tab = Table.read(data_file)

    # # Only the positional information are kept.
    # gedr3_tab.keep_columns(["iers_name",
    #                    "source_id",
    #                    "ra",
    #                    "ra_error",
    #                    "dec",
    #                    "dec_error",
    #                    "parallax",
    #                    "parallax_error",
    #                    "pmra",
    #                    "pmra_error",
    #                    "pmdec",
    #                    "pmdec_error",
    #                    "ra_dec_corr",
    #                    "ra_parallax_corr",
    #                    "ra_pmra_corr",
    #                    "ra_pmdec_corr",
    #                    "dec_parallax_corr",
    #                    "dec_pmra_corr",
    #                    "dec_pmdec_corr",
    #                    "parallax_pmra_corr",
    #                    "parallax_pmdec_corr",
    #                    "pmra_pmdec_corr",
    #                    "phot_g_mean_mag",
    #                    "phot_bp_mean_mag",
    #                    "phot_rp_mean_mag"])

    # Rename the column names
    gedr3_tab.rename_column("ra_error", "ra_err")
    gedr3_tab.rename_column("dec_error", "dec_err")
    gedr3_tab.rename_column("parallax_error", "parallax_err")
    gedr3_tab.rename_column("pmra_error", "pmra_err")
    gedr3_tab.rename_column("pmdec_error", "pmdec_err")

    # Calculate the semi-major axis of error ellipse
#     pos_err = pos_err_calc(
#         gedr3_tab["ra_err"], gedr3_tab["dec_err"], gedr3_tab["ra_dec_corr"])

    # Add the semi-major axis of error ellipse to the table
#     gedr3_tab.add_column(pos_err, name="pos_err", index=6)
#     gedr3_tab["pos_err"].unit = u.mas

    return gedr3_tab


if __name__ == "__main__":
    modify_dr2_iers()
# --------------------------------- END --------------------------------
