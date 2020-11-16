#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: read_ocars.py
"""
Created on Fri Sep 11 11:05:17 2020

@author: Neo(liuniu@smail.nju.edu.cn)

# 103 - 104  object primary type:
#              A  - AGN, including unspecified QSO in NED
#              AB - blazar (replaced by other AGN class when known)
#              AL - BL Lac type
#              AQ - quasar
#              AR - LINER-type AGN
#              AS - Seyfert galaxy
#              A1 - Seyfert 1 galaxy
#              A2 - Seyfert 2 galaxy
#              G  - radio galaxy
#              G2 - double or multiple galaxy
#              I  - IR source
#              R  - radio source
#              S  - star or stellar object
#              SN - supernova (remnant)
#              U  - UV source
#              V  - visual source
#              X  - X-ray source
#              GR - gamma-ray source

"""

import numpy as np
from astropy.table import Table
from astropy import units as u

# My modules
from .get_dir import get_data_dir


# -----------------------------  MAIN -----------------------------
def read_ocars(ocars_file=None):
    """Read the OCARS catalog

    Information about the OCARS catalog could be found at
    http://www.gaoran.ru/english/as/ac_vlbi/#OCARS.

    Returns
    -------
    ocars: an astropy.Table object
        data in the catalog

    """

    if ocars_file is None:
        data_dir = get_data_dir()
        ocars_file = "{}/ocars/ocars.csv".format(data_dir)

    ocars = Table.read(ocars_file, format="ascii.csv",
                       names=["iers_name", "iau_name",
                              "ra", "dec", "ra_err", "dec_err", "ra_dec_corr", "pos_epoch",
                              "ga_ra", "ga_dec", "gal_lon", "gal_lat", "class",
                              "z", "z_flag", "z_simbad", "z_simbad_flag", "z_sdss",
                              "u_mag", "U_mag", "B_mag", "g_mag", "V_mag", "r_mag",
                              "R_mag", "i_mag", "I_mag", "z_mag", "J_mag", "H_mag",
                              "K_mag", "G_mag"])

    # Add Correction factor to 'ra_err'
    arc_fac = np.cos(np.deg2rad(ocars["dec"]))
    ocars["ra_err"] = ocars["ra_err"] * arc_fac

    # Add unit information
    ocars["ra"].unit = u.deg
    ocars["dec"].unit = u.deg
    ocars["ra_err"].unit = u.mas
    ocars["dec_err"].unit = u.mas

    return ocars

# --------------------------------- END --------------------------------
