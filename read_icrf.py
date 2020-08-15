#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: read_icrfn.py
"""
Created on Sat Sep 29 18:15:50 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""

from astropy.table import Table, Column
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np
import sys
from sys import platform as _platform

# My modules
from .pos_err import error_ellipse_array

__all__ = ["read_icrf1", "read_icrf2", "read_icrf3"]


# -----------------------------  FUNCTIONS -----------------------------
def get_datadir():
    """Check the OS and get the data directory.

    I have several computers with different OS.
    The path to the data file is different.

    Returns
    -------
    datadir: string
        directory which stores ICRF catalogs.

    """

    # Check the type of OS
    if _platform == "linux" or _platform == "linux2":
        # linux
        datadir = "/home/neo/Astronomy/data/catalogs/icrf"
    elif _platform == "darwin":
        # MAC OS X
        datadir = "/Users/Neo/Astronomy/data/catalogs/icrf"
    elif _platform == "win32":
        # Windows
        print("Not implemented yet")
        exit()
    elif _platform == "win64":
        # Windows 64-bit
        print("Not implemented yet")
        exit()
    else:
        print("Weird! What kind of OS do you use?")
        exit()

    return datadir


def read_icrf1(icrf1file=None):
    """Read the ICRF1 catalog

    Parameter
    ---------
    icrf1file : string
        file name and path of the ICRF1 catalog

    Return
    ------
    icrf1 : an astropy.Table object
        data in the catalog
    """

    if icrf1file is None:
        datadir = get_datadir()
        icrf1file = "{}/rsc95r01.dat".format(datadir)

    # Read ICRF1 catalog
    icrf1 = Table.read(icrf1file,
                       format="ascii.fixed_width_no_header",
                       names=["icrf_name", "iers_name", "type", "si_s", "si_x",
                              "ra_err", "dec_err", "ra_dec_corr",
                              "mean_obs", "beg_obs", "end_obs",
                              "nb_sess", "nb_del"],
                       col_starts=[5, 24, 34, 35, 37, 77, 87,
                                   96, 102, 112, 122, 132, 138],
                       col_ends=[21, 32, 34, 35, 37, 84, 93,
                                 100, 110, 120, 130, 136, 143])

    # Position information
    ra_dec_str = Table.read(icrf1file,
                            format="ascii.fixed_width_no_header",
                            names=["ra_dec"], col_starts=[42], col_ends=[73])

    ra_dec = SkyCoord(ra_dec_str["ra_dec"], unit=(u.hourangle, u.deg))
    ra = Column(ra_dec.ra, name="ra")
    dec = Column(ra_dec.dec, name="dec")

    # Add source position to the table
    icrf1.add_columns([ra, dec], indexes=[3, 3])

    # Add unit information
    icrf1["ra_err"] = icrf1["ra_err"] * 15e3 * np.cos(ra_dec.dec.rad)
    icrf1["ra_err"].unit = u.mas
    icrf1["dec_err"].unit = u.arcsec
    icrf1["dec_err"] = icrf1["dec_err"].to(u.mas)
    icrf1["ra_dec_corr"] = icrf1["ra_dec_corr"].filled(0)

    # Calculate the semi-major axis of error ellipse
    pos_err, pos_err_min, pa = error_ellipse_array(
        icrf1["ra_err"], icrf1["dec_err"], icrf1["ra_dec_corr"])

    # Add the semi-major axis of error ellipse to the table
    pos_err = Column(pos_err, name="pos_err", unit=u.mas)
    pa = Column(pa, name="eepa", unit=u.deg)
    icrf1.add_columns([pos_err, pa], indexes=[9, 10])

    return icrf1


def read_icrf2(icrf2file=None):
    """Read the ICRF1 catalog

    Parameter
    ---------
    icrf2file : string
        file name and path of the ICRF2 catalog

    Return
    ------
    icrf2 : an astropy.Table object
        data in the catalog
    """

    if icrf2file is None:
        datadir = get_datadir()
        icrf2file = "{}/icrf2.dat".format(datadir)

    # Read ICRF2 catalog
    icrf2 = Table.read(icrf2file,
                       format="ascii.fixed_width_no_header",
                       names=["icrf_name", "ivs_name", "iers_name", "type",
                              "ra_err", "dec_err", "ra_dec_corr",
                              "mean_obs", "beg_obs", "end_obs",
                              "nb_sess", "nb_del"],
                       col_starts=[0, 17, 25, 35, 73, 84,
                                   94, 101, 109, 117, 125, 130],
                       col_ends=[15, 24, 33, 35, 82, 92,
                                 99, 107, 115, 123, 128, 135])

    # Position information
    ra_dec_str = Table.read(icrf2file,
                            format="ascii.fixed_width_no_header",
                            names=["ra_dec"], col_starts=[37], col_ends=[71])

    ra_dec = SkyCoord(ra_dec_str["ra_dec"], unit=(u.hourangle, u.deg))
    ra = Column(ra_dec.ra, name="ra")
    dec = Column(ra_dec.dec, name="dec")

    # Add source position to the table
    icrf2.add_columns([ra, dec], indexes=[3, 3])

    # Add unit information
    icrf2["ra_err"] = icrf2["ra_err"] * 15e3 * np.cos(ra_dec.dec.rad)
    icrf2["ra_err"].unit = u.mas
    icrf2["dec_err"].unit = u.arcsec
    icrf2["dec_err"] = icrf2["dec_err"].to(u.mas)

    # Calculate the semi-major axis of error ellipse
    pos_err, pos_err_min, pa = error_ellipse_array(
        icrf2["ra_err"], icrf2["dec_err"], icrf2["ra_dec_corr"])

    # Add the semi-major axis of error ellipse to the table
    pos_err = Column(pos_err, name="pos_err", unit=u.mas)
    pa = Column(pa, name="eepa", unit=u.deg)
    icrf2.add_columns([pos_err, pa], indexes=[9, 10])

    return icrf2


def read_icrf3(icrf3_file=None, wv="sx"):
    """Read the ICRF3 catalog.

    Parameter
    ---------
    icrf3_file : string
        file name and path of the ICRF3 catalog
    wv : string
        wavelength, could be "sx", "k", or "xka".

    Return
    ------
    icrf3 : an astropy.Table object
        data in the catalog
    """

    datadir = get_datadir()

    if icrf3_file is None:
        if wv == "sx" or wv == "SX":
            icrf3_file = "%s/icrf3sx.txt" % datadir
        elif wv == "k" or wv == "K":
            icrf3_file = "%s/icrf3k.txt" % datadir
        elif wv == "xka" or wv == "XKA":
            icrf3_file = "%s/icrf3xka.txt" % datadir
        else:
            print("wv could only be 'sx', 'k', or 'xka'.")
            sys.exit()

    icrf3 = Table.read(icrf3_file,
                       format="ascii.fixed_width", data_start=16,
                       names=["icrf_name", "iers_name", "type",
                              "ra_err", "dec_err", "ra_dec_corr",
                              "mean_obs", "beg_obs", "end_obs",
                              "nb_sess", "nb_del"],
                       col_starts=[5, 25, 35, 83, 98,
                                   108, 118, 127, 136, 145, 150],
                       col_ends=[20, 32, 35, 92, 106,
                                 114, 124, 133, 142, 148, 155])

    # Position information
    ra_dec_str = Table.read(icrf3_file,
                            format="ascii.fixed_width", data_start=16,
                            names=["ra_dec"], col_starts=[40], col_ends=[77])

    ra_dec = SkyCoord(ra_dec_str["ra_dec"], unit=(u.hourangle, u.deg))
    ra = Column(ra_dec.ra, name="ra")
    dec = Column(ra_dec.dec, name="dec")

    # Add source position to the table
    icrf3.add_columns([ra, dec], indexes=[3, 3])

    # Add unit information
    icrf3["ra_err"] = icrf3["ra_err"] * 15e3 * np.cos(ra_dec.dec.rad)
    icrf3["ra_err"].unit = u.mas
    icrf3["dec_err"].unit = u.arcsec
    icrf3["dec_err"] = icrf3["dec_err"].to(u.mas)

    # Calculate the semi-major axis of error ellipse
    pos_err, pos_err_min, pa = error_ellipse_array(
        icrf3["ra_err"], icrf3["dec_err"], icrf3["ra_dec_corr"])

    # Add the semi-major axis of error ellipse to the table
    pos_err = Column(pos_err, name="pos_err", unit=u.mas)
    pa = Column(pa, name="eepa", unit=u.deg)
    icrf3.add_columns([pos_err, pa], indexes=[9, 10])

    return icrf3


# -------------------------------- MAIN --------------------------------
if __name__ == '__main__':
    pass
# --------------------------------- END --------------------------------
