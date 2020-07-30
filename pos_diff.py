#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: pos_diff.py
"""
Created on Fri Sep 21 15:39:02 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""

from astropy.coordinates import SkyCoord
from astropy.table import Table, join, Column
from astropy import units as u
from functools import reduce
import numpy as np
from numpy import cos, sqrt
from math import hypot


__all__ = ["nor_sep_calc", "pos_diff_calc", "pa_calc",
           "radio_cat_diff_calc", ]


# -----------------------------  FUNCTIONS -----------------------------
def nor_sep_calc(dRA, dRA_err, dDC, dDC_err, C):
    '''Calculate the normalized seperation.

    Parameters
    ----------
    dRA/dDC : Right Ascension / Declination differences in micro-as
    dRA_err/dDC_err : formal uncertainty of dRA*cos(Dec)/dDC in micro-as
    C : correlation coeffient between dRA*cos(Dec) and dDC.

    Returns
    ----------
    ang_sep : angular seperation, in micro-as
    X_a / X_d : normalized coordinate differences in RA / DC, unit-less
    X : Normalized separations, unit-less.
    '''

    # Angular seperations
    ang_sep = sqrt(dRA**2 + dDC**2)
    # ang_sep = hypot(dRA, dDC) # Only used for scalar

    # Normalised coordinate difference
    X_a = dRA / dRA_err
    X_d = dDC / dDC_err

    # Normalised separation - Mignard's statistics (considering covariance)
    X = np.zeros_like(X_a)

    for i, (X_ai, X_di, Ci) in enumerate(zip(X_a, X_d, C)):
        if Ci == -1.:
            Ci = -0.999
        if Ci == 1.:
            Ci = 0.999

        wgt = np.linalg.inv(np.mat([[1, Ci], [Ci, 1]]))
        Xmat = np.mat([X_ai, X_di])
        X[i] = sqrt(reduce(np.dot, (Xmat, wgt, Xmat.T)))

    return ang_sep, X_a, X_d, X


def pos_diff_calc(RA1, RA1_err, DC1, DC1_err, Cor1,
                  RA2, RA2_err, DC2, DC2_err, Cor2,
                  arccof=None):
    '''Calculate the normalized seperation between VLBI and Gaia positions.


    Parameters
    ----------
    RA / DC : Right Ascension / Declination, degress
    e_RA / e_DC : formal uncertainty of RA * cos(Dec) / DC, mas
    Cor : correlation coeffient between RA and DC.
    arccof : cos(Dec.)

    Note: suffix 'G' stands for GaiaDR1 and I for VLBI catalog.

    Returns
    ----------
    ang_sep : angular seperation in micro-as
    X_a / X_d : normalized seperation in RA / DC, unit-less
    X : Normalized separations, unit-less.
    '''

    if arccof is None:
        arccof = cos(np.deg2rad(DC1))

    # # deg -> uas
    # dRA = (RA1 - RA2) * 3.6e9 * arccof
    # dRA_err = sqrt(RA1_err**2 + RA2_err**2)
    # dDC = (DC1 - DC2) * 3.6e9
    # dDC_err = sqrt(DC1_err**2 + DC2_err**2)

    # deg -> mas
    dRA = (RA1 - RA2) * 3.6e6 * arccof
    dRA_err = sqrt(RA1_err**2 + RA2_err**2)
    dDC = (DC1 - DC2) * 3.6e6
    dDC_err = sqrt(DC1_err**2 + DC2_err**2)

    # Correlation coefficient of combined errors
    cov = RA1_err * DC1_err * Cor1 + RA2_err * DC2_err * Cor2
    corf = cov / (dRA_err * dDC_err)

    # Normalised separation
    ang_sep, X_a, X_d, X = nor_sep_calc(dRA, dRA_err,
                                        dDC, dDC_err, corf)

    # return ang_sep, X_a, X_d, X
    return dRA, dDC, dRA_err, dDC_err, cov, ang_sep, X_a, X_d, X


def pa_calc0(dra, ddec, anticw=False):
    """Calculate positional angle from positional offset.

    Ax is used for plot and Ay for output value (positional angle).


    Parametes
    ---------
    dra : ndarray
        positional difference in R.A.(times cos(decl.))
    ddec : ndarray
        positional difference in declination
    anticw : Boolean


    Returns
    -------
    Ax : ndarray
        Angle (in degree) of positional offset vector towards to x-axis anti-clockwisely
    Ay : ndarray
        Angle (in degree) of positional offset vector towards to y-axis anti-clockwisely
    """

    Ax = np.rad2deg(np.arctan2(ddec, dra))  # clockwise
    Ay = np.rad2deg(np.arctan2(dra, ddec))  # clockwise

    if anticw:
        # anticlockwise
        Ax = np.where(Ax < 0, 360 + Ax, Ax)
        Ay = np.where(Ay < 0, -Ay, 360 - Ay)
    else:
        # clockwise
        Ax = np.where(Ax < 0, -Ax, 360 - Ax)
        Ay = np.where(Ay < 0, 360 + Ay, Ay)

    return Ax, Ay


def pa_calc1(dra, ddec, anticw=False):
    """Calculate positional angle from positional offset.

    Ax is used for plot and Ay for output value (positional angle).

    A new implementation.

    Parametes
    ---------
    dra : ndarray
        positional difference in R.A.(times cos(decl.))
    ddec : ndarray
        positional difference in declination
    anticw : Boolean


    Returns
    -------
    Ax : ndarray
        Angle (in degree) of positional offset vector towards to x-axis anti-clockwisely
    Ay : ndarray
        Angle (in degree) of positional offset vector towards to y-axis anti-clockwisely
    """

    if anticw:
        # anticlockwise
        zx = dra + 1j * ddec
        zy = ddec - 1j * dra

        Ax = np.angle(zx, deg=True)
        Ay = np.angle(zy, deg=True)
    else:
        # clockwise
        zx = dra - 1j * ddec
        zy = ddec + 1j * dra

        Ax = np.angle(zx, deg=True)
        Ay = np.angle(zy, deg=True)

    Ax = np.where(Ax < 0, 360 + Ax, Ax)
    Ay = np.where(Ay < 0, 360 + Ay, Ay)

    return Ax, Ay


def pa_calc(dra, ddec, anticw=False):
    """Calculate positional angle from positional offset.

    A new implementation.

    Parametes
    ---------
    dra : ndarray
        positional difference in R.A.(times cos(decl.))
    ddec : ndarray
        positional difference in declination
    anticw : Boolean


    Return
    -------
    PA : ndarray
        Angle (in degree) of positional offset vector towards to y-axis anti-clockwisely
    """

    cen = SkyCoord(0*u.deg, 0*u.deg, frame="icrs")
    oft = SkyCoord(dra*u.mas, ddec*u.mas, frame="icrs")

    pa = cen.position_angle(oft)
    pa = pa.to(u.deg)

    return pa.value


def radio_cat_diff_calc(table1, table2, sou_name="source_name",
                        label=["1", "2"]):
    """Calculate the positional differences between two radio source catalogs.

    These two catalogs should be stored in the astropy.table.Table object.
    And they should contain the following columns:
        -- "source_name"
        -- "ra"
        -- "dec"
        -- "ra_err"
        -- "dec_err"
        -- "ra_dec_corr"

    This function is written just for radio source catalogs.

    Parameters
    ----------
    table1/table2 : astropy.table.Table object
        Table of radi source position

    Returns
    -------
    com_sou : astropy.table.Table object
        Table of positional difference
    """

    # Table label
    label1, label2 = label

    # Keep only position-related columns
    table3 = Table(table1)
    table4 = Table(table2)

    table3.keep_columns([sou_name, "ra", "dec", "ra_err", "dec_err",
                         "ra_dec_corr", "pos_err"])
    table4.keep_columns([sou_name, "ra", "dec", "ra_err", "dec_err",
                         "ra_dec_corr", "pos_err"])

    # Cross-match by the source name
    com_sou = join(table3, table4, keys=sou_name, table_names=label)

    # Cos(decl.)
    arc_fac = cos(com_sou["dec_{}".format(label2)].to(u.rad).value)
    dra = (com_sou["ra_%s" % label1] - com_sou["ra_%s" % label2]) * arc_fac
    ddec = com_sou["dec_%s" % label1] - com_sou["dec_%s" % label2]

    dra_err = sqrt(com_sou["ra_err_%s" % label1]**2 +
                   com_sou["ra_err_%s" % label2]**2)
    ddec_err = sqrt(com_sou["dec_err_%s" % label1] ** 2 +
                    com_sou["dec_err_%s" % label2]**2)

    # Correlation coefficient of combined errors
    c1 = com_sou["ra_err_%s" % label1] * \
        com_sou["dec_err_%s" % label1] * com_sou["ra_dec_corr_%s" % label1]
    c2 = com_sou["ra_err_%s" % label2] * \
        com_sou["dec_err_%s" % label2] * com_sou["ra_dec_corr_%s" % label2]
    cov = c1 + c2
    corf = cov / (dra_err * ddec_err)

    # Add these columns
    com_sou.add_columns([dra, ddec, dra_err, ddec_err, cov],
                        names=["dra", "ddec", "dra_err", "ddec_err",
                               "dra_ddec_cov"])

    # Add unit information
    pos_unit = com_sou["ra_err_%s" % label1].unit
    com_sou["dra"].unit = u.deg
    com_sou["ddec"].unit = u.deg
    com_sou["dra"] = com_sou["dra"].to(pos_unit)
    com_sou["ddec"] = com_sou["ddec"].to(pos_unit)

    com_sou["dra"].unit = pos_unit
    com_sou["ddec"].unit = pos_unit

    com_sou["dra_err"].unit = pos_unit
    com_sou["ddec_err"].unit = pos_unit

    com_sou["dra_ddec_cov"].unit = None

    # Normalised separation
    ang_sep, X_a, X_d, X = nor_sep_calc(
        com_sou["dra"], com_sou["dra_err"],
        com_sou["ddec"], com_sou["ddec_err"], corf)

    # Direction of position offset
    pax, pay = pa_calc1(dra, ddec)
    pa = Column(pay, unit=u.deg)

    # Add these columns
    com_sou.add_columns([ang_sep, pa, X_a, X_d, X],
                        names=["ang_sep", "pa", "nor_ra", "nor_dec", "nor_sep"])

    com_sou["ang_sep"].unit = pos_unit
    com_sou["nor_ra"].unit = None
    com_sou["nor_dec"].unit = None
    com_sou["nor_sep"].unit = None

    # Remove some columns
    com_sou.remove_columns(["ra_%s" % label1, "dec_%s" % label1,
                            "ra_dec_corr_%s" % label1,
                            "ra_dec_corr_%s" % label2])

    com_sou.rename_column("ra_%s" % label2, "ra")
    com_sou.rename_column("dec_%s" % label2, "dec")

    return com_sou


def main():
    """Main function.
    """

    print("This code is just a module of functions.")


if __name__ == "__main__":
    main()
# --------------------------------- END --------------------------------
