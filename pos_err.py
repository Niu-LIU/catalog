#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: pos_err.py
"""
Created on Fri Sep 21 15:36:35 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np
from numpy import sqrt, sin, cos
from numpy.linalg import eig
import astropy


__all__ = ["pos_err_calc", "overall_error", "error_ellipse_calc"]


# -----------------------------  FUNCTIONS -----------------------------
def pos_err_calc(ra_err, dec_err, ra_dec_corr):
    """Calculate the semi-major axis of the dispersion ellipse.

    Parameters
    ----------
    ra_err/dec_err : formal uncertainty of RA/Dec, usually in micro-as
    ra_dec_corr : correlation coeffient between RA and Dec, unitless.

    Returns
    ----------
    sig_pos_max : semi-major axis of the dispersion ellipse for
                  characterising the positional uncertainty of a source;
                  same unit as ra_err/dec_err
    """

    sig_pos_max = sqrt(0.5 * (ra_err**2 + dec_err**2 +
                              sqrt((ra_err**2 - dec_err**2)**2 +
                                   (2*ra_err*dec_err*ra_dec_corr)**2)))

    return sig_pos_max


def overall_error(ra_err, dec_err, ra_dec_corr):
    """Calculate the ovrall formal uncertainty.

    ovrall formal uncertainty = sqrt(RA_err^2+Dec_err^2+C*RA_err*Dec_err)

    Parameters
    ----------
    ra_err/dec_err : formal uncertainty of RA/Dec, usually in micro-as
    ra_dec_corr : correlation coeffient between RA and Dec, unitless.

    Returns
    ----------
    overall_err : ovrall formal uncertainty;
                  same unit as ra_err/dec_err
    """

    overall_err = sqrt(ra_err**2 + dec_err**2 + ra_err*dec_err*ra_dec_corr)

    return overall_err


def eepa_calc(ra_err, dec_err, ra_dec_corr, anticw=False):
    """Calculate the ovrall formal uncertainty.

    ovrall formal uncertainty = sqrt(RA_err^2+Dec_err^2+C*RA_err*Dec_err)

    Parameters
    ----------
    ra_err/dec_err : formal uncertainty of RA/Dec, usually in micro-as
    ra_dec_corr : correlation coeffient between RA and Dec, unitless.

    Returns
    ----------
    pa : the position angle of the major axis of the error ellipse
    """

    M, m, pa = error_ellipse(ra_err, dec_err, ra_dec_corr, anticw)

    return pa


def eepa_calc2(ra_err, dec_err, ra_dec_corr, anticw=False):
    """Calculate the ovrall formal uncertainty.

    ovrall formal uncertainty = sqrt(RA_err^2+Dec_err^2+C*RA_err*Dec_err)

    Parameters
    ----------
    ra_err/dec_err : formal uncertainty of RA/Dec, usually in micro-as
    ra_dec_corr : correlation coeffient between RA and Dec, unitless.

    Returns
    ----------
    pa : the position angle of the major axis of the error ellipse
    """

    pa = np.zeros_like(ra_err)
    for i in range(len(ra_err)):
        pa[i] = eepa_calc1(ra_err[i], dec_err[i], ra_dec_corr[i], anticw)

    return pa


def ellipse_shape_calc(M, m, pa, start_from_xaxis=False):
    """Calculate the (x, y) position for an ellipse after rotation.

    Parameters
    ----------
    M : float
        major axis
    m : floar
        minor axis
    pa : float
        positional angle, degree reckoned from x- or y-axis.

    Returns
    -------
    x1/y1 : ndarray of float
    """

    t = np.linspace(0, 2 * np.pi, 360)
    x = M * cos(t)
    y = m * sin(t)

    if start_from_xaxis:
        alpha = np.deg2rad(pa)
    else:
        alpha = np.deg2rad(90-pa)

    x1 = x * cos(alpha) - y * sin(alpha)
    y1 = x * sin(alpha) + y * cos(alpha)

    return x1, y1


def error_ellipse_calc_for_single(ra_err, dec_err, ra_dec_corr, anticw=False):
    """Calculate the ovrall formal uncertainty.

    ovrall formal uncertainty = sqrt(RA_err^2+Dec_err^2+C*RA_err*Dec_err)

    Parameters
    ----------
    ra_err/dec_err : formal uncertainty of RA/Dec, usually in micro-as
    ra_dec_corr : correlation coeffient between RA and Dec, unitless.

    Returns
    ----------
    M : semi-major axis of the error ellipse
    m : semi-minor axis of the error ellipse
    pa : the position angle of the major axis of the error ellipse
    """

    cov = ra_dec_corr * ra_err * dec_err
    cov_mat = np.array([[ra_err**2, cov], [cov, dec_err**2]])

    eig_val, eig_vec = eig(cov_mat)

    if eig_val[0] > eig_val[1]:
        M2 = eig_val[0]
        m2 = eig_val[1]
        vec_M = eig_vec[:, 0]
    else:
        M2 = eig_val[1]
        m2 = eig_val[0]
        vec_M = eig_vec[:, 1]

    M, m = np.sqrt(M2), np.sqrt(m2)

    # Calculate the position angle counted anti-colockwise from the declination axis
    if vec_M[1] == 0:
        if vec_M[0] > 0:
            pa0 = 90
        else:
            pa0 = -90
    else:
        pa0 = np.rad2deg(np.arctan(vec_M[0]/vec_M[1]))

    if anticw:
        if pa0 <= 0:
            pa = -pa0
        else:
            pa = 180 - pa0
    else:
        if pa0 <= 0:
            pa = 180 + pa0
        else:
            pa = pa0

    return M, m, pa


def error_ellipse_calc_for_single2(ra_err, dec_err, ra_dec_corr):
    """Calculate the ovrall formal uncertainty.

    ovrall formal uncertainty = sqrt(RA_err^2+Dec_err^2+C*RA_err*Dec_err)

    Parameters
    ----------
    ra_err/dec_err : formal uncertainty of RA/Dec, usually in micro-as
    ra_dec_corr : correlation coeffient between RA and Dec, unitless.

    Returns
    ----------
    eema : semi-major axis of the error ellipse
    eena : semi-minor axis of the error ellipse
    pa : the position angle of the major axis of the error ellipse
    """

    a = ra_err
    b = dec_err
    r = ra_dec_corr
    rab = r * a * b

    # This angle is counter-clockwise angle measured from the X-axis
    theta = 0.5 * np.arctan2(2 * rab, (a**2-b**2))

    sx = a**2 * cos(theta)**2 + b**2 * sin(theta)**2 + \
        2 * rab * sin(theta) * cos(theta)
    sy = a**2 * sin(theta)**2 + b**2 * cos(theta)**2 - \
        2 * rab * sin(theta) * cos(theta)
    eema = sqrt(np.maximum(sx, sy))
    eena = sqrt(np.minimum(sx, sy))

    # For the position angle which is supposed be eastward from the north
    # where north is the Y-axis, i.e., clockwise angle measured from Y-axis.
    # (The RA direction is the east while the declination is the north.)
    if theta < 0:
        theta0 = 360 + theta
    if a >= b:
        pa = (360 - (theta0 - 90)) % 360
    else:
        pa = 360 - theta0

    return eema, eena, pa


def error_ellipse_calc(ra_err, dec_err, ra_dec_corr, anticw=False):
    """Calculate Parameters of error ellipse for a subset of objects

    Parameters
    ----------
    ra_err/dec_err : 1-d array
        formal uncertainty of RA/Dec
    ra_dec_corr : 1-d array
        correlation coeffient between RA and Dec, unitless.

    Returns
    ----------
    M: 1-d array
        semi-major axis of the error ellipse
    m: 1-d array
        semi-minor axis of the error ellipse
    pa: 1-d array
        the position angle of the major axis of the error ellipse
    """

    if type(ra_err) == np.ndarray:
        # Here I just check the type of ra_err and assume that ra_err, dec_err, and
        # ra_dec_corr belongs to the same data type and has the same shape.
        # I can do this just because this script is only used by myself
        M = np.zeros_like(ra_err)
        m = np.zeros_like(dec_err)
        pa = np.zeros_like(ra_dec_corr)
        for i in range(len(M)):
            M[i], m[i], pa[i] = error_ellipse_calc_for_single(
                ra_err[i], dec_err[i], ra_dec_corr[i], anticw)

    elif type(ra_err) in [astropy.table.column.Column, astropy.table.column.MaskedColumn]:
        # Or Pandas type?
        # astropy.table -> np.ndarray
        ra_err = np.array(ra_err)
        dec_err = np.array(dec_err)
        ra_dec_corr = np.array(ra_dec_corr)

        M = np.zeros_like(ra_err)
        m = np.zeros_like(dec_err)
        pa = np.zeros_like(ra_dec_corr)

        for i in range(len(M)):
            M[i], m[i], pa[i] = error_ellipse_calc_for_single(
                ra_err[i], dec_err[i], ra_dec_corr[i], anticw)

    else:
        M, m, pa = error_ellipse_calc_for_single(
            ra_err, dec_err, ra_dec_corr, anticw)

    return M, m, pa


def main():
    """Main function.
    """

    print("This code is just a module of functions.")


if __name__ == "__main__":
    main()
# --------------------------------- END --------------------------------
