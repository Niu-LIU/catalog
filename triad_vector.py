#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: triad_vector.py
"""
Created on Mon Nov 16 16:43:06 2020

@author: Neo(niu.liu@nju.edu.cn)

Generate the triad vectors of the coordinate system
"""

import numpy as np
from numpy import sin, cos


# -----------------------------  MAIN -----------------------------
def p_vector(ra, dec, unit_deg=False):
    """Return unit vector in RA

    vec{p} = [-sin(ra), cos(ra), 0]^T

    Parameters
    ----------
    ra/dec : 1-d array 
            right ascension (RA) /declination (Dec) in degree or radian
    unit_deg : Boolean
            tell if the RA/Dec is given in degree (True for degree)
    """

    if unit_deg:
        ra = np.deg2rad(ra)
        deg = np.deg2rad(deg)

    # Check if ra and dec is in same length
    if not len(ra) == len(dec):
        print("Length of ra ({}) do not equal to that of dec ({}).".format(
            len(ra), len(dec)))

    N = len(ra)
    p1 = -sin(ra)
    p2 = cos(ra)
    p3 = np.zeros(N)

    p_vec = np.concatenate((p1.reshape(N, 1),
                            p2.reshape(N, 1),
                            p3.reshape(N, 1)), axis=1)

    return p_vec


def q_vector(ra, dec, unit_deg=False):
    """Return unit vector in declination

    vec{q} = [-sin(dec)*cos(ra), -sin(dec)*sin(ra), cos(dec)]^T

    Parameters
    ----------
    ra/dec : 1-d array 
            right ascension (RA) /declination (Dec) in degree or radian
    unit_deg : Boolean
            tell if the RA/Dec is given in degree (True for degree)
    """

    if unit_deg:
        ra = np.deg2rad(ra)
        deg = np.deg2rad(deg)

    # Check if ra and dec is in same length
    if not len(ra) == len(dec):
        print("Length of ra ({}) do not equal to that of dec ({}).".format(
            len(ra), len(dec)))

    N = len(ra)
    q1 = -sin(dec) * cos(ra)
    q2 = -sin(dec) * sin(ra)
    q3 = cos(dec)

    q_vec = np.concatenate((q1.reshape(N, 1),
                            q2.reshape(N, 1),
                            q3.reshape(N, 1)), axis=1)

    return q_vec


def r_vector(ra, dec, unit_deg=False):
    """Return unit vector in declination

    vec{q} = [-sin(dec)*cos(ra), -sin(dec)*sin(ra), cos(dec)]^T

    Parameters
    ----------
    ra/dec : 1-d array 
            right ascension (RA) /declination (Dec) in degree or radian
    unit_deg : Boolean
            tell if the RA/Dec is given in degree (True for degree)
    """

    if unit_deg:
        ra = np.deg2rad(ra)
        deg = np.deg2rad(deg)

    # Check if ra and dec is in same length
    if not len(ra) == len(dec):
        print("Length of ra ({}) do not equal to that of dec ({}).".format(
            len(ra), len(dec)))

    N = len(ra)
    r1 = cos(dec) * cos(ra)
    r2 = cos(dec) * sin(ra)
    r3 = sin(dec)

    r_vec = np.concatenate((r1.reshape(N, 1),
                            r2.reshape(N, 1),
                            r3.reshape(N, 1)), axis=1)

    return r_vec
# --------------------------------- END --------------------------------
