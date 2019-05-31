#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: test.py
"""
Created on Sun Mar  3 10:53:06 2019

@author: Neo(liuniu@smail.nju.edu.cn)
"""

from astropy.table import Table, join, Column
from astropy import units as u
from functools import reduce
import numpy as np
from numpy import cos, sqrt


# -----------------------------  FUNCTIONS -----------------------------
mask1 = icrf3xka["type"] == "D"
mask2 = icrf3xka["type"] != "D"

icrf3xka_def = icrf3xka[mask1]
icrf3xka_oth = icrf3xka[mask2]


fig, (ax0, ax1) = plt.subplots(figsize=(10, 3), ncols=2)

ax0.plot(icrf3xka_oth["nb_sess"], icrf3xka_oth["pos_err"], "b*", label="Other")
ax0.plot(icrf3xka_def["nb_sess"], icrf3xka_def["pos_err"],
         "rx", label="Definingsx")
ax0.set_yscale("log")
ax0.set_xscale("log")
ax0.set_xlabel("N$_{sess}$")
ax0.legend()

ax1.plot(icrf3xka_oth["nb_del"], icrf3xka_oth["pos_err"], "b*", label="Other")
ax1.plot(icrf3xka_def["nb_del"],
         icrf3xka_def["pos_err"], "rx", label="Defining")
ax1.set_yscale("log")
ax1.set_xscale("log")
ax1.set_xlabel("N$_{delay}$")
ax1.legend()


if __name__ == '__main__':
    main()
# --------------------------------- END --------------------------------
