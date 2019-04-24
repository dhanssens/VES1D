#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
VES 1D Modeling
===============
Example file for VES1D modeling.

Mainly adapted from:
Ekinci, Y. L., Demirci, A., 2008. A Damped Least-Squares Inversion Program
for the Interpretation of Schlumberger Sounding Curves, Journal of Applied Sciences, 8, 4070-4078.
"""

# Import
import numpy as np
import VES1D
import matplotlib.pyplot as plt


###############################################################################
# Input
# -----
#

# Electrical conductivity of n layers (S/m), last layer n is assumed to be infinite
con = np.array([.1, .02, .1])  # S/m

# Thickness of n-1 layers (m), last layer n is assumed to be infinite and does not require a thickness
thick = np.array([2, 1])  # m

# Half the current (AB/2) electrode spacing (m)
ab2s = np.linspace(0, 5, 100)


###############################################################################
# Forward model
# -------------
#

# Calculate forward apparent electrical conductivities (ECa)
app_con = []
for ab2 in ab2s:
        app_con.append(VES1D.forward(con, thick, ab2))

# Calculate apparent resistivity (Ohm.m)
app_res = 1 / np.array(app_con)


###############################################################################
# Visualize
# ---------
#

# Figure
plt.figure()
plt.xlabel('AB/2 (m)')
ax1 = plt.gca()
plt.plot(ab2s, app_con, color='tab:blue')
ax1.set_ylabel('$\sigma_a$ (S/m)', color='tab:blue')
ax1.tick_params(axis='y', labelcolor='tab:blue')
ax2 = ax1.twinx()
ax2.plot(ab2s, app_res, color='tab:red')
ax2.set_ylabel('$\\rho_a$ ($\Omega$m)', color='tab:red')
ax2.tick_params(axis='y', labelcolor='tab:red')
