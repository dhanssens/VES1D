#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Import
import numpy as np


def forward(con, thick, ab2):
    """
    Calculate forward VES response with half the current electrode spacing for a 1D layered earth.

    Parameters
    ----------
    con : np.array, (n,)
        Electrical conductivity of n layers (S/m), last layer n is assumed to be infinite

    thick : np.array, (n-1,)
        Thickness of n-1 layers (m), last layer n is assumed to be infinite and does not require a thickness

    ab2 : float
        Half the current (AB/2) electrode spacing (m)

    Returns
    -------
    app_con : float
        Apparent half-space electrical conductivity (S/m)

    References
    ----------
    Ekinci, Y. L., Demirci, A., 2008. A Damped Least-Squares Inversion Program 
    for the Interpretation of Schlumberger Sounding Curves, Journal of Applied Sciences, 8, 4070-4078.
	
    Koefoed, O., 1970. A fast method for determining the layer distribution from the raised
    kernel function in geoelectrical soundings, Geophysical Prospection, 18, 564-570.

    Nyman, D. C., Landisman, M., 1977. VES Dipole-dipole filter coefficients,
    Geophysics, 42(5), 1037-1044.
    """

    # Conductivity to resistivity and number of layers
    res = 1 / con
    lays = len(res) - 1

    # Constants
    LOG = np.log(10)
    COUNTER = 1 + (2 * 13 - 2)
    UP = np.exp(0.5 * LOG / 4.438)

    # Filter integral variable
    up = ab2 * np.exp(-10 * LOG / 4.438)

    # Initialize array
    ti = np.zeros(COUNTER)

    for ii in range(COUNTER):

        # Set bottom layer equal to its resistivity
        ti1 = res[lays]

        # Recursive formula (Koefoed, 1970)
        lay = lays
        while lay > 0:
            lay -= 1
            tan_h = np.tanh(thick[lay] / up)
            ti1 = (ti1 + res[lay] * tan_h) / (1 + ti1 * tan_h / res[lay])

        # Set overlaying layer to previous
        ti[ii] = ti1

        # Update filter integral variable
        up *= UP

    # Apply point-filter weights (Nyman and Landisman, 1977)
    res_a = 105 * ti[0] - 262 * ti[2] + 416 * ti[4] - 746 * ti[6] + 1605 * ti[8] - 4390 * ti[10] + 13396 * ti[12]
    res_a += - 27841 * ti[14] + 16448 * ti[16] + 8183 * ti[18] + 2525 * ti[20] + 336 * ti[22] + 225 * ti[24]
    res_a /= 1e4

    # Resistivity to conductivity
    return 1 / res_a

