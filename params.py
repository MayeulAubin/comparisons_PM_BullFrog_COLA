#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# -------------------------------------------------------------------------------------
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# The text of the license is located in the root directory of the source package.
# -------------------------------------------------------------------------------------

__author__ = "Tristan Hoellinger"
__version__ = "0.1"
__date__ = "2024"
__license__ = "GPLv3"

"""
"""

WhichSpectrum = "EH"  # "EH" or "class"

# TODO: cosmo_planck_wo_kmax and add k_max later in compare_PM_COLA_BullFrog_parser.py
cosmo_planck = {
    "h": 0.6774,
    "Omega_m": 0.3089,
    "Omega_b": 0.0486,
    "n_s": 0.9667,
    "sigma8": 0.8159,
    "k_max": 10.0,
    "WhichSpectrum": WhichSpectrum,
}

# TODO: same as above + fix value for Omega_b
cosmo_BF_article = {
    "h": 0.677,
    "Omega_m": 0.302,
    "Omega_b": 0.001,
    "n_s": 0.968,
    "sigma8": 0.815,
    "k_max": 10.0,
    "WhichSpectrum": WhichSpectrum,
}

# cosmo = cosmo_planck
cosmo = cosmo_BF_article


def z2a(z):
    return 1.0 / (1 + z)


def cosmo_small_to_full_dict(cosmo_min):
    """Return a full cosmology dictionary from a minimal one.

    Parameters
    ----------
    cosmo_min : dict
        Minimal cosmology dictionary.

    Returns
    -------
    cosmo_full : dict
        Full cosmology dictionary.

    """
    cosmo_full = {
        "h": cosmo_min["h"],
        "Omega_r": 0.0,
        "Omega_q": 1 - cosmo_min["Omega_m"],
        "Omega_b": cosmo_min["Omega_b"],
        "Omega_m": cosmo_min["Omega_m"],
        "m_ncdm": 0.0,
        "Omega_k": 0.0,
        "tau_reio": 0.066,
        "n_s": cosmo_min["n_s"],
        "sigma8": cosmo_min["sigma8"],
        "w0_fld": -1.0,
        "wa_fld": 0.0,
        "k_max": cosmo_min["k_max"],
        "WhichSpectrum": cosmo_min["WhichSpectrum"],
    }
    return cosmo_full


BASELINE_SEEDPHASE = 300030898
