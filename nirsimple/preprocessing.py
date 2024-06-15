"""
The preprocessing module contains functions to preprocess fNIRS signals.
"""

import warnings
from os import path

import numpy as np
import pandas as pd
from scipy import interpolate


def _extinctions(wavelengths, table='wray', verbose=True):
    """
    Get molar extinction coefficients for oxygenated hemoglobin (HbO) and
    deoxygenated hemoglobin (HbR) corresponding to the wavelengths.

        Values for molar extinction coefficients in [cm-1/(moles/liter)] or
        [cm-1/M] based on the wavelength in [nm].

    Parameters
    ----------
    wavelengths : list of integers
        The two wavelengths for which to get the molar extinction
        coefficients (in nm).

    table : string
        Table to use as molar extinction coefficients.
        'wray': data from S. Wray et al., 1988
        'cope': data from M. Cope, 1991
        'gratzer': data from W.B. Gratzer and K. Kollias compiled by S. Prahl
        (https://omlc.org/spectra/hemoglobin/summary.html)
        'moaveni': data from M.K. Moaveni and J.M. Schmitt compiled by S. Prahl
        (https://omlc.org/spectra/hemoglobin/moaveni.html)
        'takatani': data from S. Takatani and M.D. Graham compiled by S. Prahl
        (https://omlc.org/spectra/hemoglobin/takatani.html)

    Returns
    -------
    ex : array
        numpy array of the extinction coefficients of shape (2, 2). Output is
        [[exHbO_1, exHbR_1], [exHbO_2, exHbR_2]], with [wl_1, wl2] as input.
    """
    citation = None
    if table == 'wray':
        citation = "S. Wray et al., 1988"
    elif table == 'cope':
        citation = "M. Cope, 1991"
    elif table == 'gratzer':
        citation = "W.B. Gratzer and K. Kollias compiled by S. Prahl"
    elif table == 'moaveni':
        citation = "M.K. Moaveni and J.M. Schmitt compiled by S. Prahl"
    elif table == 'takatani':
        citation = "S. Takatani and M.D. Graham compiled by S. Prahl"
    else:
        raise Exception("table unknown")
    ex = []
    if len(wavelengths) == 2 and wavelengths[0] != wavelengths[1]:
        ex_file = table + '.csv'
        ex_path = path.join(path.dirname(__file__), 'tables', ex_file)
        df = pd.read_csv(ex_path)
        wl = df['lambda'].to_numpy()
        hbo = df['hbo'].to_numpy()
        hbr = df['hbr'].to_numpy()
        interp_hbo = interpolate.interp1d(wl, hbo)
        interp_hbr = interpolate.interp1d(wl, hbr)
        for wavelength in wavelengths:
            try:
                ex.append([interp_hbo(wavelength), interp_hbr(wavelength)])
            except ValueError:
                raise Exception("no matching wavelength found")
    else:
        raise Exception("wavelengths should be 2 different values")

    if verbose is True:
        print("-----")
        print("Molar extinction coefficients (in cm-1/M):")
        print("{} nm | HbO: {}, HbR: {}".format(wavelengths[0], *ex[0]))
        print("{} nm | HbO: {}, HbR: {}".format(wavelengths[1], *ex[1]))
        print("(" + citation + ")")
        print("-----")

    ex = np.array(ex)
    return ex


def intensities_to_od_changes(intensities, refs=None):
    """
    Converts intensities into optical density changes. Changes are relative
    to the average intensity or a reference intensity for each channel.

        Optical density from light intensity:
        OD = log10(I_0/I_t)

        Optical density changes relative to average transmitted intensity:
        delta_OD = log10(I_0/I_t) - log10(I_0/I_average)
        delta_OD = -log10(I_t/I_average)

    Parameters
    ----------
    intensities : array
        numpy array of absolute intensities, must have the correct shape
        (channels, data points).

    refs : list of floats
        List of reference intensities to use instead of averages, length
        must be equal to the number of channels.

    Returns
    -------
    delta_od : array
        numpy array of optical density changes, relative to average
        intensities or a reference for each channel, of shape (channels, data
        points).
    """
    if refs is None:
        means = np.mean(np.absolute(intensities), axis=1)
        means = np.expand_dims(means, axis=1)
        delta_od = -np.log10(np.absolute(intensities)/means)
    else:
        references = np.expand_dims(refs, axis=1)
        delta_od = -np.log10(np.absolute(intensities)/references)
        if np.any(references <= 0):
            warnings.warn("some references are negative or equal to zero")

    if np.any(intensities <= 0):
        warnings.warn("some intensities are negative or equal to zero")
    return delta_od


def od_to_od_changes(optical_densities, refs=None):
    """
    Convert optical densities into optical density changes, relative to the
    average optical density or a reference optical density for each channel.

        Optical density changes relative to average optical density:
        delta_OD = OD - OD_average

    Parameters
    ----------
    optical_densities : array
        numpy array of optical densities, must have the correct shape
        (channels, data points).

    refs: list of floats
        List of reference optical densities to use instead of averages, length
        must be equal to the number of channels.

    Returns
    -------
    delta_od : array
        numpy array of optical density changes, relative to average optical
        densities or a reference for each channel, of shape (channels, data
        points).
    """
    if refs is None:
        means = np.mean(np.absolute(optical_densities), axis=1)
        means = np.expand_dims(means, axis=1)
        delta_od = np.absolute(optical_densities) - means
    else:
        references = np.expand_dims(refs, axis=1)
        delta_od = np.absolute(optical_densities) - references
        if np.any(references <= 0):
            warnings.warn("some references are negative or equal to zero")

    if np.any(optical_densities <= 0):
        warnings.warn("some optical densities are negative or equal to zero")
    return delta_od


def get_dpf(wavelength, age):
    """
    Get the differential pathlength factor (DPF) from wavelength and age using
    the general equation from Scholkmann & Wolf, 2013.

        General DPF equation:
        DPF = (223.3
        + 0.05624*age**0.8493
        - 5.723e-7*wavelength**3
        + 0.001245*wavelength**2
        - 0.9025*wavelength)

    Parameters
    ----------
    wavelength : float
        Wavelength in nm.

    age : float
        Participant's age in years.

    Returns
    -------
    dpf : float
        Differential pathlength factor (DPF).
    """
    dpf = (223.3
           + 0.05624*age**0.8493
           - 5.723e-7*wavelength**3
           + 0.001245*wavelength**2
           - 0.9025*wavelength)

    return dpf


def mbll(delta_od, ch_names, ch_wls, ch_dpfs, ch_distances, unit,
         table='wray'):
    """
    Apply the modified Beer-Lambert law (from Delpy et al., 1988) to optical
    density changes in order to obtain concentration changes in oxygenated
    hemoglobin (HbO) and deoxygenated hemoglobin (HbR).

        Modified Beer-Lambert law:
        Dod_wl = e_HbO_wl*Dc_HbO*l*DPF_wl + e_HbR_wl*Dc_HbR*l*DPF_wl

        With two different wavelengths we obtain a set of linear equations
        solved with matrices:
        [[Dod_1]  = [[e_HbO_1*l*DPF_1, e_HbR_1*l*DPF_1]  . [[Dc_HbO]
        [Dod_2]]    [e_HbO_2*l*DPF_2, e_HbR_2*l*DPF_2]]    [Dc_HbR]]

        Equivalent to:
        [[Dc_HbO]  = [[e_HbO_1, e_HbR_1] -1 . [[Dod_1/(l*DPF_1)]
        [DC_HbR]]    [e_HbO_2, e_HbR_2]]      [Dod_2/(l*DPF_2)]]

    Parameters
    ----------
    delta_od : array
        numpy array of optical density changes, relative to average
        intensities for each channel, of shape (channels, data points).

    ch_names : list of strings
        List of channel names.

    ch_wls : list of integers
        List of channel wavelengths (in nm).

    ch_dpfs : list of floats
        List of channel differential pathlength factors (DPF) (or partial
        pathlength factors (PPF)).

    ch_distances : list of floats
        List of channel source-detector distances.

    unit : string
        Unit for ch_distances ('cm' or 'mm').

    table : string
        Table to use as molar extinction coefficients.
        'wray': data from S. Wray et al., 1988
        'cope': data from M. Cope, 1991
        'gratzer': data from W.B. Gratzer and K. Kollias compiled by S. Prahl
        'moaveni': data from M.K. Moaveni and J.M. Schmitt compiled by S. Prahl
        'takatani': data from S. Takatani and M.D. Graham compiled by S. Prahl

    Returns
    -------
    delta_c : array
        numpy array of hemoglobin concentration changes in [moles/liter] or [M]
        for each channel, of shape (channels, data points).

    new_ch_names : list of strings
        New list of channel names.

    new_ch_types : list of strings
        New list of channel types ('hbo' for oxygenated hemoglobin and 'hbr'
        for deoxygenated hemoglobin).
    """
    if unit == 'cm':
        pass
    elif unit == 'mm':
        ch_distances = np.array(ch_distances) / 10
        ch_distances = ch_distances.tolist()
    else:
        raise Exception("unit should be cm or mm")
    delta_c = []
    new_ch_names = []
    new_ch_types = []
    for name in np.unique(ch_names):
        idx_1, idx_2 = [idx for idx, x in enumerate(ch_names) if x == name]

        sub_delta_od = [delta_od[idx_1], delta_od[idx_2]]
        sub_delta_od = np.swapaxes(sub_delta_od, 0, 1)  # timepoints first
        sub_delta_od = np.expand_dims(sub_delta_od, axis=2)

        ex = _extinctions([ch_wls[idx_1], ch_wls[idx_2]], table, verbose=False)
        ex_inv = np.linalg.inv(ex)
        ex_inv = np.tile(ex_inv, (len(sub_delta_od), 1, 1))

        pl = [ch_dpfs[idx_1]*ch_distances[idx_1],
              ch_dpfs[idx_2]*ch_distances[idx_2]]
        pl = np.tile(np.array(pl), (len(sub_delta_od), 1))
        pl = np.expand_dims(pl, axis=2)

        sub_delta_c = np.matmul(ex_inv, sub_delta_od/pl)
        sub_delta_c = np.swapaxes(sub_delta_c, 0, 1)  # channels first

        delta_c.append(sub_delta_c[0])
        delta_c.append(sub_delta_c[1])
        new_ch_names.append(name)
        new_ch_names.append(name)
        new_ch_types.append('hbo')
        new_ch_types.append('hbr')

    delta_c = np.array(delta_c)
    delta_c = np.squeeze(delta_c)
    return delta_c, new_ch_names, new_ch_types
