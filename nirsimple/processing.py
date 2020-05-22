"""
The processing module contains functions to process fNIRS signals.
"""

import numpy as np


def cbsi(delta_c, ch_names, ch_types):
    """
    Apply correlation based signal improvement (from Cui at al., 2010) to
    hemoglobin concentration changes.

        Correlation based signal improvement (CBSI):
        x_0 = (1/2)*(x-alpha*x)
        y_0 = -(1/alpha)*x_0

    Parameters
    ----------
    delta_c : array
        numpy array of hemoglobin concentration changes in [moles/liter] or [M]
        for each channel, of shape (channels, data points).

    ch_names : list of strings
        List of channel names.

    ch_types : list of strings
        List of channel types ('hbo' for oxygenated hemoglobin and 'hbr' for
        deoxygenated hemoglobin).

    Returns
    -------
    delta_c_0 : array
        numpy array of corrected activation signals in [moles/liter] or [M] for
        each channel, of shape (channels, data points).

    new_ch_names : list of strings
        New list of channel names.

    new_ch_types : list of strings
        New list of channel types ('hbo' for oxygenated hemoglobin and 'hbr'
        for deoxygenated hemoglobin).
    """
    delta_c_0 = []
    new_ch_names = []
    new_ch_types = []
    for name in np.unique(ch_names):
        idx_hbo = [idx for idx, x in enumerate(ch_names)
                   if x == name and ch_types[idx] == 'hbo']
        idx_hbr = [idx for idx, x in enumerate(ch_names)
                   if x == name and ch_types[idx] == 'hbr']

        alpha = np.std(delta_c[idx_hbo]) / np.std(delta_c[idx_hbr])
        delta_c_0_hbo = (delta_c[idx_hbo] - alpha*delta_c[idx_hbr]) / 2
        delta_c_0_hbr = -delta_c_0_hbo / alpha

        delta_c_0.append(delta_c_0_hbo)
        delta_c_0.append(delta_c_0_hbr)
        new_ch_names.append(name)
        new_ch_names.append(name)
        new_ch_types.append('hbo')
        new_ch_types.append('hbr')

    delta_c_0 = np.array(delta_c_0)
    delta_c_0 = np.squeeze(delta_c_0)
    return delta_c_0, new_ch_names, new_ch_types
