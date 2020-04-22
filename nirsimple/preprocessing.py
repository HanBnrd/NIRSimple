"""
NIRSimple
---------
NIRS processing made simple
"""

import warnings

import numpy as np


class Preprocessing():

    def optical_densities(self, intensities, refs=None):
        """
        Converts intensities into optical density changes. Changes are relative
        to the average intensity for each channel.

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
            intensities for each channel, of shape (channels, data points).
        """
        if refs is None:
            means = np.mean(np.absolute(intensities), axis=1)
            means = np.expand_dims(means, axis=1)
            delta_od = -np.log10(np.absolute(intensities)/means)
        else:
            references = np.expand_dims(refs, axis=1)
            delta_od = -np.log10(np.absolute(intensities)/references)

        if np.any(intensities <= 0):
            warnings.warn("some intensities are negative or equal to zero")
        return delta_od
