# NIRSimple

> fNIRS signal processing simplified

*NIRSimple* is a Python 3 package for fNIRS with a focus on signal processing. It gives control on preprocessing and processing, enabling to adjust parameters and use different methods from the scientific literature.

This package handles data as numpy arrays. It can be used with [MNE](https://mne.tools/stable/index.html) pipelines by creating MNE-specific objects (Raw, Epochs, ...) from numpy arrays.


[![DOI](https://img.shields.io/badge/doi-10.3389%2Ffnrgo.2023.994969-blue)](https://doi.org/10.3389/fnrgo.2023.994969)
[![License](https://img.shields.io/badge/license-MIT-lightgrey)](https://github.com/HanBnrd/NIRSimple/blob/master/LICENSE)
[![Pipeline](https://img.shields.io/github/actions/workflow/status/HanBnrd/NIRSimple/sphinx.yml?label=pipeline)](https://github.com/HanBnrd/NIRSimple)
[![PyPI version](https://img.shields.io/pypi/v/nirsimple)](https://pypi.org/project/nirsimple/)


### Features

- conversion from **light intensity** to **optical density changes**
- **differential pathlength factor** (DPF) from wavelength and age
- conversion from **optical density changes** to **hemoglobin concentration changes** with the **modified Beer-Lambert law** (options for different extinction coefficient tables)
- signal correction with **correlation based signal improvement** (CBSI)


### Documentation

The documentation can be found [here](https://hanbnrd.github.io/NIRSimple).


### Install

In a terminal or a command prompt, run:

```
pip install nirsimple
```


### Keep the package up to date

In a terminal or a command prompt, run:

```
pip install --upgrade nirsimple
```


### Example

An example of using *NIRSimple* with MNE can be found [here](https://hanbnrd.github.io/NIRSimple/examples/simple-probe.html).


### Acknowledgements

Until there is an article specifically on *NIRSimple*, please cite [this article](https://doi.org/10.3389/fnrgo.2023.994969).
