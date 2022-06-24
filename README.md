# NIRSimple

> fNIRS signal processing simplified

NIRSimple is a Python 3 package for fNIRS with a focus on signal processing. It gives control on preprocessing and processing, enabling to adjust parameters and use different methods from the scientific literature.

This package handles data as numpy arrays. It can be used with [MNE](https://mne.tools/stable/index.html) pipelines by creating MNE-specific objects (Raw, Epochs, ...) from numpy arrays.


### Features

- conversion from **light intensity** to **optical density changes**
- conversion from **optical density changes** to **hemoglobin concentration changes** with the **modified Beer-Lambert law** (options for different extinction coefficient tables)
- signal correction with **correlation based signal improvement** (CBSI)


### Installation

In a terminal or a command prompt, run:

```
pip install nirsimple
```


### Keep the package up to date

In a terminal or a command prompt, run:

```
pip install --upgrade nirsimple
```


### Tutorial

An example of using NIRSimple with MNE can be found [here](./examples/simple_probe.ipynb).
