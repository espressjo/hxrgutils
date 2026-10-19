# setup.py
from setuptools import setup, find_packages

setup(
    name="hxrgutils",
    version="1.0",
    description="Utilities for raw H2RG/H4RG data: reads, reference pixels, "
                "up-the-ramp fitting and noise diagnostics (NIRPS, SPIRou, SPIP)",
    packages=find_packages(),
    python_requires=">=3.10",
    # Minimum versions are the ones the code is developed and tested with
    # (Python 3.12). numpy, scipy and astropy are the core; the others are
    # imported lazily inside the methods that use them, but are listed so
    # that every feature works after a plain install.
    install_requires=[
        "numpy>=1.26",
        "scipy>=1.15",
        "astropy>=6.1",
        "matplotlib>=3.10",   # show(), show_fft(), dist()
        "seaborn>=0.13",      # show(), show_fft()
        "tqdm>=4.67",         # progress bar in make_ramp()
        "natsort>=8.4",       # file ordering in hxramp.create_error()
    ],
)
