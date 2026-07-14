# setup.py
from setuptools import setup, find_packages

setup(
    name="hxrgutils",
    version="1.0",
    packages=find_packages(),  # Automatically find your "mymodule"
    install_requires=["astropy>=5.3",
                      ""]
)
