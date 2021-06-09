import sys

from skbuild import setup

setup(
    name="composite-measure",
    version="1.0.0",
    description="objective measurement of speech quality simulating the ITU P.835 testing metodology",
    author="Marcin Kuropatwiński",
    license="GNU 3.0",
    packages=["composite"],
    install_requires=["cython"],
    )
