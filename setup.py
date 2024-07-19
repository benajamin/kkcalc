import os
import sys

from setuptools import find_packages, setup


def read(rel_path: str) -> str:
    here = os.path.abspath(os.path.dirname(__file__))
    with open(os.path.join(here, rel_path)) as fp:
        return fp.read()


def get_version(rel_path: str) -> str:
    for line in read(rel_path).splitlines():
        if line.startswith("__version__"):
            # __version__ = "0.9"
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    raise RuntimeError("Unable to find version string.")

with open("README.rst", "r") as fh:
    long_description = fh.read()

setup(
    name="kkcalc",
    version=get_version("kkcalc/__init__.py"),
    author="Benjamin Watts",
    author_email="benjamin.watts@gmail.com",
    description="A calculator for the Kramers Kronig transform of X-ray spectra",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url="https://github.com/benajamin/kkcalc",
    packages=find_packages(),
    package_data={'kkcalc': ['*.json']},
    install_requires=['numpy', 'scipy'],
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: zlib/libpng License  ",
        "Topic :: Scientific/Engineering :: Physics",
        "Operating System :: OS Independent",
    ),
)
