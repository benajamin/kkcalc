import setuptools

with open("README.rst", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="kkcalc",
    version="0.7.5",
    author="Benjamin Watts",
    author_email="benjamin.watts@gmail.com",
    description="A calculator for the Kramers Kronig transform of X-ray spectra",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url="https://bitbucket.org/benajamin/kkcalc",
    packages=setuptools.find_packages(),
    package_data={'kkcalc': ['*.json']},
    install_requires=['numpy', 'scipy'],
    classifiers=(
        "Programming Language :: Python :: 3.6",
        "License :: OSI Approved :: zlib/libpng License  ",
        "Topic :: Scientific/Engineering :: Physics",
        "Operating System :: OS Independent",
    ),
)