from setuptools import setup, find_packages, Extension
import os
import numpy

# Define the extension module
fractional_bs_model = Extension(
    'fractional_bs_model',
    sources=['fractional_bs_model.cpp', 'fractional_bs_model_wrap.cxx'],
    include_dirs=[
        numpy.get_include(), 
        '/usr/include/python3.10'
    ],
    extra_compile_args=['-std=c++11', '-fPIC'],
)

# Read the README.md file for the long description
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="fbs-model-prototype",
    version="1.0.0",
    author="Joel Saucedo",
    author_email="joelasaucedo@proton.me",
    description="Fractional Black-Scholes PDE Solver with calibration mechanisms for cryptocurrency options using C++ and Python bindings via SWIG.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/joel-saucedo/fbs-model-prototype",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: C++",
    ],
    packages=find_packages(),
    ext_modules=[fractional_bs_model],
    python_requires=">=3.10",
    install_requires=[
        "numpy",
        "pandas",
        "scipy"
    ],
    setup_requires=["numpy"],
)
