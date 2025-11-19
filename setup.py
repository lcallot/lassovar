"""
Setup script for lassovar Python package.
"""
from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="lassovar",
    version="0.9.0",
    author="Laurent Callot (Original R package), Python port by AI",
    author_email="l.callot@gmail.com",
    description="Estimation and forecasting of VAR model with the Lasso",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/lcallot/lassovar",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering :: Mathematics",
    ],
    python_requires=">=3.7",
    install_requires=[
        "numpy>=1.19.0",
        "pandas>=1.1.0",
        "scikit-learn>=0.23.0",
        "scipy>=1.5.0",
        "statsmodels>=0.12.0",
        "joblib>=0.16.0",
    ],
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov>=2.10",
        ],
    },
)
