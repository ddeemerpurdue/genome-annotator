'''
Setup file for genome-annotator
'''

import glob
from setuptools import setup, find_packages

setup(
    name="Genotator",
    version="1.0.0",
    author="Dane Deemer",
    author_email="ddeemer@purdue.edu",

    description="Genome Annotator",
    long_description="Add to this",
    long_description_content_type="text/markdown",

    url="https://github.com/ddeemerpurdue/genome-annotator",

    scripts=[script for script in glob.glob('bin/*')],
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    # install_requires=[
    #     'bio',
    #     'matplotlib',
    #     'numpy',
    #     'pandas',
    #     'pysam'
    # ],
    python_requires='==3.10.*',
)