from setuptools import setup

with open("README.md", 'r') as f:
    long_description = f.read()

setup(
  name="evopython",
  packages=["evopython"],
  version="0.0.1",
  license="bsd-2-clause",
  description="Feature resolution from whole-genome alignment data.",
  long_description=long_description,
  long_description_content_type="text/markdown",
  author="Steven Mick",
  author_email="smick@bu.edu",
  url="https://github.com/fiszbein-lab/evopython",
  keywords=["comparative genomics", "evopython"],
  python_requires=">=3.10.8",
  install_requires=[
    "biopython~=1.81"
  ],
  classifiers=[
    "Intended Audience :: Science/Research",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    "License :: OSI Approved :: BSD License",
    "Programming Language :: Python :: 3.10.8",
  ],
)
