#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup


def main():
    setup(
        name="spentfuelgpr",
        version="0.1",
        description="Python part of the MIsoEnrichment module",
        author="Nuclear Verification and Disarmament Group, RWTH Aachen University",
        url="https://github.com/maxschalz/miso_enrichment/",
        license="BSD-3-Clause",
        packages=["spentfuelgpr"],
        classifiers=[
            "License :: OSI Approved :: BSD-3-Clause License",
            "Programming Language :: Python :: 3",
        ],
        install_requires=["numpy", "scipy"],
    )
    return


if __name__ == "__main__":
    main()
