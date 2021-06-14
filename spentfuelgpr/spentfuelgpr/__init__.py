#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculate the composition of spent fuel based on a set of trained GPRs.

This module is part of the misoenrichment module for Cyclus, see
https://github.com/maxschalz/miso_enrichment/ .

The module is strongly based on Antonio Figueroa's work, who kindly
provided the original code. For more information, please see
https://doi.org/10.1016/j.anucene.2020.108085
or https://arxiv.org/abs/2006.12921
or https://github.com/FigueroaAC/GPs-for-SpentFuel
"""

__author__ = "Nuclear Verification and Disarmament Group, RWTH Aachen University"
__copyright__ = "Copyright 2020-2021, Nuclear Verification and Disarmament Group, RWTH Aachen University"
__credits__ = ["Antonio Figueroa", "Max Schalz"]
__license__ = "BSD-3-Clause"
__version__ = "2.0"
__maintainer__ = "Max Schalz"

from .kernel import *
from .predict_posterior import *
