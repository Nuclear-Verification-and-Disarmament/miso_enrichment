# misoenrichment: a Cyclus multi component isotope enrichment module

`misoenrichment` is a module developed for the nuclear fuel cycle simulator
(Cyclus)[fuelcycle.org]. It provides a Cyclus facility that enriches 
uranium streams composed of two or more isotopes taking into account the 
different enrichment behaviour of minor isotopes such as <sup>234</sup>U (present in
natural as well as in reprocessed uranium) or <sup>236</sup>U (present in 
reprocessed uranium from spent nuclear fuel). The tracking of minor
isotopes makes this module suitable for nuclear archaeology, see, e.g., [3].

Table of Contents
- [Getting Started](#getting-started)
- [Theoretical background](#theoretical-background)
- [References](#references)

## Getting started

## Theoretical background
The implementation of the facility itself and the interaction with Cyclus'
Dynamic Resource Exchange is based on the binary enrichment facility from 
the [Cycamore](https://github.com/cyclus/cycamore) package.

The multi-component isotope calculations are based on mainly on 1. and 2. .
1. derives the mathematical formalism of a matched abundance ratio cascade
using constant overall stage separation factors while 2. gives a new 
physically founded approach to calculating said separation factors.

## References

1. E. von Halle, _Multicomponent Isotope Separation in Matched Abundance 
  Ratio Cascades Composed of Stages With Large Separation Factors_. 
  International Technology Programs (K/ITP--131). Oak Ridge, TN (1987).
2. Houston G. Wood, _Effects of Separation Processes on Minor Uranium 
  Isotopes in Enrichment Cascades_. Science and Global Security, 16:26–36
  (2008), DOI: 10.1080/08929880802361796.
3. Steve Fetter, _Nuclear Archaeology: Verifying Declarations of 
  Fissile-material Production_. Science and Global Security, 3:237-261
  (1993).
