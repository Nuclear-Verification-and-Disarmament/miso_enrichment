# misoenrichment: a Cyclus multi component isotope enrichment module

`misoenrichment` is a module developed at the [Nuclear Verification and Disarmament Group](https://www.nvd.rwth-aachen.de/) at RWTH Aachen University for the nuclear fuel cycle simulator
[Cyclus](fuelcycle.org). It provides a Cyclus facility that enriches 
uranium streams composed of two or more isotopes taking into account the 
different enrichment behaviour of minor isotopes such as <sup>234</sup>U (present in
natural as well as in reprocessed uranium) or <sup>236</sup>U (present in 
reprocessed uranium from spent nuclear fuel). The tracking of minor
isotopes makes this module suitable for nuclear archaeology, see, e.g., Ref 3.

Table of Contents
- [Getting Started](#getting-started)
- [Theoretical background](#theoretical-background)
- [References](#references)

## Getting started
An example input file is found in `input/main.py` featuring a
`cycamore::Source` source agent, a `MIsoEnrich` enrichment facility and two
`cycamore::Sink` agents, one for enriched and one for depleted uranium. 
Note that the sink requests a binary composition of enriched uranium (90% 
<sup>235</sup>U, 10% <sup>238</sup> U) and that the enrichment facility
enriches to a level at least equal to the requested one while keeping track
 of the minor isotopes. This implies that one does _not_ need to know the 
final composition of enriched uranium beforehand (its desired assay is 
sufficient). In fact, one cannot request a material with a certain 
concentration in minor isotopes or with a constraint on the minor isotopes
concentration (e.g., to make it ASTM compliant, see Ref 4).

## Theoretical background
The implementation of the facility itself and the interaction with Cyclus'
Dynamic Resource Exchange is based on the binary enrichment facility from 
the [Cycamore](https://github.com/cyclus/cycamore) package.

The multi-component isotope calculations are based on mainly on Refs 1 and 2.
Ref 1 derives the mathematical formalism of a matched abundance ratio cascade
using constant overall stage separation factors while Ref 2 gives a new 
physically founded approach to calculating said separation factors.

## References

1. E. von Halle, _Multicomponent Isotope Separation in Matched Abundance 
  Ratio Cascades Composed of Stages With Large Separation Factors_. 
  International Technology Programs (K/ITP--131). Oak Ridge, TN (1987).
2. Houston G. Wood, _Effects of Separation Processes on Minor Uranium 
  Isotopes in Enrichment Cascades_. Science and Global Security, 16:26–36
  (2008), DOI: [10.1080/08929880802361796](https://doi.org/10.1080/08929880802361796).
3. Steve Fetter, _Nuclear Archaeology: Verifying Declarations of 
  Fissile-material Production_. Science and Global Security, 3:237-261
  (1993).
4. ASTM International. _C787-20 Standard Specification for Uranium 
  Hexafluoride for Enrichment_. West Conshohocken, PA; ASTM International, 2020. 
  DOI: [10.1520/C0787-20](https://doi.org/10.1520/C0787-20).
