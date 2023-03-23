# misoenrichment: a nuclear archaeology module for Cyclus
![GitHub](https://img.shields.io/github/license/maxschalz/miso_enrichment)
![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)]

`misoenrichment` is a module for the nuclear fuel cycle simulator
[Cyclus](http://fuelcycle.org) and is developed at the
[Nuclear Verification and Disarmament Group](https://www.nvd.rwth-aachen.de/)
at RWTH Aachen University.
It currently provides two Cyclus facilities, `MIsoEnrich` and `GprReactor`.

`MIsoEnrich` is an enrichment facility that enriches
streams composed of two or more uranium isotopes taking into account the
different enrichment behaviour of minor isotopes such as <sup>234</sup>U (present in
natural as well as in reprocessed uranium) or <sup>236</sup>U (present in
reprocessed uranium from spent nuclear fuel). The tracking of minor
isotopes makes this module suitable for nuclear archaeology, see, e.g., Ref 3.

`GprReactor` is a Cyclus reactor facility that uses Gaussian Process
Regression (GPR) to calculate the composition of the irradiated fuel depending
on various input parameters. Generally, this implementation works for any
reactor type and any input parameters. However, one needs the appropriate
GPR model (which needs to be generated using training data) and depending
on which input parameters are chosen, the source code of `GprReactor` may
need minor tweaking. Additional information on this issue will be given
in future commits.

## MIsoEnrich
### Installation
This module has some dependencies.
Python dependencies ([`scipy`](https://github.com/scipy/scipy) and
[`numpy`](https://github.com/numpy/numpy) are installed automatically via `pip`,
while the C++ dependencies
([CppOptimizationLibrary](https://github.com/PatWie/CppNumericalSolvers/tree/master),
[Eigen](https://eigen.tuxfamily.org/) and
[JSON for Modern C++](https://github.com/nlohmann/json)) are included as Git
submodules.
These need to be fetched first, as shown below:
```
$ git clone https://github.com/Nuclear-Verification-and-Disarmament/miso_enrichment.git
$ cd miso_enrichment
$ git submodule update --init --recursive  # Download C++ dependencies
$ python3 install.py  # Install misoenrichment module
$ misoenrichment_unit_tests  # Run unit tests (optional).
```

### Getting started
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

Additionally, it should be noted that the facility allows to select between
integer number of stages (default) or decimal number of stages.
The former option will exceed (undershoot) the desired product (tails) assay,
the latter option will match both assays as close as possible.

Also note that when using an integer number of stages, the facility supports
downblending of uranium.
This means that the facility tries to match the desired enrichment
level as precise as possible by first enriching the uranium (to a higher
level, as explained above) and then blending the product with uranium from
the feed. This procedure is only performed if the `use_downblending`
variable is set to `True` in the input file.

### Theoretical background
The implementation of the facility itself and the interaction with Cyclus'
Dynamic Resource Exchange is based on the binary enrichment facility from
the [Cycamore](https://github.com/cyclus/cycamore) package.

The multi-component isotope calculations are based on mainly on Refs 1 and 2.
Ref 1 derives the mathematical formalism of a matched abundance ratio cascade
using constant overall stage separation factors while Ref 2 gives a new
physically founded approach to calculating said separation factors.

## GprReactor
:rotating_light: Please note that this module does not work at the moment.
I do not manage to correctly include the JSON dependency in the build system and
will have to investigate this further or even remove this package entirely.
For updates see [issue #6](https://github.com/Nuclear-Verification-and-Disarmament/miso_enrichment/issues/6).

### Requirements
This facility needs Niels Lohmann's [JSON for Modern C++](https://json.nlohmann.me/)
library. It can be downloaded from his [GitHub repository](https://github.com/nlohmann/json)
or using one of the many package managers, see [here](https://github.com/nlohmann/json#package-managers).
Successfully tested using `conda install -c conda-forge nlohmann_json` and using `CMake`.
When using `CMake`, then the package needs to be installed globally
(i.e., under `/usr/local`), as shown in the following:
```
$ git clone https://github.com/nlohmann/json
$ cd json
$ mkdir build
$cd build
$ cmake ..
$ make
$ sudo make install
```
While it _should_ be possible to install `JSON for Modern C++` locally,
i.e. in `~/.local`, this results in `CMake` not finding `nlohmann/json.hpp`
during the `misoenrichment` installation. I will hopefully manage to
fix this in future versions.

Additionally, one needs `Python3` in combination with the `NumPy` and
`SciPy` packages.

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
