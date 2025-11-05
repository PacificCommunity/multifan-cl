# MULTIFAN-CL

[![License](https://img.shields.io/badge/License-BSD_2--Clause-green.svg)](https://opensource.org/licenses/BSD-2-Clause)
![Static Badge](https://img.shields.io/badge/Language-C%2B%2B-blue)

MULTIFAN-CL is a computer program that implements a statistical, length-based,
age-structured model for use in fisheries stock assessment.

The history of the model can be traced back to the MULTIFAN model (Fournier et
al. 1990) that provided a method of analysing time series of length-frequency
data using statistical theory to provide estimates of von Bertalanffy growth
parameters and the proportions-at-age in the length-frequency data. This served
as a basis for developing MULTIFAN-CL (Fournier et al. 1998), a full population
model that integrates fisheries catch, effort and size composition data. A later
extension of the model (Hampton and Fournier 2001) added integrated analysis of
tagging data. MULTIFAN-CL has provided the basis for undertaking stock
assessments of tuna and billfish in the Western and Central Pacific Ocean.

The [User's Guide](manual) describes the model features and documents the
release history since 2011, when Nick Davies joined the development team in the
development and maintenance of the source code and User's Guide. MULTIFAN-CL
became open-source software in 2025, released under the Simplified BSD License.

## Acknowledgements

Many agencies have provided support for the development of MULTIFAN-CL over the
years, and that support is acknowledged here. They include:

- [Australian Fisheries Management Authority](https://www.afma.gov.au) (AFMA)
- [Government of Canada](https://www.canada.ca)
- [Commonwealth Scientific and Industrial Research Organization](https://www.csiro.au) (CSIRO)
- [Food and Agriculture Organization of the United Nations](https://www.fao.org) (FAO)
- [Pacific Community](https://www.spc.int) (SPC)
- [Government of Taiwan](https://www.taiwan.gov.tw)
- [University of Hawaii Pelagic Fisheries Research Program](https://www.soest.hawaii.edu/pfrp/) (PFRP)
- [US National Marine Fisheries Service](https://www.fisheries.noaa.gov) (NMFS), Honolulu Laboratory
- [Western and Central Pacific Fisheries Commission](https://www.wcpfc.int) (WCPFC)

## References

Fournier, D.A., J. Hampton, and J.R. Sibert. 1998. MULTIFAN-CL: a length-based,
age-structured model for fisheries stock assessment, with application to South
Pacific albacore, *Thunnus alalunga*.\
[Can. J. Fish. Aquat. Sci. 55:2105-2116](https://doi.org/10.1139/f98-100).

Fournier, D.A., J.R. Sibert, J. Majkowski, and J. Hampton. 1990. MULTIFAN a
likelihood-based method for estimating growth parameters and age composition
from multiple length frequency data sets illustrated using data for southern
bluefin tuna (*Thunnus maccoyii*).\
[Can. J. Fish. Aquat. Sci. 47:301-317](https://doi.org/10.1139/f90-032).

Hampton, J. and D.A. Fournier. 2001. A spatially disaggregated, length-based,
age-structured population model of yellowfin tuna (*Thunnus albacares*) in the
western and central Pacific Ocean.\
[Mar. Freshwater Res. 52, 937-963](https://doi.org/10.1071/MF01049).

## Examples

Official stock assessments using MULTIFAN-CL are listed on the [SPC
website](https://fame.spc.int/resources/stockassessmentfiles). Since 2023, the
stock assessments are provided in a format that can be explored and run, for
example:

* [2023 Yellowfin](https://github.com/PacificCommunity/ofp-sam-yft-2023-diagnostic#readme)
* [2023 Bigeye](https://github.com/PacificCommunity/ofp-sam-bet-2023-diagnostic#readme)
* [2024 Albacore](https://github.com/PacificCommunity/ofp-sam-alb-2024-diagnostic#readme)
* [2025 Skipjack](https://github.com/PacificCommunity/ofp-sam-skj-2025-diagnostic#readme)
