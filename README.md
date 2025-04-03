# UrbanTethysChloris

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://EPFL-ENAC.github.io/UrbanTethysChloris.jl/stable)
[![In development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://EPFL-ENAC.github.io/UrbanTethysChloris.jl/dev)
[![Build Status](https://github.com/EPFL-ENAC/UrbanTethysChloris.jl/workflows/Test/badge.svg)](https://github.com/EPFL-ENAC/UrbanTethysChloris.jl/actions)
[![Test workflow status](https://github.com/EPFL-ENAC/UrbanTethysChloris.jl/actions/workflows/Test.yml/badge.svg?branch=main)](https://github.com/EPFL-ENAC/UrbanTethysChloris.jl/actions/workflows/Test.yml?query=branch%3Amain)
[![Lint workflow Status](https://github.com/EPFL-ENAC/UrbanTethysChloris.jl/actions/workflows/Lint.yml/badge.svg?branch=main)](https://github.com/EPFL-ENAC/UrbanTethysChloris.jl/actions/workflows/Lint.yml?query=branch%3Amain)
[![Docs workflow Status](https://github.com/EPFL-ENAC/UrbanTethysChloris.jl/actions/workflows/Docs.yml/badge.svg?branch=main)](https://github.com/EPFL-ENAC/UrbanTethysChloris.jl/actions/workflows/Docs.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/EPFL-ENAC/UrbanTethysChloris.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/EPFL-ENAC/UrbanTethysChloris.jl)
[![DOI](https://zenodo.org/badge/DOI/FIXME)](https://doi.org/FIXME)
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](CODE_OF_CONDUCT.md)
[![All Contributors](https://img.shields.io/github/all-contributors/EPFL-ENAC/UrbanTethysChloris.jl?labelColor=5e1ec7&color=c0ffee&style=flat-square)](#contributors)
[![BestieTemplate](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/JuliaBesties/BestieTemplate.jl/main/docs/src/assets/badge.json)](https://github.com/JuliaBesties/BestieTemplate.jl)

## How to Cite

If you use UrbanTethysChloris.jl in your work, please cite using the reference given in [CITATION.cff](https://github.com/EPFL-ENAC/UrbanTethysChloris.jl/blob/main/CITATION.cff).

## Contributing

If you want to make contributions of any kind, please first that a look into our [contributing guide directly on GitHub](docs/src/90-contributing.md) or the [contributing page on the website](https://EPFL-ENAC.github.io/UrbanTethysChloris.jl/dev/90-contributing/)

### Creating example input data

To create example input data based on the Zurich example in the [original MATLAB code](https://github.com/NaikaMeili/UTC_ModelCode), run the `data-raw/create_input_data` script. The script will automatically create a `data` subfolder and download a MAT file, `ForcingData_ZH2010.mat`, as well as create an example NetCDF and YAML input for the model. Lastly, smaller files will also be created based on the Zurich example for testing purposes.

First, instantiate the environment and manually add the MAT package:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.add("MAT")'
```

Then, run the script:

```bash
julia --project=. data-raw/create_input_data.jl
```

---

### Contributors

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->
