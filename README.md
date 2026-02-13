# UrbanTethysChloris

<!-- [![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://EPFL-ENAC.github.io/UrbanTethysChloris.jl/stable) -->
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

## Example

8 examples across the globe are provided alongside this package, in `/examples`. Before running the scripts, the Julia environments for the examples must be instantiated with

```bash
julia --project=examples -e 'using Pkg; Pkg.instantiate();'
```

For each example, two scripts are provided: `create_<location_name>_data.jl` and `<location_name>_example.jl`.

The first script will create the NetCDF and YAML files containing the parameters and forcing inputs necessary to run the simulation.
The first script of a given location can be run from the command line as

```bash
julia --project=examples examples/<location_name>/create_<location_name>_data.jl
```

Subsequent scripts can be run in a similar fashion, by changing the script name.

The second script shows how to run the simulation, given the parameters and forcing inputs, as well as visualizing the results. This last script is best run interactively in the Julia REPL to interact with the results.

---

### Contributors

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<table>
  <tbody>
    <tr>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/sphamba"><img src="https://avatars.githubusercontent.com/u/17217484?v=4?s=100" width="100px;" alt="Son Pham-Ba"/><br /><sub><b>Son Pham-Ba</b></sub></a><br /><a href="#projectManagement-sphamba" title="Project Management">📆</a> <a href="#review-sphamba" title="Reviewed Pull Requests">👀</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/hsolleder"><img src="https://avatars.githubusercontent.com/u/9566930?v=4?s=100" width="100px;" alt="H Solleder"/><br /><sub><b>H Solleder</b></sub></a><br /><a href="#code-hsolleder" title="Code">💻</a> <a href="#doc-hsolleder" title="Documentation">📖</a> <a href="#review-hsolleder" title="Reviewed Pull Requests">👀</a></td>
    </tr>
  </tbody>
</table>

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->
