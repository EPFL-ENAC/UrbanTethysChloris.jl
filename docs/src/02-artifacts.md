# Artifacts

UrbanTethysChloris.jl is checked against its original MATLAB counterpart, by comparing the output of the Julia functions with that of the MATLAB functions. This ensures that the Julia implementation is consistent with the expected behavior defined in the original MATLAB code. The comparisons are performed via the Julia artifact system and two repositories associated with this package

* [URBES-utc-test-data](https://github.com/EPFL-ENAC/URBES-utc-test-data) creates the translation test data and the associated `.tar.gz`
* [UTC_BEM_ModelCode](https://github.com/NaikaMeili/UTC_BEM_ModelCode) contains the original MATLAB code for the Urban Tethys & Chloris model, including the building energy model.

At present, the creation of the artifact is not automated in the UrbanTethysChloris.jl package. Users must manually create and manage the artifacts as needed by updating the Artifact.toml file with the new sha256 hashes of the `.tar.gz` file. To do so, first download the latest version of the file from the [URBES-utc-test-data](https://github.com/EPFL-ENAC/URBES-utc-test-data) repository. Then, run the following code in the Julia REPL

```julia
using Tar, Inflate, SHA

filename = "path/to/your/file.tar.gz"
println("sha256: ", bytes2hex(open(sha256, filename)))
println("git-tree-sha1: ", Tar.tree_hash(IOBuffer(inflate_gzip(filename))))
```
