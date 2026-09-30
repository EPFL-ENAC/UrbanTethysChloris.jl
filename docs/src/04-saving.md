# [Saving the simulation results](@id saving)

Running a simulation with `output_filename` streams all selected inputs and outputs to a
single NetCDF file during the simulation:

```julia
model, forcing = create_model(Float64, "path/to/forcing.nc", "path/to/parameters.yaml")
initialize!(model, forcing)
manager, view_factor, view_factor_point = run_simulation(
    model,
    forcing;
    NN = 720,
    output_level = extended_outputs,
    output_filename = "path/to/results.nc",
    compress = false,
    buffer_len = 256,
)
```

`run_simulation` returns a `(manager, ViewFactor, ViewFactorPoint)` triple; `manager` is
`nothing` when `output_filename` is not set, and the result file is closed automatically
when the simulation finishes. The results can be read back per group with
[`unpack_results`](@ref):

```julia
results = unpack_results("path/to/results.nc", 1:NN, "temperature.tempvec", "humidity.Results2m")
T2m = results[:Results2m][:T2m]
```

## Specifying variables to save

The saved variables are never hardcoded: they are enumerated from the `outputs_to_save`
traits of the model components at the given `output_level`. To restrict or extend the
selection for one run, pass `spec_filter` (a function mapping the
`Vector{VariableSpec}` list to the desired subset) and/or `extra` (additional
`TethysChlorisCore.VariableSpec`s):

```julia
using TethysChlorisCore: VariableSpec, output_specs

# only moisture-related variables
spec_filter = specs -> filter(s -> occursin("waterflux", string(s.group)), specs)

# one extra diagnostic
extra = [VariableSpec(:MyAux, :temperature, :extra, Symbol[], :mean, x -> x.variables.temperature.tempvec.T2m)]
```

## File layout

The file contains flat subsystem groups with dot-joined names plus the root `time` axis:

```text
results.nc
├── dims: time (= max(NN, forcing steps), fixed), rooflayers, groundlayers, ...
├── time                       (root coordinate variable)
├── /temperature.tempvec       (hourly: buffered, chunked along the fixed time axis)
├── /waterflux.Owater          (hourly: one row per simulated step)
├── /soil.roof                 (static: parameter leaves, written once)
├── /forcing.meteorological    (static: full forcing time series written at creation)
└── /initial.OwaterInitial     (static: initial soil moisture, written at creation)
```

* `time` is a fixed axis with one row per simulation step (or forcing step, whichever is
  longer), so partial writes are always safe,
* model variables are streamed at `hourly_storage` (one row per timestep, buffered for at
  most `buffer_len` timesteps and flushed chunk-wise),
* `static_storage` components (forcing and parameters) are written once when the manager
  is created,
* an unlimited `days` axis appears only when `daily_storage` outputs are defined.

See the output manager docstrings for the function signatures
(`initialize_outputs`, `save_outputs!`, `finalize_outputs!`).
