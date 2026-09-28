using NCDatasets: NCDataset, defDim, defGroup, defVar, sync
using TethysChlorisCore:
    VariableSpec,
    output_specs,
    add_variable!,
    chunksizes,
    axis_names,
    StreamTask,
    make_task,
    store_task!,
    flush_tasks!,
    storage_frequency,
    hourly_storage,
    daily_storage,
    static_storage,
    AllOutputs
using ..UrbanTethysChloris: AbstractModel, UrbanTethysChloris
using ..UrbanTethysChloris.ModelComponents: ModelVariableSet, ParameterSet

"""
    OutputManager{FT}

Streamed NetCDF output manager for the (plot-scale) UrbanTethysChloris simulations,
mirroring the `OutputManager1D` of TethysChloris.jl.

The `model` field is deliberately untyped.

The variables are never hardcoded: they are enumerated from the `outputs_to_save` traits
(`output_specs`), merged with caller-provided `extra` specs and passed through a
caller-provided `spec_filter`. Their location is determined by the component
`storage_frequency`:
- `hourly_storage` variables are streamed with the buffered, fixed `time` axis (defined
  under their subsystem group); precomputed tasks extract and buffer the values, no
  reflection or dictionary lookups happen at save time;
- `daily_storage` variables are buffered across days into day-slot buffers shaped
  `(days_buf_len, ...)`, flushed chunk-wise every `days_buf_len` days (the state of the
  last hourly step of each day is kept — per-step stores overwrite the slot);
- `static_storage` variables (forcing and parameters) are written once when the manager is
  created.

# Fields
- `ds::NCDataset`: output dataset
- `model`: the simulated model
- `filename::String`: output path
- `hourly_tasks::Vector{StreamTask}`: precomputed extraction+buffer tasks for the specs
  streamed with the buffered `time` axis
- `daily_tasks::Vector{StreamTask}`: precomputed extraction+buffer tasks for the specs
  written to the `days` axis
- `buffer_len::Int`: number of timesteps buffered between flushes in the `time`
  dimension (aligned chunks)
- `days_buf_len::Int`: number of days buffered between chunk flushes in the `days`
  dimension (aligned chunks, equals the `days` chunk size)
- `tidx::Int`: next row index of the `time` dimension
- `daylast::Int`: last day index written to the `days` axis
- `ft::Type{<:AbstractFloat}`: floating-point type of the simulated model
"""
Base.@kwdef mutable struct OutputManager{FT}
    ds::NCDataset
    model::Any
    filename::String
    hourly_tasks::Vector{StreamTask}
    daily_tasks::Vector{StreamTask}
    buffer_len::Int
    days_buf_len::Int
    tidx::Int
    daylast::Int
    ft::Type{FT}
end

"""
    initialize_outputs(
        model::AbstractModel,
        output_level,
        NN::Int;
        filename::AbstractString,
        forcing = nothing,
        compress = false,
        buffer_len::Int = 256,
        time_chunk::Int = 32,
        extra::Vector{VariableSpec} = VariableSpec[],
        spec_filter = identity,
    ) -> OutputManager

Create the streamed NetCDF output file for an UrbanTethysChloris simulation.

All dimensions are derived from the model components (`rooflayers` ←
`model.parameters.soil.roof.ms`, `groundlayers` ← `model.parameters.soil.ground.ms`) plus
the fixed `time` axis (`max(NN, forcing_steps)` rows) and the unlimited `days` (daily)
axis when daily outputs exist. The root `time` coordinate is written once covering the
fixed axis; the root `days` coordinate grows with a ranged write per day as daily saves
progress. Per-variable axis dimensions are discovered from the extracted shapes with the
model attributes as ordered axis names.

The variables come from the `outputs_to_save` traits of the model variable components
(enumerated per component, since their fields are nested one level deeper than in
TethysChloris.jl), of the parameters (flat components at the `ParameterSet` level, nested
components individually), and of the forcing (falling back to `AllOutputs` when the level
defines no forcing outputs), merged with user-provided `extra` `VariableSpec`s and passed
through `spec_filter`. `forcing` defaults to `model.forcing`; pass the full-timestep
`forcing` returned by `create_model` to store the complete time series once at creation
time.

`buffer_len` timesteps are buffered before flush in the `time` dimension; daily values
are buffered `time_chunk` days and flushed chunk-wise, so the chunked HDF5 slab covers a
full `days` chunk.
"""
function initialize_outputs(
    model::AbstractModel,
    output_level,
    NN::Int;
    filename::AbstractString,
    forcing=nothing,
    compress=false,
    buffer_len::Int=256,
    time_chunk::Int=32,
    extra::Vector{VariableSpec}=VariableSpec[],
    spec_filter=identity,
)
    FT = eltype(model.forcing.meteorological.Tatm)
    full_forcing = isnothing(forcing) ? model.forcing : forcing
    isfile(filename) && throw(ArgumentError("Output file $filename already exists."))

    specs, forcing_steps = _probe_specs!(
        _build_variables(model, full_forcing, output_level, extra, spec_filter),
        model,
        full_forcing,
    )
    partition = _partition_by_frequency(specs, model)

    ds = NCDataset(filename, "c"; maskingvalue=FT(NaN))
    ds.attrib["UrbanTethysChloris_version"] = string(pkgversion(UrbanTethysChloris))
    ds.attrib["created_at"] = string(Dates.now())

    time_len = max(NN, forcing_steps)
    _def_dimensions!(ds, specs, model, _known_axes(model), time_len)
    defVar(ds, "time", Float64, (:time,))
    haskey(ds.dim, "days") && defVar(ds, "days", Float64, (:days,))
    ds["time"][1:time_len] = collect(1.0:time_len)
    vars_by_group = _def_variables!(
        ds, specs, model, full_forcing, FT; compress, time_chunk
    )

    getvar = (group, name) -> vars_by_group[group][name]

    hourly_tasks = Vector{StreamTask}(undef, length(partition.hourly))
    for (idx, spec) in enumerate(partition.hourly)
        hourly_tasks[idx] = make_task(
            getvar,
            _group_path(spec),
            spec.name,
            _spec_instance(spec.source, model),
            spec.extractor,
            buffer_len,
            FT,
        )
    end
    daily_tasks = Vector{StreamTask}(undef, length(partition.daily))
    for (idx, spec) in enumerate(partition.daily)
        daily_tasks[idx] = make_task(
            getvar,
            _group_path(spec),
            spec.name,
            _spec_instance(spec.source, model),
            spec.extractor,
            time_chunk,
            FT,
        )
    end

    return OutputManager{FT}(;
        ds=ds,
        model=model,
        filename=String(filename),
        hourly_tasks=hourly_tasks,
        daily_tasks=daily_tasks,
        buffer_len=buffer_len,
        days_buf_len=time_chunk,
        tidx=1,
        daylast=0,
        ft=FT,
    )
end

"""
    save_outputs!(mgr::OutputManager, model, day)

Save one timestep of outputs, called once per timestep from `run_simulation`.

`hourly_storage` variables extract through precomputed tasks (no reflection, no
dictionary lookups) into the step buffers, flushed chunk-wise at aligned chunks of
`buffer_len` steps in the `time` dimension.
`daily_storage` variables extract into the day-slot buffers, overwritten each step so
the state of the last hourly step of the day wins; every `days_buf_len` days the
previous, full `days` chunk is flushed in a single slab write.
The root `days` coordinate variable is extended with a single-row ranged write before
the daily values of the new day are stored, so that the unlimited `days` axis is at
least `day` rows long.
"""
function save_outputs!(mgr::OutputManager, model, day)
    slot = mod(mgr.tidx - 1, mgr.buffer_len) + 1

    if !isempty(mgr.daily_tasks)
        dslot = mod(day - 1, mgr.days_buf_len) + 1
        if day > mgr.daylast
            if mod(day - 1, mgr.days_buf_len) == 0 && mgr.daylast > 0
                flush_tasks!(mgr.daily_tasks, (day - mgr.days_buf_len):(day - 1))
            end
            first = mgr.daylast + 1
            mgr.ds["days"][first:day] = Float64.(first:day)
            mgr.daylast = day
        end
        for t in mgr.daily_tasks
            x = _unwrap_scalar(t.extract(t.source))
            isnothing(x) && continue
            store_task!(t, x, dslot)
        end
    end

    for t in mgr.hourly_tasks
        x = _unwrap_scalar(t.extract(t.source))
        isnothing(x) && continue
        store_task!(t, x, slot)
    end

    if slot == mgr.buffer_len
        flush_tasks!(mgr.hourly_tasks, (mgr.tidx - mgr.buffer_len + 1):mgr.tidx)
        sync(mgr.ds)
    end

    mgr.tidx += 1
    return nothing
end

"""
    finalize_outputs!(mgr::OutputManager)

Flush the still-open trailing `days` chunk of the daily tasks, flush the remaining
buffered steps and close the output dataset.
"""
function finalize_outputs!(mgr::OutputManager)
    if !isempty(mgr.daily_tasks)
        daily_warm = fld(mgr.daylast - 1, mgr.days_buf_len) * mgr.days_buf_len
        flush_tasks!(mgr.daily_tasks, (daily_warm + 1):mgr.daylast)
    end
    warm = fld(mgr.tidx - 1, mgr.buffer_len) * mgr.buffer_len
    flush_tasks!(mgr.hourly_tasks, (warm + 1):(mgr.tidx - 1))
    sync(mgr.ds)
    close(mgr.ds)
    return nothing
end

################################################################################
# Variable selection
################################################################################

function _spec_instance(source::Symbol, model, full_forcing=model.forcing)
    source === :forcing && return full_forcing
    source === :parameters && return model.parameters
    source in fieldnames(ModelVariableSet) && return getfield(model.variables, source)
    source in fieldnames(ParameterSet) && return getfield(model.parameters, source)
    return model
end

_group_path(spec) = string(spec.source) * "." * string(spec.group)

_unwrap_scalar(v) = v isa AbstractArray{<:Any,0} ? only(v) : v

function _build_variables(
    model, full_forcing, output_level, extra::Vector{VariableSpec}, spec_filter
)
    OL = typeof(output_level)

    specs = VariableSpec[]
    for field in fieldnames(ModelVariableSet)
        append!(specs, output_specs(fieldtype(ModelVariableSet, field), OL; source=field))
    end
    append!(specs, output_specs(ParameterSet, OL; source=:parameters))
    for field in fieldnames(ParameterSet)
        append!(specs, output_specs(fieldtype(ParameterSet, field), OL; source=field))
    end
    forcing_specs = output_specs(typeof(full_forcing), OL; source=:forcing)
    isempty(forcing_specs) &&
        (forcing_specs = output_specs(typeof(full_forcing), AllOutputs; source=:forcing))
    append!(specs, forcing_specs)
    append!(specs, extra)

    return _dedupe(spec_filter(specs))
end

function _dedupe(vars::Vector{VariableSpec})
    specs = VariableSpec[]
    seen = Set{Tuple{Symbol,Symbol}}()
    for spec in vars
        key = (spec.name, spec.group)
        key in seen && continue
        push!(seen, key)
        push!(specs, spec)
    end
    return specs
end

function _spec_freq(spec, model)
    spec.source in fieldnames(ModelVariableSet) &&
        return storage_frequency(fieldtype(ModelVariableSet, spec.source))
    return static_storage
end

function _partition_by_frequency(specs::Vector{VariableSpec}, model)
    hourly = VariableSpec[]
    daily = VariableSpec[]
    for spec in specs
        freq = _spec_freq(spec, model)
        freq === hourly_storage && push!(hourly, spec)
        freq === daily_storage && push!(daily, spec)
    end
    return (hourly=hourly, daily=daily)
end

function _known_axes(model)
    return (
        rooflayers=model.parameters.soil.roof.ms,
        groundlayers=model.parameters.soil.ground.ms,
    )
end

function _probe_specs!(vars::Vector{VariableSpec}, model, full_forcing)
    known = _known_axes(model)
    probed = VariableSpec[]
    forcing_steps = 0
    for spec in vars
        spec, steps = _probe_spec!(spec, model, full_forcing, known)
        isnothing(spec) || push!(probed, spec)
        forcing_steps = max(forcing_steps, steps)
    end
    return probed, forcing_steps
end

function _probe_spec!(spec::VariableSpec, model, full_forcing, known)
    instance = _spec_instance(spec.source, model, full_forcing)
    sample = _unwrap_scalar(spec.extractor(instance))
    sample === nothing && return nothing, 0
    sample isa AbstractString && return nothing, 0

    lead = spec.source === :forcing ? [:time] : Symbol[]
    if sample isa Real
        spec.dims = lead
        return spec, 0
    elseif sample isa AbstractArray
        eltype(sample) <: Real || return nothing, 0
        ndims(sample) > 3 && throw(
            ArgumentError(
                "Cannot map output $(spec.name) ($(ndims(sample)) dimensions) to a NetCDF variable",
            ),
        )
        isempty(sample) && return nothing, 0
        if spec.source === :forcing
            spec.dims = vcat(lead, axis_names(size(sample)[2:end]; known=known))
        else
            spec.dims = axis_names(size(sample); known=known)
        end
        return spec, length(sample)
    else
        throw(
            ArgumentError(
                "Cannot map output $(spec.name) ($(typeof(sample))) to a NetCDF variable"
            ),
        )
    end
end

################################################################################
# NetCDF structure
################################################################################

function _dimsize(d::Symbol, known::NamedTuple)
    (d === :time || d === :days) && return Inf
    d in propertynames(known) && return Int(getfield(known, d))
    s = string(d)
    startswith(s, "axis") && return parse(Int, s[5:end])
    throw(ArgumentError("Unknown dimension $d"))
end

# The `time` dimension is fixed at `time_len` rows: partial writes into an unlimited
# dimension that was extended by another variable are silently truncated to one element
# per row by the NetCDF backend (NetCDF error -57 territory otherwise).
function _def_dimensions!(
    ds, specs::Vector{VariableSpec}, model, known::NamedTuple, time_len::Int
)
    needed = Set{Symbol}()
    for spec in specs
        union!(needed, spec.dims)
    end
    any(s -> _spec_freq(s, model) === hourly_storage, specs) && push!(needed, :time)
    any(s -> _spec_freq(s, model) === daily_storage, specs) && push!(needed, :days)
    for d in needed
        size = d === :time ? time_len : _dimsize(d, known)
        defDim(ds, String(d), size)
    end
    return nothing
end

function _ensure_group!(ds, cache::Dict{String,Any}, path::AbstractString)
    haskey(cache, path) && return cache[path]
    new_group = defGroup(ds, path)
    cache[path] = new_group
    return new_group
end

function _def_variables!(
    ds, specs::Vector{VariableSpec}, model, full_forcing, FT; compress, time_chunk
)
    known = _known_axes(model)
    vars_by_group = Dict{String,Dict{Symbol,Any}}()
    group_cache = Dict{String,Any}()

    for spec in specs
        path = _group_path(spec)
        group = _ensure_group!(ds, group_cache, path)
        dict = get!(() -> Dict{Symbol,Any}(), vars_by_group, path)
        instance = _spec_instance(spec.source, model, full_forcing)
        x = spec.extractor(instance)
        freq = _spec_freq(spec, model)
        dims = if freq === static_storage
            Tuple(spec.dims)
        elseif freq === hourly_storage
            (:time, spec.dims...)
        else
            (:days, spec.dims...)
        end
        vtype = if freq === static_storage && x isa Real
            nothing
        else
            (freq === static_storage ? eltype(x) : FT)
        end

        if freq === static_storage && x isa Real
            dict[spec.name] = defVar(group, spec.name, FT(x))
        else
            var = add_variable!(
                group,
                spec.name,
                vtype,
                dims;
                chunksizes=chunksizes(_dimspecs(dims, known); time_chunk=time_chunk),
                compress=compress,
            )
            if freq === static_storage
                _assign_static_block!(var, x)
            end
            dict[spec.name] = var
        end
    end

    return vars_by_group
end

function _dimspecs(dims::Tuple, known::NamedTuple)
    return [(d, _dimsize(d, known)) for d in dims]
end

function _assign_static_block!(var, x::AbstractArray)
    var[1:size(x, 1), ntuple(_ -> Colon(), ndims(x) - 1)...] = x
    return nothing
end
