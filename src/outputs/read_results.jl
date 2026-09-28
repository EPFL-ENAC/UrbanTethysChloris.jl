using NCDatasets: NCDataset, dimnames

"""
    unpack_results(ncpath, steps, groups...) -> Dict{Symbol,Dict{Symbol,Any}}

Read the given `groups` (dot-joined `source.group` names, e.g.
`"temperature.tempvec"`) of a streamed results file, materializing every variable as a
plain array with missing values replaced by `NaN`. Variables with a `time` dimension are
sliced to `steps`; static variables are returned as saved. The outer dictionary is keyed
by the last component of each group name, mirroring the keys of the former in-RAM
results dictionary.
"""
function unpack_results(
    ncpath::AbstractString, steps::UnitRange{Int}, groups::AbstractString...
)
    ds = NCDataset(ncpath)
    results = Dict{Symbol,Dict{Symbol,Any}}()
    for path in groups
        group = ds.group[path]
        out = Dict{Symbol,Any}()
        for name in keys(group)
            var = group[name]
            arr = Array(var)
            if "time" in dimnames(var)
                arr = arr[steps, ntuple(_ -> Colon(), ndims(arr) - 1)...]
            elseif ndims(arr) == 0
                out[Symbol(name)] = _missing_scalar(arr)
                continue
            end
            out[Symbol(name)] = replace(arr, missing => NaN)
        end
        key = Symbol(last(split(path, ".")))
        haskey(results, key) && throw(
            ArgumentError(
                "Duplicate results group key $key for $path; already read from another group",
            ),
        )
        results[key] = out
    end
    close(ds)
    return results
end

_missing_scalar(v) = (r=only(v); ismissing(r) ? NaN : r)
