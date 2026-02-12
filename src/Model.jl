abstract type AbstractModelMetadata end
Base.@kwdef mutable struct ModelMetadata <: AbstractModelMetadata
    const ncdfpath::String
    const yamlpath::String
end

abstract type AbstractModel{FT<:AbstractFloat} end
const AM = AbstractModel

abstract type Abstract1PModel{FT<:AbstractFloat} <: AbstractModel{FT} end
Base.@kwdef mutable struct Model{FT<:AbstractFloat} <: AM{FT}
    variables::ModelComponents.ModelVariables.AbstractModelVariableSet{FT}
    const parameters::ModelComponents.Parameters.ParameterSet{FT}
    forcing::ModelComponents.ForcingInputs.ForcingInputSet{FT,0}
    const metadata::AbstractModelMetadata
end

"""
    Model(::Type{FT}) where {FT<:AbstractFloat}

Initialize a TethysChloris model by loading data from a NetCDF file and a YAML file.

# Arguments
- `FT`: The floating-point type to use (e.g., `Float32` or `Float64`).
- `netcdfpath`: The path to the NetCDF file containing model input data.
- `yamlpath`: The path to the YAML file containing model parameters.

# Returns
- A `Model{FT, MR, MG}` instance initialized with the provided data.
"""
function create_model(
    ::Type{FT}, netcdfpath::AbstractString, yamlpath::AbstractString
) where {FT<:AbstractFloat}
    yamldata = YAML.load_file(yamlpath; dicttype=Dict{String,Any});

    model, forcing = NCDataset(netcdfpath) do netcdfdata
        _create_model(FT, netcdfdata, yamldata, abspath(netcdfpath), abspath(yamlpath))
    end

    return model, forcing
end

function _create_model(
    ::Type{FT},
    netcdfdata::NCDataset,
    yamldata::Dict,
    netcdfpath::AbstractString="",
    yamlpath::AbstractString="",
) where {FT<:AbstractFloat}
    parameters = initialize_parameter_set(FT, yamldata)
    forcing = ForcingInputSet(FT, TimeSeries(), netcdfdata, parameters.location)
    variables = ModelVariableSet(FT, parameters.soil)
    metadata = ModelMetadata(abspath(netcdfpath), abspath(yamlpath))

    model = Model{FT}(;
        variables=variables, parameters=parameters, forcing=forcing[1], metadata=metadata
    )

    return model, forcing
end

"""
    initialize!(model::T, args...) where {T<:AbstractModel}

Provide initial conditions for the model.
"""
function initialize!(model::T, args...) where {T<:AbstractModel}
    return nothing
end

function initialize!(
    model::Model{FT}, forcing::ModelComponents.ForcingInputs.ForcingInputSet{FT,1}
) where {FT<:AbstractFloat}

    # Initialize all vegetation and soil temperatures with atmospheric temperature
    initialize!(model.variables.temperature.tempvec, model.forcing.meteorological.Tatm)

    # Initializes all dampening temperature with the nanmean T_atm
    initialize!(model.variables.temperature.tempdamp, forcing.meteorological.Tatm)

    # Initialize humidity canyon specific variable with atmospheric specific humidity
    initialize!(model.variables.humidity.Humidity, model.forcing.meteorological)

    # Initialize the soil moisture at field capacity
    initialize!(
        model.variables.waterflux,
        model.parameters.soil,
        model.parameters.vegetation,
        FT(400),
    )

    # Initialize building temperature variables with atmospheric temperature
    initialize!(
        model.variables.buildingenergymodel.TempVecB,
        model.forcing.meteorological.Tatm,
        model.forcing.meteorological.q_atm,
    )

    return nothing
end

function initialize!(
    x::ModelComponents.ModelVariables.TempVec{FT}, Tatm::FT
) where {FT<:AbstractFloat}
    # set all layers off all fields to mean_Tatm
    # T2m is technically initialize as Results2m.T2m[1] for the first timestep, but this is
    # actually MeteoData.Tatm[1], so we use Tatm here for simplicity
    for field in fieldnames(typeof(x))
        setproperty!(x, field, Tatm)
    end

    return nothing
end

function initialize!(
    x::ModelComponents.ModelVariables.TempVecB{FT}, Tatm::FT, AtmSpecific::FT
) where {FT<:AbstractFloat}
    temperature_fields = [
        :Tceiling, :Tinwallsun, :Tinwallshd, :Twindows, :Tinground, :Tintmass, :Tbin
    ]
    for var in temperature_fields
        setproperty!(x, var, Tatm)
    end

    humidity_fields = [:qbin]
    for var in humidity_fields
        setproperty!(x, var, AtmSpecific)
    end

    return nothing
end

function initialize!(
    x::ModelComponents.ModelVariables.Humidity{FT},
    meteo::ModelComponents.ForcingInputs.MeteorologicalInputs{FT},
) where {FT<:AbstractFloat}
    # q2m is technically initialize as MeteoData.q_atm but this corresponds to
    # AtmSpecific, so we use AtmSpecific here for simplicity
    setproperty!(x, :CanyonSpecific, meteo.q_atm)
    setproperty!(x, :q2m, meteo.q_atm)
    setproperty!(x, :AtmRelative, meteo.AtmRelative)
    setproperty!(x, :AtmSpecific, meteo.AtmSpecific)
    setproperty!(x, :AtmVapourPre, meteo.AtmVapourPre)
    setproperty!(x, :AtmRelativeSat, meteo.AtmRelativeSat)
    setproperty!(x, :AtmSpecificSat, meteo.AtmSpecificSat)
    setproperty!(x, :AtmVapourPreSat, meteo.AtmVapourPreSat)

    return nothing
end

function initialize!(
    x::ModelComponents.ModelVariables.WaterFluxVariables{FT,MR,MG},
    soil::ModelComponents.Parameters.SoilParameters{FT},
    vegetation::ModelComponents.Parameters.VegetationParameters{FT},
    initial_value::FT,
) where {FT<:AbstractFloat,MR,MG}
    soil_values = calculate_soil_values(soil, vegetation)

    # Initialize the water volumes in soil at field capacity
    initialize!(x.Vwater, soil_values)

    # Initialize the soil moisture at field capacity
    initialize!(x.Owater, soil_values)

    initialize!(x.CiCO2Leaf, initial_value)

    return nothing
end

function calculate_soil_values(
    soil::ModelComponents.Parameters.SoilParameters{FT},
    vegetation::ModelComponents.Parameters.VegetationParameters{FT},
) where {FT<:AbstractFloat}
    roof_soil_values = calculate_soil_values(soil.roof, vegetation.roof, vegetation.roof)

    ground_soil_values = calculate_soil_values(
        soil.ground, vegetation.tree, vegetation.ground
    )

    return (; roof=roof_soil_values, ground=ground_soil_values)
end

function calculate_soil_values(
    soil::ModelComponents.Parameters.VegetatedSoilParameters{FT},
    tree::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
    ground::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
) where {FT<:AbstractFloat}
    _, _, _, Osat, Ohy, _, _, _, _, _, O33 = Soil.soil_parameters_total(
        soil.Pcla,
        soil.Psan,
        soil.Porg,
        soil.Kfc,
        soil.Phy,
        soil.SPAR,
        soil.Kbot,
        tree.CASE_ROOT,
        ground.CASE_ROOT,
        tree.ZR95,
        ground.ZR95,
        tree.ZR50,
        ground.ZR50,
        tree.ZRmax,
        ground.ZRmax,
        soil.Zs,
    )

    unique_Osat = unique(Osat)
    @assert length(unique_Osat) == 1 "Osat should be unique after preprocessing"
    unique_Ohy = unique(Ohy)
    @assert length(unique_Ohy) == 1 "Ohy should be unique after preprocessing"
    unique_O33 = unique(O33)
    @assert length(unique_O33) == 1 "O33 should be unique after preprocessing"

    return (;
        Osat=unique_Osat[], Ohy=unique_Ohy[], O33=unique_O33[], ms=soil.ms, dz=diff(soil.Zs)
    )
end

function initialize!(
    x::ModelComponents.ModelVariables.Vwater{FT,MR,MG}, soil_values::NamedTuple
) where {FT<:AbstractFloat,MR,MG}
    roof_init = soil_values.roof.O33 * soil_values.roof.dz
    roof_fields = [:VRoofSoilVeg]
    for field in roof_fields
        setproperty!(x, field, roof_init)
    end

    ground_init = soil_values.ground.O33 * soil_values.ground.dz
    ground_fields = [:VGroundSoilImp, :VGroundSoilBare, :VGroundSoilVeg, :VGroundSoilTot]
    for field in ground_fields
        setproperty!(x, field, ground_init)
    end

    return nothing
end

function initialize!(
    x::ModelComponents.ModelVariables.Owater{FT,MR,MG}, soil_values::NamedTuple
) where {FT<:AbstractFloat,MR,MG}
    x.OwRoofSoilVeg[:] .= soil_values.roof.O33

    ground_fields = [
        :OwGroundSoilImp, :OwGroundSoilBare, :OwGroundSoilVeg, :OwGroundSoilTot
    ]
    for field in ground_fields
        # set all ground soil layers to ground.O33
        setproperty!(x, field, fill(soil_values.ground.O33, length(getproperty(x, field))))
    end
    x.OwGroundSoilImp[1:2] .= FT(NaN)

    return nothing
end

function initialize!(
    x::ModelComponents.ModelVariables.TempDamp{FT}, Tatm::Vector{FT}
) where {FT<:AbstractFloat}
    mean_Tatm = NaNMath.mean(Tatm)
    # set all layers off all fields to mean_Tatm
    for field in fieldnames(typeof(x))
        setproperty!(x, field, mean_Tatm)
    end

    return nothing
end

function initialize!(
    x::ModelComponents.ModelVariables.CiCO2Leaf{FT}, initial_value::FT
) where {FT<:AbstractFloat}
    # set all layers off all fields to initial_value
    for field in fieldnames(typeof(x))
        setproperty!(x, field, initial_value)
    end

    return nothing
end
