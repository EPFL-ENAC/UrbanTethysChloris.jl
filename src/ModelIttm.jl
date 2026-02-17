
Base.@kwdef mutable struct Meteotm1{FT<:AbstractFloat}
    SWRin::SVector{2,FT}
    Rain::SVector{2,FT}
end

"""
    ModelIttm{FT<:AbstractFloat,MR,MG}

Model state at previous timestep (t-1), including both model variables and meteorological forcing.

This composite type stores all variables that need to be maintained from the previous timestep
for use in the current timestep's calculations. It consolidates what were previously individual
`_ittm` variables scattered throughout the simulation loop.

# Type Parameters
- `FT`: Floating-point type (e.g., Float64)
- `MR`: Number of roof soil layers
- `MG`: Number of ground soil layers

# Fields
- `tempvec::TempVec{FT}`: Temperature vector containing surface temperatures [K]
- `humidity::Humidity{FT}`: Canyon and atmospheric humidity variables
- `tempvecb::TempVecB{FT}`: Building interior temperatures [K]
- `interception::Interception{FT}`: Interception on different urban surfaces [mm]
- `exwater::ExWater{FT,MR,MG}`: Extractable water for plants from soil [mm m²/m² ground h]
- `vwater::Vwater{FT,MR,MG}`: Water volume in soil layers [mm]
- `owater::Owater{FT,MR,MG}`: Soil moisture in different soil layers [-]
- `soilpotw::SoilPotW{FT}`: Soil water potential for plants [MPa]
- `cico2leaf::CiCO2Leaf{FT}`: Intercellular CO2 concentration in leaves [umolCO2/mol]
- `runon::Runon{FT}`: Runon variables for urban area [mm/time step]
- `tempdamp::TempDamp{FT}`: Temperature dampening fields [K]
- `qinlat::Qinlat{FT,MG}`: Lateral soil water flux variables [mm/h]
- `resistance::Resistance{FT}`: Aerodynamic and stomatal resistances [s/m]
- `meteo::Meteotm1{FT}`: Meteorological forcing at t-1 and t
"""
Base.@kwdef mutable struct ModelIttm{FT<:AbstractFloat,MR,MG}
    tempvec::ModelComponents.ModelVariables.TempVec{FT}
    humidity::ModelComponents.ModelVariables.Humidity{FT}
    tempvecb::ModelComponents.ModelVariables.TempVecB{FT}
    interception::ModelComponents.ModelVariables.Interception{FT}
    exwater::ModelComponents.ModelVariables.ExWater{FT,MR,MG}
    vwater::ModelComponents.ModelVariables.Vwater{FT,MR,MG}
    owater::ModelComponents.ModelVariables.Owater{FT,MR,MG}
    soilpotw::ModelComponents.ModelVariables.SoilPotW{FT}
    cico2leaf::ModelComponents.ModelVariables.CiCO2Leaf{FT}
    runon::ModelComponents.ModelVariables.Runon{FT}
    tempdamp::ModelComponents.ModelVariables.TempDamp{FT}
    qinlat::ModelComponents.ModelVariables.Qinlat{FT,MG}
    resistance::ModelComponents.ModelVariables.Resistance{FT}
    meteo::Meteotm1{FT}
end

"""
    ModelIttm(model::Model{FT}) where {FT<:AbstractFloat}

Initialize ModelIttm from a Model instance using deepcopy of current state.

# Arguments
- `model::Model{FT}`: The model instance to initialize from

# Returns
- `ModelIttm{FT,MR,MG}`: A new ModelIttm instance with deepcopied values from the model
"""
function ModelIttm(model::Model{FT}) where {FT<:AbstractFloat}
    MR = model.parameters.soil.roof.ms
    MG = model.parameters.soil.ground.ms
    return ModelIttm{FT,MR,MG}(;
        tempvec=deepcopy(model.variables.temperature.tempvec),
        humidity=deepcopy(model.variables.humidity.Humidity),
        tempvecb=deepcopy(model.variables.buildingenergymodel.TempVecB),
        interception=deepcopy(model.variables.waterflux.Interception),
        exwater=deepcopy(model.variables.waterflux.ExWater),
        vwater=deepcopy(model.variables.waterflux.Vwater),
        owater=deepcopy(model.variables.waterflux.Owater),
        soilpotw=deepcopy(model.variables.waterflux.SoilPotW),
        cico2leaf=deepcopy(model.variables.waterflux.CiCO2Leaf),
        runon=deepcopy(model.variables.waterflux.Runon),
        tempdamp=deepcopy(model.variables.temperature.tempdamp),
        qinlat=deepcopy(model.variables.waterflux.Qinlat),
        resistance=deepcopy(model.variables.environmentalconditions.resistance),
        meteo=Meteotm1(model.forcing.meteorological),
    )
end

"""
    update!(ittm::ModelIttm{FT,MR,MG}, model::Model{FT}) where {FT<:AbstractFloat,MR,MG}

Update ModelIttm with current values from the Model.

This function copies the current state of all model variables into the ModelIttm structure,
effectively storing the state at timestep t as the new t-1 state for the next iteration.

# Arguments
- `ittm::ModelIttm{FT,MR,MG}`: The ModelIttm instance to update
- `model::Model{FT}`: The model instance to copy values from
"""
function update!(
    ittm::ModelIttm{FT,MR,MG}, model::Model{FT}
) where {FT<:AbstractFloat,MR,MG}
    update!(ittm.tempvecb, model.variables.buildingenergymodel.TempVecB)
    update!(ittm.tempvec, model.variables.temperature.tempvec)
    update!(ittm.humidity, model.variables.humidity.Humidity)
    update!(ittm.interception, model.variables.waterflux.Interception)
    update!(ittm.exwater, model.variables.waterflux.ExWater)
    update!(ittm.vwater, model.variables.waterflux.Vwater)
    update!(ittm.owater, model.variables.waterflux.Owater)
    update!(ittm.soilpotw, model.variables.waterflux.SoilPotW)
    update!(ittm.cico2leaf, model.variables.waterflux.CiCO2Leaf)
    update!(ittm.tempdamp, model.variables.temperature.tempdamp)
    update!(ittm.runon, model.variables.waterflux.Runon)
    update!(ittm.qinlat, model.variables.waterflux.Qinlat)
    update!(ittm.resistance, model.variables.environmentalconditions.resistance)

    return nothing
end

function Meteotm1(
    x::ModelComponents.ForcingInputs.MeteorologicalInputs{FT,0}
) where {FT<:AbstractFloat}
    SWRin = x.SAB1_in + x.SAB2_in + x.SAD1_in + x.SAD2_in
    return Meteotm1{FT}(;
        SWRin=SVector{2,FT}(SWRin, SWRin), Rain=SVector{2,FT}(x.Rain, x.Rain)
    )
end

function update!(
    y::Meteotm1{FT}, x::ModelComponents.ForcingInputs.MeteorologicalInputs{FT,0}
) where {FT<:AbstractFloat}
    SWRin = x.SAB1_in + x.SAB2_in + x.SAD1_in + x.SAD2_in
    y.SWRin = SVector{2,FT}(y.SWRin[2], SWRin)
    y.Rain = SVector{2,FT}(y.Rain[2], x.Rain)

    return nothing
end
