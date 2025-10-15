# Subtype for the set of subsets of energy balance variables
abstract type AbstractEnergyBalanceVariables{FT<:AbstractFloat} <:
              AbstractModelVariableSet{FT} end

# Subtype for each subset of energy balance variables
abstract type AbstractEnergyBalanceVariablesSubset{FT<:AbstractFloat,N} <:
              AbstractModelVariables{FT} end

abstract type AbstractWBRoof{FT<:AbstractFloat,N} <:
              AbstractEnergyBalanceVariablesSubset{FT,N} end
abstract type AbstractWBCanyonIndv{FT<:AbstractFloat,N} <:
              AbstractEnergyBalanceVariablesSubset{FT,N} end
abstract type AbstractWBCanyonTot{FT<:AbstractFloat,N} <:
              AbstractEnergyBalanceVariablesSubset{FT,N} end
abstract type AbstractEB{FT<:AbstractFloat,N} <: AbstractEnergyBalanceVariablesSubset{FT,N} end
abstract type AbstractSolverVariables{FT<:AbstractFloat,N,Np} <:
              AbstractEnergyBalanceVariablesSubset{FT,N} end

"""
    WBRoof{FT<:AbstractFloat, N} <: AbstractWBRoof{FT,N}

Water balance checks for roof surfaces and soil.

# Fields
- `WBRoofImp`: Water balance check for impervious roof
- `WBRoofVegInVeg`: Water balance check for vegetation interception on roof
- `WBRoofVegInGround`: Water balance check for ground interception under roof vegetation
- `WBRoofVegSoil`: Water balance check for soil under roof vegetation
- `WBRoofVeg`: Water balance check for overall vegetated roof
- `WBRoofTot`: Water balance check for total roof
"""
Base.@kwdef struct WBRoof{FT<:AbstractFloat,N} <: AbstractWBRoof{FT,N}
    WBRoofImp::Array{FT,N}
    WBRoofVegInVeg::Array{FT,N}
    WBRoofVegInGround::Array{FT,N}
    WBRoofVegSoil::Array{FT,N}
    WBRoofVeg::Array{FT,N}
    WBRoofTot::Array{FT,N}
end

"""
    WBCanyonIndv{FT<:AbstractFloat, N} <: AbstractWBCanyonIndv{FT,N}

Individual water balance checks for canyon components.

# Fields
- `WB_In_tree`: Water balance check for tree interception
- `WB_In_gveg`: Water balance check for ground vegetation interception
- `WB_In_gimp`: Water balance check for impervious ground interception
- `WB_In_gbare`: Water balance check for bare ground interception
- `WB_Pond_gveg`: Water balance check for ponding under ground vegetation
- `WB_Soil_gimp`: Water balance check for impervious ground soil
- `WB_Soil_gbare`: Water balance check for bare ground soil
- `WB_Soil_gveg`: Water balance check for vegetated ground soil
"""
Base.@kwdef struct WBCanyonIndv{FT<:AbstractFloat,N} <: AbstractWBCanyonIndv{FT,N}
    WB_In_tree::Array{FT,N}
    WB_In_gveg::Array{FT,N}
    WB_In_gimp::Array{FT,N}
    WB_In_gbare::Array{FT,N}
    WB_Pond_gveg::Array{FT,N}
    WB_Soil_gimp::Array{FT,N}
    WB_Soil_gbare::Array{FT,N}
    WB_Soil_gveg::Array{FT,N}
end

"""
    WBCanyonTot{FT<:AbstractFloat, N} <: AbstractWBCanyonTot{FT,N}

Total water balance checks for canyon components and levels.

# Fields
- `WBsurf_tree`: Water balance check for tree surface
- `WBsurf_imp`: Water balance check for impervious surface
- `WBsurf_bare`: Water balance check for bare surface
- `WBsurf_veg`: Water balance check for vegetated surface
- `WBsoil_imp`: Water balance check for impervious soil
- `WBsoil_bare`: Water balance check for bare soil
- `WBsoil_veg`: Water balance check for vegetated soil
- `WBimp_tot`: Water balance check for total impervious area
- `WBbare_tot`: Water balance check for total bare area
- `WBveg_tot`: Water balance check for total vegetated area
- `WBcanyon_flux`: Water balance check for canyon flux
- `WBtree_level`: Water balance check at tree level
- `WBground_level`: Water balance check at ground level
- `WBsoil_level`: Water balance check at soil level
- `WBcanyon_level`: Water balance check at canyon level
"""
Base.@kwdef struct WBCanyonTot{FT<:AbstractFloat,N} <: AbstractWBCanyonTot{FT,N}
    WBsurf_tree::Array{FT,N}
    WBsurf_imp::Array{FT,N}
    WBsurf_bare::Array{FT,N}
    WBsurf_veg::Array{FT,N}
    WBsoil_imp::Array{FT,N}
    WBsoil_bare::Array{FT,N}
    WBsoil_veg::Array{FT,N}
    WBimp_tot::Array{FT,N}
    WBbare_tot::Array{FT,N}
    WBveg_tot::Array{FT,N}
    WBcanyon_flux::Array{FT,N}
    WBtree_level::Array{FT,N}
    WBground_level::Array{FT,N}
    WBsoil_level::Array{FT,N}
    WBcanyon_level::Array{FT,N}
end

"""
    EB{FT<:AbstractFloat, N} <: AbstractEB{FT,N}

Energy balance checks for different urban components.

# Fields
- `EBRoofImp`: Energy balance check for impervious roof [W/m² impervious roof area]
- `EBRoofVeg`: Energy balance check for vegetated roof [W/m² vegetated roof area]
- `EBGroundImp`: Energy balance check for impervious ground [W/m² impervious ground area]
- `EBGroundBare`: Energy balance check for bare ground [W/m² bare ground area]
- `EBGroundVeg`: Energy balance check for vegetated ground [W/m² vegetated ground area]
- `EBTree`: Energy balance check for trees [W/m² horizontally projected tree area]
- `EBWallSun`: Energy balance check for sunlit wall [W/m² vertical wall area]
- `EBWallShade`: Energy balance check for shaded wall [W/m² vertical wall area]
- `EBWallSunInt`: Energy balance check for sunlit wall interior [W/m² vertical wall area]
- `EBWallShadeInt`: Energy balance check for shaded wall interior [W/m² vertical wall area]
- `EBCanyonT`: Energy balance check for canyon temperature [W/m² canyon area]
- `EBCanyonQ`: Energy balance check for canyon humidity [kg/kg]
"""
Base.@kwdef struct EB{FT<:AbstractFloat,N} <: AbstractEB{FT,N}
    EBRoofImp::Array{FT,N}
    EBRoofVeg::Array{FT,N}
    EBGroundImp::Array{FT,N}
    EBGroundBare::Array{FT,N}
    EBGroundVeg::Array{FT,N}
    EBTree::Array{FT,N}
    EBWallSun::Array{FT,N}
    EBWallShade::Array{FT,N}
    EBWallSunInt::Array{FT,N}
    EBWallShadeInt::Array{FT,N}
    EBCanyonT::Array{FT,N}
    EBCanyonQ::Array{FT,N}
end

"""
    SolverVariables{FT<:AbstractFloat, N, Np} <: AbstractModelVariables{FT}

Optical properties for indoor building surfaces.

# Fields
- `Success` Boolean indicating convergence of solution of energy balance
- `ValuesEB`: Energy balance closure for the different equations [W/m²]
- `Tsolver`: Temperatures and humidity of different canyon faces and air [K], [kg/kg]
- `YfunctionOutput`: Solver function outputs
"""
Base.@kwdef struct SolverVariables{FT<:AbstractFloat,N,Np} <:
                   AbstractSolverVariables{FT,N,Np}
    Success::Array{Bool,N}
    ValuesEB::Array{FT,Np}
    Tsolver::Array{FT,Np}
    YfunctionOutput::Array{FT,Np}
end

function get_vector_fields(obj::SolverVariables{FT,0,1}) where {FT<:AbstractFloat}
    return [:Success]
end

function Base.getproperty(
    obj::SolverVariables{FT,0,1}, field::Symbol
) where {FT<:AbstractFloat}
    if field in get_vector_fields(obj)
        return getfield(obj, field)[]
    else
        return getfield(obj, field)
    end
end

function initialize_solver_variables(
    ::Type{FT}, N::Int, hours::Int=1
) where {FT<:AbstractFloat}
    return initialize(FT, SolverVariables, Dict{String,Any}(), (FT, N, N+1), hours)
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT}, ::Type{T}, data::Dict{String,Any}, params::Tuple, hours::Int
) where {FT<:AbstractFloat,T<:AbstractSolverVariables}
    processed = Dict{String,Any}()

    dimensions = get_dimensions(T, data, params, hours)

    for (var, dims) in dimensions
        if var != "Success"
            processed[var] = zeros(FT, dims)
        else
            processed[var] = fill(false, dims)
        end
    end

    return processed
end

function get_dimensions(
    ::Type{T}, data::Dict{String,Any}, params::Tuple, hours::Int
) where {T<:AbstractSolverVariables}
    if params[2] ∉ [0, 1]
        throw(ArgumentError("Only N=0 and N=1 are currently supported"))
    end

    if params[2] == 0
        return Dict(
            "Success" => (),
            "ValuesEB" => (22,),
            "Tsolver" => (22,),
            "YfunctionOutput" => (22,),
        )
    else
        return Dict(
            "Success" => (hours,),
            "ValuesEB" => (hours, 22),
            "Tsolver" => (hours, 22),
            "YfunctionOutput" => (hours, 22),
        )
    end
end

"""
    EnergyBalanceVariables{FT<:AbstractFloat, N} <: AbstractEnergyBalanceVariables{FT}

Container for all energy balance variable components.

# Fields
- `WBRoof`: Water balance checks for roof surfaces and soil
- `WBCanyonIndv`: Individual water balance checks for canyon components
- `WBCanyonTot`: Total water balance checks for canyon components and levels
- `EB`: Energy balance checks for different urban components
- `Solver`: Energy balance solver variables
"""
Base.@kwdef struct EnergyBalanceVariables{FT<:AbstractFloat,N,Np} <:
                   AbstractEnergyBalanceVariables{FT}
    WBRoof::WBRoof{FT,N}
    WBCanyonIndv::WBCanyonIndv{FT,N}
    WBCanyonTot::WBCanyonTot{FT,N}
    EB::EB{FT,N}
    Solver::SolverVariables{FT,N,Np}
end

# Base getproperty methods for scalar access
function Base.getproperty(
    obj::T, field::Symbol
) where {FT<:AbstractFloat,T<:AbstractEnergyBalanceVariablesSubset{FT,0}}
    return getfield(obj, field)[]
end

# Initialization functions for individual components
function initialize_wbroof(::Type{FT}, N::Int, hours::Int=1) where {FT<:AbstractFloat}
    return initialize(FT, WBRoof, Dict{String,Any}(), (FT, N), hours)
end

function initialize_wbcanyon_indv(
    ::Type{FT}, N::Int, hours::Int=1
) where {FT<:AbstractFloat}
    return initialize(FT, WBCanyonIndv, Dict{String,Any}(), (FT, N), hours)
end

function initialize_wbcanyon_tot(::Type{FT}, N::Int, hours::Int=1) where {FT<:AbstractFloat}
    return initialize(FT, WBCanyonTot, Dict{String,Any}(), (FT, N), hours)
end

function initialize_eb(::Type{FT}, N::Int, hours::Int=1) where {FT<:AbstractFloat}
    return initialize(FT, EB, Dict{String,Any}(), (FT, N), hours)
end

# Main initialization function for EnergyBalanceVariables
function initialize_energy_balance_variables(
    ::Type{FT}, N::Int, hours::Int=1
) where {FT<:AbstractFloat}
    return initialize(FT, EnergyBalanceVariables, Dict{String,Any}(), (FT, N, N+1), hours)
end

# Main preprocess_fields method for EnergyBalanceVariables
function TethysChlorisCore.preprocess_fields(
    ::Type{FT},
    ::Type{EnergyBalanceVariables},
    data::Dict{String,Any},
    params::Tuple,
    hours::Int,
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()

    # Initialize each component
    processed["WBRoof"] = initialize_wbroof(FT, params[2], hours)
    processed["WBCanyonIndv"] = initialize_wbcanyon_indv(FT, params[2], hours)
    processed["WBCanyonTot"] = initialize_wbcanyon_tot(FT, params[2], hours)
    processed["EB"] = initialize_eb(FT, params[2], hours)
    processed["Solver"] = initialize_solver_variables(FT, params[2], hours)

    return processed
end
