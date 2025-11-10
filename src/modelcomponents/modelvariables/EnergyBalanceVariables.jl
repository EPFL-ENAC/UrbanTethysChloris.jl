"""
    WBRoof{FT<:AbstractFloat, N} <: Abstract1PModelVariables{FT,N}

Water balance checks for roof surfaces and soil.

# Fields
- `WBRoofImp`: Water balance check for impervious roof
- `WBRoofVegInVeg`: Water balance check for vegetation interception on roof
- `WBRoofVegInGround`: Water balance check for ground interception under roof vegetation
- `WBRoofVegSoil`: Water balance check for soil under roof vegetation
- `WBRoofVeg`: Water balance check for overall vegetated roof
- `WBRoofTot`: Water balance check for total roof
"""
Base.@kwdef mutable struct WBRoof{FT<:AbstractFloat,N} <: Abstract1PModelVariables{FT,N}
    WBRoofImp::Array{FT,N}
    WBRoofVegInVeg::Array{FT,N}
    WBRoofVegInGround::Array{FT,N}
    WBRoofVegSoil::Array{FT,N}
    WBRoofVeg::Array{FT,N}
    WBRoofTot::Array{FT,N}
end

function initialize_wbroof(::Type{FT}, ::TimeSlice) where {FT<:AbstractFloat}
    return initialize(FT, WBRoof, Dict{String,Any}(), (FT, dimension_value(TimeSlice())))
end

function initialize_wbroof(::Type{FT}, ::TimeSeries, hours::Int) where {FT<:AbstractFloat}
    return initialize(
        FT, WBRoof, Dict{String,Any}(), (FT, dimension_value(TimeSeries())), hours
    )
end

"""
    WBCanyonIndv{FT<:AbstractFloat, N} <: Abstract1PModelVariables{FT,N}

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
Base.@kwdef mutable struct WBCanyonIndv{FT<:AbstractFloat,N} <:
                           Abstract1PModelVariables{FT,N}
    WB_In_tree::Array{FT,N}
    WB_In_gveg::Array{FT,N}
    WB_In_gimp::Array{FT,N}
    WB_In_gbare::Array{FT,N}
    WB_Pond_gveg::Array{FT,N}
    WB_Soil_gimp::Array{FT,N}
    WB_Soil_gbare::Array{FT,N}
    WB_Soil_gveg::Array{FT,N}
end

function initialize_wbcanyon_indv(::Type{FT}, ::TimeSlice) where {FT<:AbstractFloat}
    return initialize(
        FT, WBCanyonIndv, Dict{String,Any}(), (FT, dimension_value(TimeSlice()))
    )
end

function initialize_wbcanyon_indv(
    ::Type{FT}, ::TimeSeries, hours::Int
) where {FT<:AbstractFloat}
    return initialize(
        FT, WBCanyonIndv, Dict{String,Any}(), (FT, dimension_value(TimeSeries())), hours
    )
end

"""
    WBCanyonTot{FT<:AbstractFloat, N} <: Abstract1PModelVariables{FT,N}

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
Base.@kwdef mutable struct WBCanyonTot{FT<:AbstractFloat,N} <:
                           Abstract1PModelVariables{FT,N}
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

function initialize_wbcanyon_tot(::Type{FT}, ::TimeSlice) where {FT<:AbstractFloat}
    return initialize(
        FT, WBCanyonTot, Dict{String,Any}(), (FT, dimension_value(TimeSlice()))
    )
end

function initialize_wbcanyon_tot(
    ::Type{FT}, ::TimeSeries, hours::Int
) where {FT<:AbstractFloat}
    return initialize(
        FT, WBCanyonTot, Dict{String,Any}(), (FT, dimension_value(TimeSeries())), hours
    )
end

"""
    EB{FT<:AbstractFloat, N} <: Abstract1PModelVariables{FT,N}

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
Base.@kwdef mutable struct EB{FT<:AbstractFloat,N} <: Abstract1PModelVariables{FT,N}
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

function initialize_eb(::Type{FT}, ::TimeSlice) where {FT<:AbstractFloat}
    return initialize(FT, EB, Dict{String,Any}(), (FT, dimension_value(TimeSlice())))
end

function initialize_eb(::Type{FT}, ::TimeSeries, hours::Int) where {FT<:AbstractFloat}
    return initialize(
        FT, EB, Dict{String,Any}(), (FT, dimension_value(TimeSeries())), hours
    )
end

"""
    SolverVariables{FT<:AbstractFloat, N, Np1} <: Abstract2PModelVariables{FT,N,Np1}

Optical properties for indoor building surfaces.

# Fields
- `Success` Boolean indicating convergence of solution of energy balance
- `ValuesEB`: Energy balance closure for the different equations [W/m²]
- `Tsolver`: Temperatures and humidity of different canyon faces and air [K], [kg/kg]
- `YfunctionOutput`: Solver function outputs
"""
Base.@kwdef mutable struct SolverVariables{FT<:AbstractFloat,N,Np1} <:
                           Abstract2PModelVariables{FT,N,Np1}
    Success::Array{Bool,N}
    ValuesEB::Array{FT,Np1}
    Tsolver::Array{FT,Np1}
    YfunctionOutput::Array{FT,Np1}
end

function get_vector_fields(obj::SolverVariables)
    return [:Success]
end

function initialize_solver_variables(::Type{FT}, ::TimeSlice) where {FT<:AbstractFloat}
    N = dimension_value(TimeSlice())
    return initialize(FT, SolverVariables, Dict{String,Any}(), (FT, N, N+1))
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT}, ::Type{SolverVariables}, data::Dict{String,Any}, params::Tuple
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()

    dimensions = Dict(
        "Success" => (), "ValuesEB" => (22,), "Tsolver" => (22,), "YfunctionOutput" => (22,)
    )

    for (var, dims) in dimensions
        if var != "Success"
            processed[var] = zeros(FT, dims)
        else
            processed[var] = fill(false, dims)
        end
    end

    return processed
end

function initialize_solver_variables(
    ::Type{FT}, ::TimeSeries, hours::Int
) where {FT<:AbstractFloat}
    N = dimension_value(TimeSeries())
    return initialize(FT, SolverVariables, Dict{String,Any}(), (FT, N, N+1), hours)
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT}, ::Type{SolverVariables}, data::Dict{String,Any}, params::Tuple, hours::Int
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()

    dimensions = Dict(
        "Success" => (hours,),
        "ValuesEB" => (hours, 22),
        "Tsolver" => (hours, 22),
        "YfunctionOutput" => (hours, 22),
    )

    for (var, dims) in dimensions
        if var != "Success"
            processed[var] = zeros(FT, dims)
        else
            processed[var] = fill(false, dims)
        end
    end

    return processed
end

"""
    EnergyBalanceVariables{FT<:AbstractFloat, N} <: Abstract1PModelVariablesSet{FT, N}

Container for all energy balance variable components.

# Fields
- `WBRoof`: Water balance checks for roof surfaces and soil
- `WBCanyonIndv`: Individual water balance checks for canyon components
- `WBCanyonTot`: Total water balance checks for canyon components and levels
- `EB`: Energy balance checks for different urban components
- `Solver`: Energy balance solver variables
"""
Base.@kwdef struct EnergyBalanceVariables{FT<:AbstractFloat,N,Np1} <:
                   Abstract2PModelVariablesSet{FT,N,Np1}
    WBRoof::WBRoof{FT,N}
    WBCanyonIndv::WBCanyonIndv{FT,N}
    WBCanyonTot::WBCanyonTot{FT,N}
    EB::EB{FT,N}
    Solver::SolverVariables{FT,N,Np1}
end

function initialize_energy_balance_variables(
    ::Type{FT}, ::TimeSlice
) where {FT<:AbstractFloat}
    N = dimension_value(TimeSlice())
    return initialize(FT, EnergyBalanceVariables, Dict{String,Any}(), (FT, N, N+1))
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT}, ::Type{EnergyBalanceVariables}, data::Dict{String,Any}, params::Tuple
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()
    processed["WBRoof"] = initialize_wbroof(FT, dimensionality_type(params[2]))
    processed["WBCanyonIndv"] = initialize_wbcanyon_indv(FT, dimensionality_type(params[2]))
    processed["WBCanyonTot"] = initialize_wbcanyon_tot(FT, dimensionality_type(params[2]))
    processed["EB"] = initialize_eb(FT, dimensionality_type(params[2]))
    processed["Solver"] = initialize_solver_variables(FT, dimensionality_type(params[2]))
    return processed
end

function initialize_energy_balance_variables(
    ::Type{FT}, ::TimeSeries, hours::Int
) where {FT<:AbstractFloat}
    N = dimension_value(TimeSeries())
    return initialize(FT, EnergyBalanceVariables, Dict{String,Any}(), (FT, N, N+1), hours)
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT},
    ::Type{EnergyBalanceVariables},
    data::Dict{String,Any},
    params::Tuple,
    hours::Int,
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()
    processed["WBRoof"] = initialize_wbroof(FT, dimensionality_type(params[2]), hours)
    processed["WBCanyonIndv"] = initialize_wbcanyon_indv(
        FT, dimensionality_type(params[2]), hours
    )
    processed["WBCanyonTot"] = initialize_wbcanyon_tot(
        FT, dimensionality_type(params[2]), hours
    )
    processed["EB"] = initialize_eb(FT, dimensionality_type(params[2]), hours)
    processed["Solver"] = initialize_solver_variables(
        FT, dimensionality_type(params[2]), hours
    )
    return processed
end
