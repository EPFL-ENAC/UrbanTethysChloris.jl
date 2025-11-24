"""
    WBRoof{FT<:AbstractFloat} <: AbstractModelVariables{FT}

Water balance checks for roof surfaces and soil.

# Fields
- `WBRoofImp`: Water balance check for impervious roof
- `WBRoofVegInVeg`: Water balance check for vegetation interception on roof
- `WBRoofVegInGround`: Water balance check for ground interception under roof vegetation
- `WBRoofVegSoil`: Water balance check for soil under roof vegetation
- `WBRoofVeg`: Water balance check for overall vegetated roof
- `WBRoofTot`: Water balance check for total roof
"""
Base.@kwdef mutable struct WBRoof{FT<:AbstractFloat} <: AbstractModelVariables{FT}
    WBRoofImp::FT
    WBRoofVegInVeg::FT
    WBRoofVegInGround::FT
    WBRoofVegSoil::FT
    WBRoofVeg::FT
    WBRoofTot::FT
end

function WBRoof(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, WBRoof, Dict{String,Any}())
end

"""
    WBCanyonIndv{FT<:AbstractFloat} <: AbstractModelVariables{FT}

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
Base.@kwdef mutable struct WBCanyonIndv{FT<:AbstractFloat} <: AbstractModelVariables{FT}
    WB_In_tree::FT
    WB_In_gveg::FT
    WB_In_gimp::FT
    WB_In_gbare::FT
    WB_Pond_gveg::FT
    WB_Soil_gimp::FT
    WB_Soil_gbare::FT
    WB_Soil_gveg::FT
end

function WBCanyonIndv(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, WBCanyonIndv, Dict{String,Any}())
end

"""
    WBCanyonTot{FT<:AbstractFloat} <: AbstractModelVariables{FT}

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
Base.@kwdef mutable struct WBCanyonTot{FT<:AbstractFloat} <: AbstractModelVariables{FT}
    WBsurf_tree::FT
    WBsurf_imp::FT
    WBsurf_bare::FT
    WBsurf_veg::FT
    WBsoil_imp::FT
    WBsoil_bare::FT
    WBsoil_veg::FT
    WBimp_tot::FT
    WBbare_tot::FT
    WBveg_tot::FT
    WBcanyon_flux::FT
    WBtree_level::FT
    WBground_level::FT
    WBsoil_level::FT
    WBcanyon_level::FT
end

function WBCanyonTot(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, WBCanyonTot, Dict{String,Any}())
end

"""
    EB{FT<:AbstractFloat} <: AbstractModelVariables{FT}

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
Base.@kwdef mutable struct EB{FT<:AbstractFloat} <: AbstractModelVariables{FT}
    EBRoofImp::FT
    EBRoofVeg::FT
    EBGroundImp::FT
    EBGroundBare::FT
    EBGroundVeg::FT
    EBTree::FT
    EBWallSun::FT
    EBWallShade::FT
    EBWallSunInt::FT
    EBWallShadeInt::FT
    EBCanyonT::FT
    EBCanyonQ::FT
end

function EB(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, EB, Dict{String,Any}())
end

"""
    SolverVariables{FT<:AbstractFloat} <: AbstractModelVariables{FT}

Optical properties for indoor building surfaces.

# Fields
- `Success` Boolean indicating convergence of solution of energy balance
- `ValuesEB`: Energy balance closure for the different equations [W/m²]
- `Tsolver`: Temperatures and humidity of different canyon faces and air [K], [kg/kg]
- `YfunctionOutput`: Solver function outputs
"""
Base.@kwdef mutable struct SolverVariables{FT<:AbstractFloat} <: AbstractModelVariables{FT}
    Success::Bool
    ValuesEB::MVector{22,FT}
    Tsolver::MVector{22,FT}
    YfunctionOutput::MVector{22,FT}
end

function SolverVariables(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, SolverVariables, Dict{String,Any}())
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT}, ::Type{SolverVariables}, data::Dict{String,Any}, params::Tuple
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()

    for field in String.(fieldnames(SolverVariables))
        if field != "Success"
            processed[field] = zero(MVector{22,FT})
        else
            processed[field] = false
        end
    end

    return processed
end

"""
    EnergyBalanceVariables{FT<:AbstractFloat} <: AbstractModelVariableSet{FT}

Container for all energy balance variable components.

# Fields
- `WBRoof`: Water balance checks for roof surfaces and soil
- `WBCanyonIndv`: Individual water balance checks for canyon components
- `WBCanyonTot`: Total water balance checks for canyon components and levels
- `EB`: Energy balance checks for different urban components
- `Solver`: Energy balance solver variables
"""
Base.@kwdef struct EnergyBalanceVariables{FT<:AbstractFloat} <: AbstractModelVariableSet{FT}
    WBRoof::WBRoof{FT}
    WBCanyonIndv::WBCanyonIndv{FT}
    WBCanyonTot::WBCanyonTot{FT}
    EB::EB{FT}
    Solver::SolverVariables{FT}
end

function EnergyBalanceVariables(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, EnergyBalanceVariables, Dict{String,Any}())
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT}, ::Type{EnergyBalanceVariables}, data::Dict{String,Any}, params::Tuple
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()
    processed["WBRoof"] = WBRoof(FT)
    processed["WBCanyonIndv"] = WBCanyonIndv(FT)
    processed["WBCanyonTot"] = WBCanyonTot(FT)
    processed["EB"] = EB(FT)
    processed["Solver"] = SolverVariables(FT)
    return processed
end
