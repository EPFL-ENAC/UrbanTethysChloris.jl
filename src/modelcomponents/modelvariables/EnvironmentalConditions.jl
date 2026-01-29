"""
    Wind{FT<:AbstractFloat} <: AbstractModelVariables{FT}

# Fields
- `u_Hcan`: Wind speed at canyon calculation height (hdisp + canyon roughness height) [m/s]
- `u_Zref_und`: Wind speed at undercanopy reference height [m/s]
- `u_ZPerson`: Wind speed at person height [m/s]
"""
Base.@kwdef mutable struct Wind{FT<:AbstractFloat} <: AbstractModelVariables{FT}
    u_Hcan::FT
    u_Zref_und::FT
    u_ZPerson::FT
end

function Wind(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, Wind, Dict{String,Any}())
end

"""
    LAITimeSeries{FT<:AbstractFloat} <: AbstractModelVariables{FT}

# Fields
- `LAI_R`: LAI of roof vegetation [-]
- `LAI_G`: LAI of ground vegetation [-]
- `LAI_T`: LAI of tree vegetation [-]
"""
Base.@kwdef mutable struct LAITimeSeries{FT<:AbstractFloat} <: AbstractModelVariables{FT}
    LAI_R::FT
    LAI_G::FT
    LAI_T::FT
end

function LAITimeSeries(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, LAITimeSeries, Dict{String,Any}())
end

"""
    Resistance{FT<:AbstractFloat} <: AbstractModelVariables{FT}

# Fields
- `raRooftoAtm`: Aerodynamic resistance ra from roof to atmosphere [s/m]
- `raCanyontoAtmOrig`: Original aerodynamic resistance (without enhancement term) from canyon to atmosphere [s/m]
- `rap_LRoof`: Undercanopy resistance rap_L roof [s/m]
- `rb_LRoof`: Leaf boundary resistance rb_L roof [s/m]
- `r_soilRoof`: Soil resistance rb_soil roof [s/m]
- `rs_sunRoof`: Stomata resistance sunlit vegetation rs_sun_roof [s/m]
- `rs_shdRoof`: Stomata resistance shaded vegetation rs_shd_roof [s/m]
- `raCanyontoAtm`: Aerodynamic resistance ra from canyon to atmosphere [s/m]
- `rap_can`: Aerodynamic urban undercanopy resistance from zom_und to the canyon displacement height plus canyon roughness length [s/m]
- `rap_Htree_In`: Vertical aerodynamic resistance applied from tree height to canyon displacement height plus momentum roughness length height [s/m]
- `rb_HGround`: Leaf boundary resistance rb_H ground [s/m]
- `rb_LGround`: Leaf boundary resistance rb_L ground [s/m]
- `r_soilGroundbare`: Soil resistance from bare soil ground [s/m]
- `r_soilGroundveg`: Soil resistance from vegetated ground [s/m]
- `alp_soilGroundbare`: Factor accounting for soil moisture on r_soil of bare ground [-]
- `alp_soilGroundveg`: Factor accounting for soil moisture on r_soil of vegetated ground [-]
- `rs_sunGround`: Stomata resistance sunlit vegetation rs_sun_ground [s/m]
- `rs_shdGround`: Stomata resistance shaded vegetation rs_shd_ground [s/m]
- `rs_sunTree`: Stomata resistance sunlit vegetation rs_sun_tree [s/m]
- `rs_shdTree`: Stomata resistance shaded vegetation rs_shd_ground [s/m]
- `RES_w1`: Horizontal aerodynamic resistance from canyon wall to air for sunlit wall [s/m]
- `RES_w2`: Horizontal aerodynamic resistance from canyon wall to air for shaded wall [s/m]
- `rap_W1_In`: Vertical aerodynamic resistance applied from 2m height to canyon displacement height plus momentum roughness length height for sunlit wall [s/m]
- `rap_W2_In`: Vertical aerodynamic resistance applied from 2m height to canyon displacement height plus momentum roughness length height for shaded wall [s/m]
- `rap_Zp1`: Vertical aerodynamic resistance from the ground to 2m height [s/m]
"""
Base.@kwdef mutable struct Resistance{FT<:AbstractFloat} <: AbstractModelVariables{FT}
    raRooftoAtm::FT
    raCanyontoAtmOrig::FT
    rap_LRoof::FT
    rb_LRoof::FT
    r_soilRoof::FT
    rs_sunRoof::FT
    rs_shdRoof::FT
    raCanyontoAtm::FT
    rap_can::FT
    rap_Htree_In::FT
    rb_HGround::FT
    rb_LGround::FT
    r_soilGroundbare::FT
    r_soilGroundveg::FT
    alp_soilGroundbare::FT
    alp_soilGroundveg::FT
    rs_sunGround::FT
    rs_shdGround::FT
    rs_sunTree::FT
    rs_shdTree::FT
    RES_w1::FT
    RES_w2::FT
    rap_W1_In::FT
    rap_W2_In::FT
    rap_Zp1::FT
end

function Resistance(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, Resistance, Dict{String,Any}())
end

function update!(x::Resistance{FT}, y::Resistance{FT}) where {FT<:AbstractFloat}
    x.raRooftoAtm = y.raRooftoAtm
    x.raCanyontoAtmOrig = y.raCanyontoAtmOrig
    x.rap_LRoof = y.rap_LRoof
    x.rb_LRoof = y.rb_LRoof
    x.r_soilRoof = y.r_soilRoof
    x.rs_sunRoof = y.rs_sunRoof
    x.rs_shdRoof = y.rs_shdRoof
    x.raCanyontoAtm = y.raCanyontoAtm
    x.rap_can = y.rap_can
    x.rap_Htree_In = y.rap_Htree_In
    x.rb_HGround = y.rb_HGround
    x.rb_LGround = y.rb_LGround
    x.r_soilGroundbare = y.r_soilGroundbare
    x.r_soilGroundveg = y.r_soilGroundveg
    x.alp_soilGroundbare = y.alp_soilGroundbare
    x.alp_soilGroundveg = y.alp_soilGroundveg
    x.rs_sunGround = y.rs_sunGround
    x.rs_shdGround = y.rs_shdGround
    x.rs_sunTree = y.rs_sunTree
    x.rs_shdTree = y.rs_shdTree
    x.RES_w1 = y.RES_w1
    x.RES_w2 = y.RES_w2
    x.rap_W1_In = y.rap_W1_In
    x.rap_W2_In = y.rap_W2_In
    x.rap_Zp1 = y.rap_Zp1
end

"""
    EnvironmentalConditionSet{FT<:AbstractFloat} <: AbstractModelVariableSet{FT}

Environmental condition set including wind, LAI time series, and resistances.
"""
Base.@kwdef struct EnvironmentalConditions{FT<:AbstractFloat} <:
                   AbstractModelVariableSet{FT}
    wind::Wind{FT}
    LAI_time_series::LAITimeSeries{FT}
    resistance::Resistance{FT}
end

function EnvironmentalConditions(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, EnvironmentalConditions, Dict{String,Any}())
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT}, ::Type{EnvironmentalConditions}, data::Dict{String,Any}, params::Tuple
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()
    processed["wind"] = Wind(FT)
    processed["LAI_time_series"] = LAITimeSeries(FT)
    processed["resistance"] = Resistance(FT)
    return processed
end

function ModelComponents.outputs_to_save(
    ::Type{EnvironmentalConditions}, ::Type{EssentialOutputs}
)
    return (:wind,)
end

function ModelComponents.outputs_to_save(
    ::Type{EnvironmentalConditions}, ::Type{ExtendedOutputs}
)
    return (:resistance,)
end
