abstract type AbstractEnvironmentalConditions{FT<:AbstractFloat} <:
              AbstractModelVariables{FT} end

"""
    EnvironmentalConditions{FT<:AbstractFloat, N} <: AbstractModelVariables{FT}

Environmental conditions including wind speed, LAI time series, and resistances.

# Fields
## Wind Fields
- `u_Hcan`: Wind speed at canyon calculation height (hdisp + canyon roughness height) [m/s]
- `u_Zref_und`: Wind speed at undercanopy reference height [m/s]
- `u_ZPerson`: Wind speed at person height [m/s]

## LAI Time Series Fields
- `LAI_R`: LAI of roof vegetation [-]
- `LAI_G`: LAI of ground vegetation [-]
- `LAI_T`: LAI of tree vegetation [-]

## Resistance Fields
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
Base.@kwdef struct EnvironmentalConditions{FT<:AbstractFloat,N} <:
                   AbstractEnvironmentalConditions{FT}
    # Wind fields
    u_Hcan::Array{FT,N}
    u_Zref_und::Array{FT,N}
    u_ZPerson::Array{FT,N}

    # LAI Time Series
    LAI_R::Array{FT,N}
    LAI_G::Array{FT,N}
    LAI_T::Array{FT,N}

    # Resistance fields
    raRooftoAtm::Array{FT,N}
    raCanyontoAtmOrig::Array{FT,N}
    rap_LRoof::Array{FT,N}
    rb_LRoof::Array{FT,N}
    r_soilRoof::Array{FT,N}
    rs_sunRoof::Array{FT,N}
    rs_shdRoof::Array{FT,N}
    raCanyontoAtm::Array{FT,N}
    rap_can::Array{FT,N}
    rap_Htree_In::Array{FT,N}
    rb_HGround::Array{FT,N}
    rb_LGround::Array{FT,N}
    r_soilGroundbare::Array{FT,N}
    r_soilGroundveg::Array{FT,N}
    alp_soilGroundbare::Array{FT,N}
    alp_soilGroundveg::Array{FT,N}
    rs_sunGround::Array{FT,N}
    rs_shdGround::Array{FT,N}
    rs_sunTree::Array{FT,N}
    rs_shdTree::Array{FT,N}
    RES_w1::Array{FT,N}
    RES_w2::Array{FT,N}
    rap_W1_In::Array{FT,N}
    rap_W2_In::Array{FT,N}
    rap_Zp1::Array{FT,N}
end

function Base.getproperty(
    obj::EnvironmentalConditions{FT,0}, field::Symbol
) where {FT<:AbstractFloat}
    return getfield(obj, field)[]
end

function initialize_environmental_conditions(
    ::Type{FT}, N::Int, hours::Int=1
) where {FT<:AbstractFloat}
    return initialize(FT, EnvironmentalConditions, Dict{String,Any}(), (FT, N), hours)
end
