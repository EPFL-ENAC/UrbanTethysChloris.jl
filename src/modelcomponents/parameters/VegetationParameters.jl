"""
    HeightDependentVegetationParameters{FT<:AbstractFloat} <: AbstractParameters{FT}

Parameters for the location-specific (tree, roof, ground) vegetation.

# Fields
- `LAI::FT`: Leaf area index [-]
- `SAI::FT`: Stem area index [-].
- `hc::FT`: Canopy height [m].
- `h_disp::FT`:  [-].
- `d_leaf::FT`:  [-].

- `CASE_ROOT::Int`: Type of root profile [-].
- `ZR95::FT`: Root depth, 95th percentile [mm].
- `ZR50::FT`: Root depth, 50th percentile [mm].
- `ZRmax::FT`: Root depth, maximum [mm]
- `Rrootl::FT`: Root length index [m root m-2 PFT].
- `PsiL50::FT`: Water potential at 50 % of leaf hydraulic conductivity [MPa].
- `PsiX50::FT`: Water potential at 50 % of xylem hydraulic conductivity and limit for water extraction from soil [MPa].

- `FI::FT`: Intrinsec quantum efficency [µmol CO2 µmol-1 photons].
- `Do::FT`: Empirical coefficient that expresses the value of vapor pressure deficit at which f(∆e) = 0.5 [Pa].
- `a1::FT`: Empirical parameter linking net assimilaton `AnC`` to stomatal conductance `g_{s,CO_2}` [-].
- `go::FT`: Minimum/cuticular stomatal conductance [mol CO2 m-2 leaf s-1].
- `CT::Int`: Photosyntesis Typology for Plants, Photosynthetic pathway C3 or C4 [-]
- `DSE::FT`: Activation Energy - Plant Dependent, Activation Energy in Photosynthesis for Rubisco Capacity [kJ/mol].
- `Ha::FT`: Entropy factor - Plant Dependent, Activation energy [kJ / mol K].
- `gmes::FT`: Mesophyll conductance, not used [mol CO2/s m2].
- `rjv::FT`: Scaling factor between Jmax and `V_{c,max}` [µmol equivalent µmol-1 CO2].
- `Kopt::FT`: Optical depth of direct beam per unit plant area (?) [-].
- `Knit::FT`: Canopy nitrogen decay coefficient [-].
- `Vmax::FT`: Maximum Rubisco capacity at 25◦C leaf scale [µmol CO2 m-2 s-1]
- `mSl::FT`:  [-].
- `e_rel::FT`: Relative Efficiency of the photosynthesis apparatus due to Age/Day-length [-].
- `e_relN::FT`: Relative efficiency of the photosynthesis apparatus due to N limitations [-].
- `Psi_sto_00::FT`: Soil water potential at the beginning of stomatal closure [MPa].
- `Psi_sto_50::FT`: Soil water potential at 50 % stomatal closure [MPa].
- `Sl::FT`: Specific leaf area of biomass [m^2 /gC] [-].
- `SPARTREE::Int`: Tree root distribution: 1 = Tree roots can access all water in the soil
(imp, bare, veg) equally; 2 = If the tree crown is smaller than the combined vegetated and
bare fraction, then the trees only transpire from these fractions. Otherwise, they also
transpire from the impervious ground fraction. [-].
"""
Base.@kwdef struct HeightDependentVegetationParameters{FT<:AbstractFloat} <:
                   AbstractHeightDependentParameters{FT}
    # General
    LAI::FT
    SAI::FT
    hc::FT
    h_disp::FT
    d_leaf::FT

    # Tree water uptake
    CASE_ROOT::Int
    ZR95::Vector{FT}
    ZR50::Vector{FT}
    ZRmax::Vector{FT}
    Rrootl::Vector{FT}
    PsiL50::Vector{FT}
    PsiX50::Vector{FT}

    # Photosynthesis and transpiration
    FI::FT
    Do::FT
    a1::FT
    go::FT
    CT::Int
    DSE::FT
    Ha::FT
    gmes::FT
    rjv::FT
    Kopt::FT
    Knit::FT
    Vmax::FT
    mSl::FT
    e_rel::FT
    e_relN::FT
    Psi_sto_00::FT
    Psi_sto_50::FT
    Sl::FT

    SPARTREE::Int = 1
end

function initialize_heightdependent_vegetationparameters(
    ::Type{FT}, data::Dict{String,Any}
) where {FT<:AbstractFloat}
    return initialize(FT, HeightDependentVegetationParameters, data, (FT,))
end

function TethysChlorisCore.get_optional_fields(::Type{HeightDependentVegetationParameters})
    return [:SPARTREE]
end

function TethysChlorisCore.validate_fields(
    ::Type{HeightDependentVegetationParameters}, data::Dict{String,Any}
)
    check_extraneous_fields(
        HeightDependentVegetationParameters,
        data,
        String.(fieldnames(HeightDependentVegetationParameters)),
    )

    if data["LAI"] <= 0.0
        throw(ArgumentError("LAI must be > 0"))
    end
end

"""
    VegetationParameters{FT<:AbstractFloat} <: AbstractParameters{FT}

Container for vegetation parameters for different urban surface components.

# Fields
- `roof::HeightDependentVegetationParameters{FT}`: Vegetation parameters for roof vegetation
- `ground::HeightDependentVegetationParameters{FT}`: Vegetation parameters for ground-level vegetation
- `tree::HeightDependentVegetationParameters{FT}`: Vegetation parameters for trees
"""
Base.@kwdef struct VegetationParameters{FT<:AbstractFloat} <: AbstractParameters{FT}
    roof::HeightDependentVegetationParameters{FT}
    ground::HeightDependentVegetationParameters{FT}
    tree::HeightDependentVegetationParameters{FT}
end

function initialize_vegetationparameters(
    ::Type{FT}, data::Dict{String,Any}
) where {FT<:AbstractFloat}
    return initialize(FT, VegetationParameters, data, (FT,))
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT}, ::Type{VegetationParameters}, data::Dict{String,Any}, params::Tuple
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()

    for (key, value) in data
        processed[key] = initialize(FT, HeightDependentVegetationParameters, value)
    end

    return processed
end
