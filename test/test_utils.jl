module TestUtils

using NCDatasets
using YAML
using Dates
using UrbanTethysChloris.ModelComponents.ForcingInputs: MeteorologicalInputs
using UrbanTethysChloris.ModelComponents.Parameters:
    BuildingEnergyModelParameters,
    HVACParameters,
    HeightDependentVegetationParameters,
    IndoorOpticalProperties,
    LocationProperties,
    LocationSpecificSurfaceFractions,
    LocationSpecificThermalProperties,
    SimpleOpticalProperties,
    ThermalBuilding,
    TreeThermalProperties,
    UrbanGeometryParameters,
    VegetatedOpticalProperties,
    VegetatedSoilParameters,
    WallSoilParameters,
    WindowParameters
using MAT
using LazyArtifacts
using Artifacts
using JSON

"""
    convert_special_values(data)

Convert special string values in the input to their appropriate numerical representations.

This function processes the input (typically from JSON data) and converts
representations of special floating point values to their Julia equivalents.

# Arguments
- `data`: Input containing the data to convert

# Returns
- `data`: The processed input with special values converted

# Special conversions:
- Values above 1e308 (`realmax` in MATLAB) → `Inf` (positive infinity)
- `nothing` → `NaN` (not a number)
- Values below -1e308 (`realmin` in MATLAB) → `-Inf` (negative infinity)

# Example
```julia
data = Dict("x" => 1.797693134862316e+308, "y" => nothing, "z" => -1.797693134862316e+308)
result = convert_special_values(data)
# Returns Dict("x" => Inf, "y" => NaN, "z" => -Inf)
```
"""
function convert_special_values(data::Dict)
    for (key, value) in data
        data[key] = convert_special_values(value)
    end
    return data
end

function convert_special_values(data::Number)
    if data > 1e308
        return Inf
    elseif data < -1e308
        return -Inf
    end
    return data
end

function convert_special_values(data::Nothing)
    return NaN
end

function convert_special_values(data::Vector)
    return convert_special_values.(data)
end

function convert_special_values(data)
    return data
end

"""
    convert_to_float_type(x, T::Type)

Recursively convert numerical values in nested structures to the specified float type.

# Arguments
- `x`: Value to convert
- `T::Type`: Target float type (e.g., Float64)

# Returns
- Converted value of appropriate type

# Notes
- Integers and vectors of integers are kept as is
- Other numbers are converted to the target float type
- Arrays are converted element-wise unless they contain integers
- Dictionaries are processed recursively
- Other types are left unchanged
"""
function convert_to_float_type(x::Integer, T::Type)
    return x
end

function convert_to_float_type(x::Vector{<:Integer}, T::Type)
    return x
end

function convert_to_float_type(x::Number, T::Type)
    return convert(T, x)
end

function convert_to_float_type(x::AbstractArray, T::Type)
    eltype(x) <: Integer && return x
    return convert.(T, x)
end

function convert_to_float_type(d::Dict, T::Type)
    return Dict(k => convert_to_float_type(v, T) for (k, v) in d)
end

convert_to_float_type(x, T::Type) = x

"""
    load_matlab_data(path::String)

Load MATLAB data from MAT files containing input and output variables.

This function reads paired input/output MAT files that were generated from MATLAB data.
It converts the file paths to point to corresponding MAT files and reads both input
and output data, converting special values as needed.

# Arguments
- `path::String`: Base path to the MAT files. If path ends in ".jl", it will be converted to ".json"

# Returns
- `Tuple{Any,Any}`: A tuple containing:
  - `inputVars`: Processed input variables from the input MAT file
  - `outputVars`: Processed output variables from the output MAT file

# Notes
- Files are expected to be in "utcdata" artifact
- Numbers are parsed as Float64

"""
function load_matlab_data(func_name::String)
    FT = Float64

    function json_parser(path::String)
        json = open(path, "r") do f
            JSON.parse(f)
        end
        return json
    end

    dir = artifact"utcdata"

    input_path = joinpath(dir, "inputs", func_name)
    output_path = joinpath(dir, "outputs", func_name)

    input_json = json_parser(input_path)
    output_json = json_parser(output_path)

    inputVars = convert_special_values(input_json)
    outputVars = convert_special_values(output_json)

    return convert_to_float_type(inputVars, FT), convert_to_float_type(outputVars, FT)
end

"""
    load_test_netcdf(path::String = joinpath(@__DIR__, "data", "input_data.nc"))

Load NetCDF input data for testing.

# Arguments
- `path`: Path to NetCDF file containing test data, defaults to data/input_data.nc

# Returns
- `NCDataset`: Dataset containing test input data
"""
function load_test_netcdf(path::String=joinpath(@__DIR__, "data", "input_data.nc"))
    return NCDataset(path)
end

"""
    load_test_parameters(path::String = joinpath(@__DIR__, "data", "parameters.yaml"))

Load and preprocess model parameters from YAML file for testing.

# Arguments
- `path`: Path to YAML parameter file, defaults to data/parameters.yaml

# Returns
- `Dict`: Preprocessed parameter dictionary with:
"""
function load_test_parameters(path::String=joinpath(@__DIR__, "data", "parameters.yaml"))
    return YAML.load_file(path; dicttype=Dict{String,Any})
end

function create_building_energy_model_parameters(
    ::Type{FT};
    indoor_optical::IndoorOpticalProperties{FT}=create_indoor_optical_properties(FT),
    thermal::ThermalBuilding{FT}=create_thermal_building(FT),
    windows::WindowParameters{FT}=create_window_parameters(FT),
    hvac::HVACParameters{FT}=create_hvac_parameters(FT),
) where {FT<:AbstractFloat}
    return BuildingEnergyModelParameters{FT}(;
        indoor_optical=indoor_optical, thermal=thermal, windows=windows, hvac=hvac
    )
end

function create_hvac_parameters(
    ::Type{FT};
    ACon::Bool=false,
    AC_onCool::Bool=false,
    AC_onDehum::Bool=false,
    MasterOn::Bool=true,
    Heatingon::Bool=false,
    TsetpointCooling::FT=FT(0.0),
    TsetpointHeating::FT=FT(0.0),
    RHsetpointCooling::FT=FT(0.0),
    RHsetpointHeating::FT=FT(0.0),
    ACH::FT=FT(0.0),
    COPAC::FT=FT(0.0),
    COPHeat::FT=FT(0.0),
    f_ACLatentToQ::FT=FT(0.0),
    q_RHspCooling::FT=FT(0.0),
) where {FT<:AbstractFloat}
    return HVACParameters{FT}(;
        ACon=ACon,
        AC_onCool=AC_onCool,
        AC_onDehum=AC_onDehum,
        MasterOn=MasterOn,
        Heatingon=Heatingon,
        TsetpointCooling=TsetpointCooling,
        TsetpointHeating=TsetpointHeating,
        RHsetpointCooling=RHsetpointCooling,
        RHsetpointHeating=RHsetpointHeating,
        ACH=ACH,
        COPAC=COPAC,
        COPHeat=COPHeat,
        f_ACLatentToQ=f_ACLatentToQ,
        q_RHspCooling=q_RHspCooling,
    )
end

function create_height_dependent_vegetation_parameters(
    ::Type{FT};
    LAI::FT=FT(0.0),
    SAI::FT=FT(0.0),
    hc::FT=FT(0.0),
    h_disp::FT=FT(0.0),
    d_leaf::FT=FT(0.0),
    CASE_ROOT::Int=0,
    ZR95::Vector{FT}=[0.0],
    ZR50::Vector{FT}=[0.0],
    ZRmax::Vector{FT}=[0.0],
    Rrootl::Vector{FT}=[0.0],
    PsiL50::Vector{FT}=[0.0],
    PsiX50::Vector{FT}=[0.0],
    FI::FT=FT(0.0),
    Do::FT=FT(0.0),
    a1::FT=FT(0.0),
    go::FT=FT(0.0),
    CT::Int=0,
    DSE::FT=FT(0.0),
    Ha::FT=FT(0.0),
    gmes::FT=FT(0.0),
    rjv::FT=FT(0.0),
    Kopt::FT=FT(0.0),
    Knit::FT=FT(0.0),
    Vmax::FT=FT(0.0),
    mSl::FT=FT(0.0),
    e_rel::FT=FT(0.0),
    e_relN::FT=FT(0.0),
    Psi_sto_00::FT=FT(0.0),
    Psi_sto_50::FT=FT(0.0),
    Sl::FT=FT(0.0),
    SPARTREE::Int=1,
) where {FT<:AbstractFloat}
    return HeightDependentVegetationParameters{FT}(;
        LAI=LAI,
        SAI=SAI,
        hc=hc,
        h_disp=h_disp,
        d_leaf=d_leaf,
        CASE_ROOT=CASE_ROOT,
        ZR95=ZR95,
        ZR50=ZR50,
        ZRmax=ZRmax,
        Rrootl=Rrootl,
        PsiL50=PsiL50,
        PsiX50=PsiX50,
        FI=FI,
        Do=Do,
        a1=a1,
        go=go,
        CT=CT,
        DSE=DSE,
        Ha=Ha,
        gmes=gmes,
        rjv=rjv,
        Kopt=Kopt,
        Knit=Knit,
        Vmax=Vmax,
        mSl=mSl,
        e_rel=e_rel,
        e_relN=e_relN,
        Psi_sto_00=Psi_sto_00,
        Psi_sto_50=Psi_sto_50,
        Sl=Sl,
        SPARTREE=SPARTREE,
    )
end

function create_indoor_optical_properties(
    ::Type{FT};
    abc::FT=FT(0.0),
    abw::FT=FT(0.0),
    abg::FT=FT(0.0),
    abm::FT=FT(0.0),
    ec::FT=FT(0.0),
    eg::FT=FT(0.0),
    ew::FT=FT(0.0),
    em::FT=FT(0.0),
) where {FT<:AbstractFloat}
    return IndoorOpticalProperties{FT}(;
        abc=abc, abw=abw, abg=abg, abm=abm, ec=ec, eg=eg, ew=ew, em=em
    )
end

function create_location_properties(
    ::Type{FT};
    phi::FT=FT(0.0),
    lambda::FT=FT(0.0),
    theta_canyon::FT=FT(0.0),
    DeltaGMT::FT=FT(0.0),
) where {FT<:AbstractFloat}
    return LocationProperties{FT}(;
        phi=phi, lambda=lambda, theta_canyon=theta_canyon, DeltaGMT=DeltaGMT
    )
end

function create_location_specific_surface_fractions(
    ::Type{FT};
    fveg::FT=FT(0.0),
    fbare::FT=FT(NaN),
    fimp::FT=FT(0.0),
    Per_runoff::FT=FT(0.0),
) where {FT<:AbstractFloat}
    return LocationSpecificSurfaceFractions{FT}(;
        fveg=fveg, fbare=fbare, fimp=fimp, Per_runoff=Per_runoff
    )
end

function create_location_specific_thermal_properties(
    ::Type{FT}; lan_dry::FT=FT(0.0), cv_s::FT=FT(0.0)
) where {FT<:AbstractFloat}
    return LocationSpecificThermalProperties{FT}(; lan_dry=lan_dry, cv_s=cv_s)
end

function create_meteorological_inputs(
    ::Type{FT};
    LWR_in::Vector{FT}=FT[],
    SAB1_in::Vector{FT}=FT[],
    SAB2_in::Vector{FT}=FT[],
    SAD1_in::Vector{FT}=FT[],
    SAD2_in::Vector{FT}=FT[],
    Tatm::Vector{FT}=FT[],
    Uatm::Vector{FT}=FT[],
    Pre::Vector{FT}=FT[],
    Rain::Vector{FT}=FT[],
    rel_hum::Vector{FT}=FT[],
    datetime::Vector{DateTime}=DateTime[],
    esat_Tatm::Vector{FT}=FT[],
    ea::Vector{FT}=FT[],
    q_atm::Vector{FT}=FT[],
    qSat_atm::Vector{FT}=FT[],
    SW_dir::Vector{FT}=FT[],
    SW_diff::Vector{FT}=FT[],
    Zatm::FT=FT(0.0),
    Catm_CO2::FT=FT(0.0),
    Catm_O2::FT=FT(0.0),
    SunDSM_MRT::FT=FT(0.0),
    cp_atm::Vector{FT}=FT[],
    rho_atm::Vector{FT}=FT[],
) where {FT<:AbstractFloat}
    return MeteorologicalInputs{FT}(;
        LWR_in=LWR_in,
        SAB1_in=SAB1_in,
        SAB2_in=SAB2_in,
        SAD1_in=SAD1_in,
        SAD2_in=SAD2_in,
        Tatm=Tatm,
        Uatm=Uatm,
        Pre=Pre,
        Rain=Rain,
        rel_hum=rel_hum,
        datetime=datetime,
        esat_Tatm=esat_Tatm,
        ea=ea,
        q_atm=q_atm,
        qSat_atm=qSat_atm,
        SW_dir=SW_dir,
        SW_diff=SW_diff,
        Zatm=Zatm,
        Catm_CO2=Catm_CO2,
        Catm_O2=Catm_O2,
        SunDSM_MRT=SunDSM_MRT,
        cp_atm=cp_atm,
        rho_atm=rho_atm,
    )
end

function create_person_parameters(
    ::Type{FT};
    PositionPx::FT=FT(0.0),
    PositionPz::FT=FT(0.0),
    PersonWidth::FT=FT(0.0),
    PersonHeight::FT=FT(0.0),
    HeightWind::FT=FT(0.0),
) where {FT<:AbstractFloat}
    return PersonParameters{FT}(;
        PositionPx=PositionPx,
        PositionPz=PositionPz,
        PersonWidth=PersonWidth,
        PersonHeight=PersonHeight,
        HeightWind=HeightWind,
    )
end

function create_simple_optical_properties(
    ::Type{FT}; albedo::FT=FT(0.0), emissivity::FT=FT(0.0)
) where {FT<:AbstractFloat}
    return SimpleOpticalProperties{FT}(; albedo=albedo, emissivity=emissivity)
end

function create_thermal_building(
    ::Type{FT};
    IntMassOn::Bool=false,
    FloorHeight::FT=FT(0.0),
    dzFloor::FT=FT(0.0),
    dzWall::FT=FT(0.0),
    lan_ground_floor::FT=FT(0.0),
    cv_ground_floor::FT=FT(0.0),
    lan_floor_IntMass::FT=FT(0.0),
    cv_floor_IntMass::FT=FT(0.0),
    lan_wall_IntMass::FT=FT(0.0),
    cv_wall_IntMass::FT=FT(0.0),
) where {FT<:AbstractFloat}
    return ThermalBuilding{FT}(;
        IntMassOn=IntMassOn,
        FloorHeight=FloorHeight,
        dzFloor=dzFloor,
        dzWall=dzWall,
        lan_ground_floor=lan_ground_floor,
        cv_ground_floor=cv_ground_floor,
        lan_floor_IntMass=lan_floor_IntMass,
        cv_floor_IntMass=cv_floor_IntMass,
        lan_wall_IntMass=lan_wall_IntMass,
        cv_wall_IntMass=cv_wall_IntMass,
    )
end

function create_tree_thermal_properties(
    ::Type{FT}; Cthermal_leaf::FT=FT(0.0)
) where {FT<:AbstractFloat}
    return TreeThermalProperties{FT}(; Cthermal_leaf=Cthermal_leaf)
end

function create_urban_geometry_parameters(
    ::Type{FT};
    Height_canyon::FT=FT(0.0),
    Width_canyon::FT=FT(0.0),
    Width_roof::FT=FT(0.0),
    Height_tree::FT=FT(0.0),
    Radius_tree::FT=FT(0.0),
    Distance_tree::FT=FT(0.0),
    Hcan_max::FT=FT(0.0),
    Hcan_std::FT=FT(0.0),
    trees::Bool=false,
    ftree::FT=FT(0.0),
    hcanyon::FT=FT(0.0),
    wcanyon::FT=FT(0.0),
    wroof::FT=FT(0.0),
    htree::FT=FT(0.0),
    radius_tree::FT=FT(0.0),
    distance_tree::FT=FT(0.0),
    ratio::FT=FT(0.0),
    wcanyon_norm::FT=FT(0.0),
    wroof_norm::FT=FT(0.0),
) where {FT<:AbstractFloat}
    return UrbanGeometryParameters{FT}(;
        Height_canyon=Height_canyon,
        Width_canyon=Width_canyon,
        Width_roof=Width_roof,
        Height_tree=Height_tree,
        Radius_tree=Radius_tree,
        Distance_tree=Distance_tree,
        Hcan_max=Hcan_max,
        Hcan_std=Hcan_std,
        trees=trees,
        ftree=ftree,
        hcanyon=hcanyon,
        wcanyon=wcanyon,
        wroof=wroof,
        htree=htree,
        radius_tree=radius_tree,
        distance_tree=distance_tree,
        ratio=ratio,
        wcanyon_norm=wcanyon_norm,
        wroof_norm=wroof_norm,
    )
end

function create_vegetated_optical_properties(
    ::Type{FT};
    aveg::FT=FT(0.0),
    aimp::FT=FT(0.0),
    abare::FT=FT(NaN),
    albedo::FT=FT(0.0),
    eveg::FT=FT(0.0),
    eimp::FT=FT(0.0),
    ebare::FT=FT(NaN),
    emissivity::FT=FT(0.0),
) where {FT<:AbstractFloat}
    return VegetatedOpticalProperties{FT}(;
        aveg=aveg,
        aimp=aimp,
        abare=abare,
        albedo=albedo,
        eveg=eveg,
        eimp=eimp,
        ebare=ebare,
        emissivity=emissivity,
    )
end

function create_vegetated_soil_parameters(
    ::Type{FT};
    Pcla::FT=FT(0.0),
    Psan::FT=FT(0.0),
    Porg::FT=FT(0.0),
    In_max_imp::FT=FT(0.0),
    In_max_ground::FT=FT(NaN),
    In_max_underveg::FT=FT(NaN),
    In_max_bare::FT=FT(NaN),
    Sp_In::FT=FT(0.0),
    Kimp::FT=FT(0.0),
    Kfc::FT=FT(0.0),
    Phy::FT=FT(0.0),
    SPAR::Int=0,
    Kbot::FT=FT(0.0),
    dz1::FT=FT(NaN),
    dz2::FT=FT(NaN),
    Zs::Vector{FT}=FT[],
    ms::Int=0,
    FixSM::Bool=false,
    FixSM_LayerStart::Int=0,
    FixSM_LayerEnd::Int=0,
) where {FT<:AbstractFloat}
    return VegetatedSoilParameters{FT}(;
        Pcla=Pcla,
        Psan=Psan,
        Porg=Porg,
        In_max_imp=In_max_imp,
        In_max_ground=In_max_ground,
        In_max_underveg=In_max_underveg,
        In_max_bare=In_max_bare,
        Sp_In=Sp_In,
        Kimp=Kimp,
        Kfc=Kfc,
        Phy=Phy,
        SPAR=SPAR,
        Kbot=Kbot,
        dz1=dz1,
        dz2=dz2,
        Zs=Zs,
        ms=ms,
        FixSM=FixSM,
        FixSM_LayerStart=FixSM_LayerStart,
        FixSM_LayerEnd=FixSM_LayerEnd,
    )
end

function create_wall_soil_parameters(
    ::Type{FT}; dz1::FT=FT(0.0), dz2::FT=FT(0.0)
) where {FT<:AbstractFloat}
    return WallSoilParameters{FT}(; dz1=dz1, dz2=dz2)
end

function create_window_parameters(
    ::Type{FT};
    WindowsOn::Int=0,
    GlazingRatio::FT=FT(0.0),
    Uvalue::FT=FT(0.0),
    lan_windows::FT=FT(0.0),
    cv_glass::FT=FT(0.0),
    dztot::FT=FT(0.0),
    SHGC::FT=FT(0.0),
    SolarTransmittance::FT=FT(0.0),
    SolarAbsorptivity::FT=FT(0.0),
    SolarAlbedo::FT=FT(0.0),
) where {FT<:AbstractFloat}
    return WindowParameters{FT}(;
        WindowsOn=WindowsOn,
        GlazingRatio=GlazingRatio,
        Uvalue=Uvalue,
        lan_windows=lan_windows,
        cv_glass=cv_glass,
        dztot=dztot,
        SHGC=SHGC,
        SolarTransmittance=SolarTransmittance,
        SolarAbsorptivity=SolarAbsorptivity,
        SolarAlbedo=SolarAlbedo,
    )
end

end
