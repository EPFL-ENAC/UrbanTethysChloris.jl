using NCDatasets
using Dates
using YAML
using MAT

FT = Float64

# Check if data directory exists, if not create it
data_dir = joinpath(@__DIR__, "..", "data")
!isdir(data_dir) && mkdir(data_dir)

# https://github.com/NaikaMeili/UTC_BEM_ModelCode/raw/61af9eeeca7c0fbe6ae19a8d78f4f481c45826aa/UTC_Model/+data_functions/TMYNewDelhi_RadPart.mat

# Define files and their GitHub URLs
repo_url = "https://github.com/NaikaMeili/UTC_BEM_ModelCode/raw/61af9eeeca7c0fbe6ae19a8d78f4f481c45826aa/UTC_Model"
files = Dict(
    "TMYNewDelhi_RadPart.mat" => repo_url * "/+data_functions/TMYNewDelhi_RadPart.mat",
    "ViewFactor_NDLCZ3.mat" => repo_url * "/+data_functions/ViewFactor_ND_LCZ3.mat",
)

# Check each file and download if missing
for (file, url) in files
    filepath = joinpath(data_dir, file)
    if !isfile(filepath)
        @info "Downloading $file..."
        download(url, filepath)
    end
end

# Create the parameters file
data = Dict{String,Any}()

data["urbangeometry"] = Dict{String,Any}(
    "Height_canyon" => 6.5,
    "Width_canyon" => 5.8,
    "Width_roof" => 7.1,
    "Hcan_max" => NaN,
    "Hcan_std" => NaN,
    "trees" => true,
    "ftree" => 1.0,
)

data["urbangeometry"]["Radius_tree"] = data["urbangeometry"]["Width_canyon"] * 0.2 / 4
data["urbangeometry"]["Height_tree"] = 5.0 - data["urbangeometry"]["Radius_tree"]

if data["urbangeometry"]["Radius_tree"] + 2 < data["urbangeometry"]["Width_canyon"]
    data["urbangeometry"]["Distance_tree"] = 1 + data["urbangeometry"]["Radius_tree"]
else
    dt =
        (data["urbangeometry"]["Width_canyon"] - (4*data["urbangeometry"]["Radius_tree"])) /
        3
    data["urbangeometry"]["Distance_tree"] = dt + data["urbangeometry"]["Radius_tree"]
end

data["surfacefractions"] = Dict{String,Any}()
data["surfacefractions"]["roof"] = Dict{String,Any}(
    "fveg" => 0.0, "fimp" => 1.0, "Per_runoff" => 1.0
)
data["surfacefractions"]["ground"] = Dict{String,Any}(
    "fveg" => 0.23, "fbare" => 0.1, "Per_runoff" => 0.5
)
data["surfacefractions"]["ground"]["fimp"] =
    1.0 - data["surfacefractions"]["ground"]["fveg"] -
    data["surfacefractions"]["ground"]["fbare"]

# Vegetation
data["vegetation"] = Dict{String,Any}()
## Roof
data["vegetation"]["roof"] = Dict{String,Any}(
    "LAI" => 4.0,
    "SAI" => 0.001,
    "hc" => 0.05,
    "d_leaf" => 2.0,
    "CASE_ROOT" => 1,
    "ZR95" => [95.0],
    "ZR50" => [NaN],
    "ZRmax" => [NaN],
    "Rrootl" => [3000.0],
    "PsiL50" => [-4.0],
    "PsiX50" => [-4.5],
    "FI" => 0.081,
    "Do" => 1000.0,
    "a1" => 6.0,
    "go" => 0.01,
    "CT" => 3,
    "DSE" => 0.656,
    "Ha" => 55.0,
    "gmes" => Inf,
    "rjv" => 2.4,
    "Kopt" => 0.5,
    "Knit" => 0.15,
    "Vmax" => 96.0,
    "mSl" => 0.0,
    "e_rel" => 1.0,
    "e_relN" => 1.0,
    "Psi_sto_00" => -0.5,
    "Psi_sto_50" => -3.0,
    "Sl" => 0.035,
)

data["vegetation"]["roof"]["h_disp"] = 2.0 / 3.0 * data["vegetation"]["roof"]["hc"]

## Ground
data["vegetation"]["ground"] = Dict{String,Any}(
    "LAI" => 2.5,
    "SAI" => 0.001,
    "hc" => 0.05,
    "d_leaf" => 2.0,
    "CASE_ROOT" => 1,
    "ZR95" => [300.0],
    "ZR50" => [NaN],
    "ZRmax" => [NaN],
    "Rrootl" => [4000.0],
    "PsiL50" => [-2.0],
    "PsiX50" => [-5.5],
    "FI" => 0.004,
    "Do" => 2000.0,
    "a1" => 5.0,
    "go" => 0.01,
    "CT" => 4,
    "DSE" => 0.649,
    "Ha" => 72.0,
    "gmes" => Inf,
    "rjv" => 2.1,
    "Kopt" => 0.5,
    "Knit" => 0.3,
    "Vmax" => 54.0,
    "mSl" => 0.0,
    "e_rel" => 1.0,
    "e_relN" => 1.0,
    "Psi_sto_00" => -0.5,
    "Psi_sto_50" => -1.6,
    "Sl" => 0.025,
)

data["vegetation"]["ground"]["h_disp"] = 2.0 / 3.0 * data["vegetation"]["ground"]["hc"]

## Tree
data["vegetation"]["tree"] = Dict{String,Any}(
    "LAI" => 3.0,
    "SAI" => 0.2,
    "hc" => 0.05,
    "d_leaf" => 5.0,
    "CASE_ROOT" => 1,
    "ZR95" => [1500.0],
    "ZR50" => [NaN],
    "ZRmax" => [NaN],
    "Rrootl" => [2200.0],
    "PsiL50" => [-2.8],
    "PsiX50" => [-4.5],
    "FI" => 0.081,
    "Do" => 2000.0,
    "a1" => 9.0,
    "go" => 0.01,
    "CT" => 3,
    "DSE" => 0.649,
    "Ha" => 72.0,
    "gmes" => Inf,
    "rjv" => 2.2,
    "Kopt" => 0.5,
    "Knit" => 0.4,
    "Vmax" => 49.0,
    "mSl" => 0.0,
    "e_rel" => 1.0,
    "e_relN" => 1.0,
    "Psi_sto_00" => -0.9,
    "Psi_sto_50" => -1.7,
    "Sl" => 0.02,
    "hc" => NaN,
)

# Thermal properties
data["thermal"] = Dict{String,Any}()
data["thermal"]["roof"] = Dict{String,Any}("lan_dry" => 0.1, "cv_s" => 1.26e6)

data["thermal"]["ground"] = Dict{String,Any}("lan_dry" => 1.5, "cv_s" => 1.5e6)

data["thermal"]["wall"] = Dict{String,Any}("lan_dry" => 0.28, "cv_s" => 1.7e6)
data["thermal"]["tree"] = Dict{String,Any}("Cthermal_leaf" => 640.0)

data["optical"] = Dict{String,Any}()
data["optical"]["wall"] = Dict{String,Any}("albedo" => 0.3, "emissivity" => 0.97)
LAI_T = data["vegetation"]["tree"]["LAI"]
SAI_T = data["vegetation"]["tree"]["SAI"]
data["optical"]["tree"] = Dict{String,Any}(
    "albedo" => 0.2, "emissivity" => 1 - exp(-(LAI_T + SAI_T))
)
LAI_R = data["vegetation"]["roof"]["LAI"]
SAI_R = data["vegetation"]["roof"]["SAI"]
data["optical"]["roof"] = Dict{String,Any}(
    "aveg" => 0.2, "aimp" => 0.2, "eveg" => 1 - exp(-(LAI_R + SAI_R)), "eimp" => 0.97
)

LAI_G = data["vegetation"]["ground"]["LAI"]
SAI_G = data["vegetation"]["ground"]["SAI"]
data["optical"]["ground"] = Dict{String,Any}(
    "aveg" => 0.2,
    "abare" => 0.2,
    "aimp" => 0.1,
    "eveg" => 1 - exp(-(LAI_G + SAI_G)),
    "ebare" => 0.95,
    "eimp" => 0.95,
)

data["soil"] = Dict{String,Any}()

# Roof soil parameters
data["soil"]["roof"] = Dict{String,Any}(
    "Pcla" => 0.20,
    "Psan" => 0.40,
    "Porg" => 0.025,
    "In_max_imp" => 0.25,
    "In_max_ground" => 10.0,
    "Sp_In" => 0.2,
    "Kimp" => 0.0,
    "Kfc" => 0.2,
    "Phy" => 10000.0,
    "SPAR" => 2,
    "Kbot" => NaN,
    "dz1" => 0.105,
    "dz2" => 0.105,
    "Zs" => [0.0, 10.0, 20.0, 50.0, 100.0],
    "FixSM" => true,
    "FixSM_LayerStart" => 1,
    "FixSM_LayerEnd" => 4,
)

# Ground soil parameters
data["soil"]["ground"] = Dict{String,Any}(
    "Pcla" => 0.20,
    "Psan" => 0.40,
    "Porg" => 0.025,
    "In_max_imp" => 0.25,
    "In_max_underveg" => 10.0,
    "In_max_bare" => 10.0,
    "Sp_In" => 0.2,
    "Kimp" => 0.001,
    "Kfc" => 0.2,
    "Phy" => 10000.0,
    "SPAR" => 2,
    "Kbot" => NaN,
    "Zs" => [
        0.0,
        10.0,
        20.0,
        50.0,
        100.0,
        150.0,
        200.0,
        300.0,
        400.0,
        600.0,
        800.0,
        1000.0,
        1500.0,
        2000.0,
    ],
    "FixSM" => true,
    "FixSM_LayerStart" => 11,
    "FixSM_LayerEnd" => 13,
)

# Wall soil parameters
data["soil"]["wall"] = Dict{String,Any}("dz1" => 0.11, "dz2" => 0.11)

# Tree interception parameter
data["soil"]["Sp_In_T"] = 0.2

data["building_energy"] = Dict{String,Any}()
data["building_energy"]["indoor_optical"] = Dict{String,Any}(
    "abc" => 0.3,
    "abw" => 0.3,
    "abg" => 0.3,
    "abm" => 0.3,
    "ec" => 0.95,
    "eg" => 0.95,
    "ew" => 0.95,
    "em" => 0.95,
)
data["building_energy"]["thermal"] = Dict{String,Any}(
    "IntMassOn" => false,
    "FloorHeight" => 3.0,
    "dzFloor" => 0.2,
    "dzWall" => 0.2,
    "lan_ground_floor" => 1.2,
    "cv_ground_floor" => 1.5e6,
    "lan_floor_IntMass" => 0.67,
    "cv_floor_IntMass" => 1.0e6,
    "lan_wall_IntMass" => 0.67,
    "cv_wall_IntMass" => 1.0e6,
)
data["building_energy"]["windows"] = Dict{String,Any}(
    "WindowsOn" => 1,
    "GlazingRatio" => 0.15,
    "Uvalue" => 4.95,
    "lan_windows" => NaN,
    "cv_glass" => 2.1e6,
    "dztot" => 0.02,
    "SHGC" => 0.8,
    "SolarTransmittance" => 0.6,
    "SolarAbsorptivity" => 0.0,
    "SolarAlbedo" => 0.4,
)
data["building_energy"]["hvac"] = Dict{String,Any}(
    "ACon" => true,
    "Heatingon" => false,
    "TsetpointCooling" => 298.15,
    "TsetpointHeating" => 293.15,
    "RHsetpointCooling" => 60.0,
    "RHsetpointHeating" => NaN,
    "ACH" => 0.5,
    "COPAC" => 3.26,
    "COPHeat" => 0.9,
    "f_ACLatentToQ" => 1.0,
)

# Add person parameters
data["person"] = Dict{String,Any}(
    "PositionPx" => data["urbangeometry"]["Width_canyon"] / 2,
    "PositionPz" => 1.1,
    "PersonWidth" => 0.06 / 2,
    "PersonHeight" => 0.22 / 2,
    "HeightWind" => 1.1,
)

data["location"] = Dict{String,Any}(
    "phi" => 28.6, "lambda" => 77.1, "theta_canyon" => deg2rad(45), "DeltaGMT" => 5.0
)

YAML.write_file(joinpath(@__DIR__, "..", "data", "newdelhi_parameters.yaml"), data)

## NetCDF section
input_data = matread(joinpath(@__DIR__, "..", "data", "TMYNewDelhi_RadPart.mat"))
input_data["Time"] = [
    DateTime(
        input_data["Time"][i, 1],
        input_data["Time"][i, 2],
        input_data["Time"][i, 3],
        input_data["Time"][i, 4],
        input_data["Time"][i, 5],
    ) for i in axes(input_data["Time"], 1)
]
input_data["RelativeHumidity"] ./= 100.0

filename = "newdelhi_data.nc"
filepath = joinpath(@__DIR__, "..", "data", filename)

isfile(filepath) && rm(filepath)

ds = NCDataset(filepath, "c") do ds
    defDim(ds, "hours", length(input_data["Time"]))
    defVar(ds, "datetime", input_data["Time"], ("hours",))

    # Meteorological inputs
    defVar(ds, "LWR_in", input_data["LWRin"], ("hours",))
    defVar(ds, "SAB1_in", input_data["SAB1"], ("hours",))
    defVar(ds, "SAB2_in", input_data["SAB2"], ("hours",))
    defVar(ds, "SAD1_in", input_data["SAD1"], ("hours",))
    defVar(ds, "SAD2_in", input_data["SAD2"], ("hours",))
    defVar(ds, "Tatm", input_data["Tatm"] .+ 273.15, ("hours",))
    defVar(ds, "Uatm", input_data["Windspeed"], ("hours",))
    defVar(ds, "Pre", input_data["Pressure_Pa"], ("hours",))
    defVar(ds, "Rain", input_data["Rain"], ("hours",))
    defVar(ds, "rel_hum", input_data["RelativeHumidity"], ("hours",))
    defVar(ds, "Zatm", 30.0, ())
    defVar(ds, "Catm_CO2", 400.0, ())
    defVar(ds, "Catm_O2", 210000.0, ())
    defVar(ds, "SunDSM_MRT", NaN, ())

    # Anthropogenic inputs
    defVar(ds, "Tbmin", 20.0, ())
    defVar(ds, "Tbmax", 25.0, ())
    defVar(ds, "Qf_canyon", 0.0, ())

    # Sun position inputs
    defVar(ds, "t_bef", 0.5, ())
    defVar(ds, "t_aft", 0.5, ())
end
