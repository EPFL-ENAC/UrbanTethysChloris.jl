using Test
using UrbanTethysChloris.ModelComponents.Parameters:
    VegetationParameters,
    HeightDependentVegetationParameters,
    initialize_vegetationparameters,
    initialize_heightdependent_vegetationparameters

FT = Float64

# Base test case for height dependent parameters
base_heightdep_data = Dict{String,Any}(
    "LAI" => 2.0,
    "SAI" => 0.5,
    "hc" => 1.0,
    "h_disp" => 0.7,
    "d_leaf" => 0.02,
    "CASE_ROOT" => 1,
    "ZR95" => [0.5],
    "ZR50" => [0.3],
    "ZRmax" => [1.0],
    "Rrootl" => [5e-4],
    "PsiL50" => [-2.0],
    "PsiX50" => [-4.0],
    "FI" => 0.08,
    "Do" => 1000.0,
    "a1" => 10.0,
    "go" => 0.01,
    "CT" => 3,
    "DSE" => 0.649,
    "Ha" => 72000.0,
    "gmes" => 0.02,
    "rjv" => 0.01,
    "Kopt" => 0.5,
    "Knit" => 0.5,
    "Vmax" => 50.0,
    "mSl" => 5.0,
    "e_rel" => 0.0001,
    "e_relN" => 0.3,
    "Psi_sto_00" => -0.5,
    "Psi_sto_50" => -2.0,
    "Sl" => 10.0,
    "SPARTREE" => 1,
)

# Base test case for vegetation parameters
base_data = Dict{String,Any}(
    "roof" => base_heightdep_data,
    "ground" => base_heightdep_data,
    "tree" => base_heightdep_data,
)

# Test empty dictionary
@test_throws ArgumentError initialize_vegetationparameters(FT, Dict{String,Any}())

# Test missing component
for key in ["roof", "ground", "tree"]
    invalid_data = copy(base_data)
    delete!(invalid_data, key)
    @test_throws ArgumentError initialize_vegetationparameters(FT, invalid_data)
end

additional_data = copy(base_data)
additional_data["toast"] = base_heightdep_data
@test_throws ArgumentError initialize_vegetationparameters(FT, additional_data)

# Test successful initialization
@test_nowarn initialize_vegetationparameters(FT, base_data)

# Test height dependent parameters separately
@test_throws ArgumentError initialize_heightdependent_vegetationparameters(
    FT, Dict{String,Any}()
)
@test_nowarn initialize_heightdependent_vegetationparameters(FT, base_heightdep_data)

# Test that SPARTREE is optional
missing_spartree = copy(base_heightdep_data)
delete!(missing_spartree, "SPARTREE")
@test_nowarn initialize_heightdependent_vegetationparameters(FT, missing_spartree)
