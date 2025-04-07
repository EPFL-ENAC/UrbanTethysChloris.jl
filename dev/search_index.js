var documenterSearchIndex = {"docs":
[{"location":"91-developer/#dev_docs","page":"Developer documentation","title":"Developer documentation","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"note: Contributing guidelines\nIf you haven't, please read the Contributing guidelines first.","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"If you want to make contributions to this package that involves code, then this guide is for you.","category":"page"},{"location":"91-developer/#First-time-clone","page":"Developer documentation","title":"First time clone","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"tip: If you have writing rights\nIf you have writing rights, you don't have to fork. Instead, simply clone and skip ahead. Whenever upstream is mentioned, use origin instead.","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"If this is the first time you work with this repository, follow the instructions below to clone the repository.","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Fork this repo\nClone your repo (this will create a git remote called origin)\nAdd this repo as a remote:\ngit remote add upstream https://github.com/EPFL-ENAC/UrbanTethysChloris.jl","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"This will ensure that you have two remotes in your git: origin and upstream. You will create branches and push to origin, and you will fetch and update your local main branch from upstream.","category":"page"},{"location":"91-developer/#Linting-and-formatting","page":"Developer documentation","title":"Linting and formatting","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Install a plugin on your editor to use EditorConfig. This will ensure that your editor is configured with important formatting settings.","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"We use https://pre-commit.com to run the linters and formatters. In particular, the Julia code is formatted using JuliaFormatter.jl, so please install it globally first:","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"julia> # Press ]\npkg> activate\npkg> add JuliaFormatter","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"To install pre-commit, we recommend using pipx as follows:","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"# Install pipx following the link\npipx install pre-commit","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"With pre-commit installed, activate it as a pre-commit hook:","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"pre-commit install","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"To run the linting and formatting manually, enter the command below:","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"pre-commit run -a","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Now, you can only commit if all the pre-commit tests pass.","category":"page"},{"location":"91-developer/#Testing","page":"Developer documentation","title":"Testing","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"As with most Julia packages, you can just open Julia in the repository folder, activate the environment, and run test:","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"julia> # press ]\npkg> activate .\npkg> test","category":"page"},{"location":"91-developer/#Working-on-a-new-issue","page":"Developer documentation","title":"Working on a new issue","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"We try to keep a linear history in this repo, so it is important to keep your branches up-to-date.","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Fetch from the remote and fast-forward your local main\ngit fetch upstream\ngit switch main\ngit merge --ff-only upstream/main\nBranch from main to address the issue (see below for naming)\ngit switch -c 42-add-answer-universe\nPush the new local branch to your personal remote repository\ngit push -u origin 42-add-answer-universe\nCreate a pull request to merge your remote branch into the org main.","category":"page"},{"location":"91-developer/#Branch-naming","page":"Developer documentation","title":"Branch naming","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"If there is an associated issue, add the issue number.\nIf there is no associated issue, and the changes are small, add a prefix such as \"typo\", \"hotfix\", \"small-refactor\", according to the type of update.\nIf the changes are not small and there is no associated issue, then create the issue first, so we can properly discuss the changes.\nUse dash separated imperative wording related to the issue (e.g., 14-add-tests, 15-fix-model, 16-remove-obsolete-files).","category":"page"},{"location":"91-developer/#Commit-message","page":"Developer documentation","title":"Commit message","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Use imperative or present tense, for instance: Add feature or Fix bug.\nHave informative titles.\nWhen necessary, add a body with details.\nIf there are breaking changes, add the information to the commit message.","category":"page"},{"location":"91-developer/#Before-creating-a-pull-request","page":"Developer documentation","title":"Before creating a pull request","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"tip: Atomic git commits\nTry to create \"atomic git commits\" (recommended reading: The Utopic Git History).","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Make sure the tests pass.\nMake sure the pre-commit tests pass.\nFetch any main updates from upstream and rebase your branch, if necessary:\ngit fetch upstream\ngit rebase upstream/main BRANCH_NAME\nThen you can open a pull request and work with the reviewer to address any issues.","category":"page"},{"location":"91-developer/#Building-and-viewing-the-documentation-locally","page":"Developer documentation","title":"Building and viewing the documentation locally","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Following the latest suggestions, we recommend using LiveServer to build the documentation. Here is how you do it:","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Run julia --project=docs to open Julia in the environment of the docs.\nIf this is the first time building the docs\nPress ] to enter pkg mode\nRun pkg> dev . to use the development version of your package\nPress backspace to leave pkg mode\nRun julia> using LiveServer\nRun julia> servedocs()","category":"page"},{"location":"91-developer/#Making-a-new-release","page":"Developer documentation","title":"Making a new release","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"To create a new release, you can follow these simple steps:","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Create a branch release-x.y.z\nUpdate version in Project.toml\nUpdate the CHANGELOG.md:\nRename the section \"Unreleased\" to \"[x.y.z] - yyyy-mm-dd\" (i.e., version under brackets, dash, and date in ISO format)\nAdd a new section on top of it named \"Unreleased\"\nAdd a new link in the bottom for version \"x.y.z\"\nChange the \"[unreleased]\" link to use the latest version - end of line, vx.y.z ... HEAD.\nCreate a commit \"Release vx.y.z\", push, create a PR, wait for it to pass, merge the PR.\nGo back to main screen and click on the latest commit (link: https://github.com/EPFL-ENAC/UrbanTethysChloris.jl/commit/main)\nAt the bottom, write @JuliaRegistrator register","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"After that, you only need to wait and verify:","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Wait for the bot to comment (should take < 1m) with a link to a PR to the registry\nFollow the link and wait for a comment on the auto-merge\nThe comment should said all is well and auto-merge should occur shortly\nAfter the merge happens, TagBot will trigger and create a new GitHub tag. Check on https://github.com/EPFL-ENAC/UrbanTethysChloris.jl/releases\nAfter the release is create, a \"docs\" GitHub action will start for the tag.\nAfter it passes, a deploy action will run.\nAfter that runs, the stable docs should be updated. Check them and look for the version number.","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"","category":"page"},{"location":"95-reference/#reference","page":"Reference","title":"Reference","text":"","category":"section"},{"location":"95-reference/#Contents","page":"Reference","title":"Contents","text":"","category":"section"},{"location":"95-reference/","page":"Reference","title":"Reference","text":"Pages = [\"95-reference.md\"]","category":"page"},{"location":"95-reference/#Index","page":"Reference","title":"Index","text":"","category":"section"},{"location":"95-reference/","page":"Reference","title":"Reference","text":"Pages = [\"95-reference.md\"]","category":"page"},{"location":"95-reference/","page":"Reference","title":"Reference","text":"Modules = [UrbanTethysChloris]","category":"page"},{"location":"95-reference/#TethysChlorisCore.validate_fields-Union{Tuple{T}, Tuple{Type{T}, Dict{String, Any}}} where T<:TethysChlorisCore.ModelComponents.AbstractModelComponent","page":"Reference","title":"TethysChlorisCore.validate_fields","text":"TethysChlorisCore.validate_fields(::Type{T}, data::Dict{String,Any}) where {T<:AbstractModelComponent}\n\nValidate that a dictionary contains only the required fields for a given model component type.\n\nArguments\n\nT: Type of model component to validate fields against\ndata: Dictionary containing fields to validate\n\nThrows\n\nArgumentError: If any extraneous fields are found in the data dictionary\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#UrbanTethysChloris.check_extraneous_fields-Union{Tuple{N}, Tuple{T}, Tuple{Type{T}, Dict{String, Any}, Union{NTuple{N, String}, Vector{String}}}} where {T<:TethysChlorisCore.ModelComponents.AbstractModelComponent, N}","page":"Reference","title":"UrbanTethysChloris.check_extraneous_fields","text":"check_extraneous_fields(::Type{T}, data::Dict{String,Any}) where {T<:AbstractModelComponent}\n\nCheck if the input dictionary contains any keys that are not part of the required fields for type T.\n\nArguments\n\nT: Type of model component to check fields against\ndata: Dictionary containing fields to validate\nrequired_fields: List of required fields for the model component. This is typically obtained from get_required_fields(T) or fieldnames(T)\n\nThrows\n\nArgumentError: If any extraneous keys are found in the data dictionary\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#UrbanTethysChloris.get_calculated_fields-Union{Tuple{Type{T}}, Tuple{T}} where T<:TethysChlorisCore.ModelComponents.AbstractModelComponent","page":"Reference","title":"UrbanTethysChloris.get_calculated_fields","text":"get_calculated_fields(::Type{T}) where {T<:AbstractModelComponent}\n\nGet a list of calculated fields for a given model component type. Calculated fields are those that are computed based on other fields and should not be provided when creating a model component. By default, returns an empty list. Components should override this method if they have calculated fields.\n\nArguments\n\nT: Type of model component to get calculated fields for\n\nReturns\n\nVector{Symbol}: List of calculated field names\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#UrbanTethysChloris.get_optional_fields-Union{Tuple{Type{T}}, Tuple{T}} where T<:TethysChlorisCore.ModelComponents.AbstractModelComponent","page":"Reference","title":"UrbanTethysChloris.get_optional_fields","text":"get_optional_fields(::Type{T}) where {T<:AbstractModelComponent}\n\nGet a list of optional fields for a given model component type. Optional fields are those that can be omitted when creating a model component. By default, returns an empty list. Components should override this method if they have optional fields.\n\nArguments\n\nT: Type of model component to get optional fields for\n\nReturns\n\nVector{Symbol}: List of optional field names\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#Model-Components","page":"Reference","title":"Model Components","text":"","category":"section"},{"location":"95-reference/#Parameters","page":"Reference","title":"Parameters","text":"","category":"section"},{"location":"95-reference/","page":"Reference","title":"Reference","text":"UrbanTethysChloris.ModelComponents.Parameters.HeightDependentVegetationParameters\nUrbanTethysChloris.ModelComponents.Parameters.SurfaceFractions\nUrbanTethysChloris.ModelComponents.Parameters.LocationSpecificThermalProperties\nUrbanTethysChloris.ModelComponents.Parameters.UrbanGeometryParameters\nUrbanTethysChloris.ModelComponents.Parameters.VegetationParameters\nUrbanTethysChloris.ModelComponents.Parameters.ThermalProperties\nUrbanTethysChloris.ModelComponents.Parameters.OpticalProperties\nUrbanTethysChloris.ModelComponents.Parameters.SimpleOpticalProperties\nUrbanTethysChloris.ModelComponents.Parameters.VegetatedOpticalProperties\nUrbanTethysChloris.ModelComponents.Parameters.BuildingEnergyModelParameters\nUrbanTethysChloris.ModelComponents.Parameters.IndoorOpticalProperties\nUrbanTethysChloris.ModelComponents.Parameters.ThermalBuilding\nUrbanTethysChloris.ModelComponents.Parameters.WindowParameters\nUrbanTethysChloris.ModelComponents.Parameters.HVACParameters\nUrbanTethysChloris.ModelComponents.Parameters.SoilParameters\nUrbanTethysChloris.ModelComponents.Parameters.VegetatedSoilParameters\nUrbanTethysChloris.ModelComponents.Parameters.PersonParameters\nUrbanTethysChloris.ModelComponents.Parameters.ParameterSet\nUrbanTethysChloris.ModelComponents.Parameters.TreeThermalProperties","category":"page"},{"location":"95-reference/#UrbanTethysChloris.ModelComponents.Parameters.HeightDependentVegetationParameters","page":"Reference","title":"UrbanTethysChloris.ModelComponents.Parameters.HeightDependentVegetationParameters","text":"HeightDependentVegetationParameters{FT<:AbstractFloat} <: AbstractParameters{FT}\n\nParameters for the location-specific (tree, roof, ground) vegetation.\n\nFields\n\nLAI::FT: Leaf area index [-]\nSAI::FT: Stem area index [-].\nhc::FT: Canopy height [m].\nh_disp::FT:  [-].\nd_leaf::FT:  [-].\nCASE_ROOT::Int: Type of root profile [-].\nZR95::FT: Root depth, 95th percentile [mm].\nZR50::FT: Root depth, 50th percentile [mm].\nZRmax::FT: Root depth, maximum [mm]\nRrootl::FT: Root length index [m root m-2 PFT].\nPsiL50::FT: Water potential at 50 % of leaf hydraulic conductivity [MPa].\nPsiX50::FT: Water potential at 50 % of xylem hydraulic conductivity and limit for water extraction from soil [MPa].\nFI::FT: Intrinsec quantum efficency [µmol CO2 µmol-1 photons].\nDo::FT: Empirical coefficient that expresses the value of vapor pressure deficit at which f(∆e) = 0.5 [Pa].\na1::FT: Empirical parameter linking net assimilaton AnCto stomatal conductanceg{s,CO2}` [-].\ngo::FT: Minimum/cuticular stomatal conductance [mol CO2 m-2 leaf s-1].\nCT::Int: Photosyntesis Typology for Plants, Photosynthetic pathway C3 or C4 [-]\nDSE::FT: Activation Energy - Plant Dependent, Activation Energy in Photosynthesis for Rubisco Capacity [kJ/mol].\nHa::FT: Entropy factor - Plant Dependent, Activation energy [kJ / mol K].\ngmes::FT: Mesophyll conductance, not used [mol CO2/s m2].\nrjv::FT: Scaling factor between Jmax and V_{c,max} [µmol equivalent µmol-1 CO2].\nKopt::FT: Optical depth of direct beam per unit plant area (?) [-].\nKnit::FT: Canopy nitrogen decay coefficient [-].\nVmax::FT: Maximum Rubisco capacity at 25◦C leaf scale [µmol CO2 m-2 s-1]\nmSl::FT:  [-].\ne_rel::FT: Relative Efficiency of the photosynthesis apparatus due to Age/Day-length [-].\ne_relN::FT: Relative efficiency of the photosynthesis apparatus due to N limitations [-].\nPsi_sto_00::FT: Soil water potential at the beginning of stomatal closure [MPa].\nPsi_sto_50::FT: Soil water potential at 50 % stomatal closure [MPa].\nSl::FT: Specific leaf area of biomass [m^2 /gC] [-].\nSPARTREE::Int: Tree root distribution: 1 = Tree roots can access all water in the soil\n\n(imp, bare, veg) equally; 2 = If the tree crown is smaller than the combined vegetated and bare fraction, then the trees only transpire from these fractions. Otherwise, they also transpire from the impervious ground fraction. [-].\n\n\n\n\n\n","category":"type"},{"location":"95-reference/#UrbanTethysChloris.ModelComponents.Parameters.SurfaceFractions","page":"Reference","title":"UrbanTethysChloris.ModelComponents.Parameters.SurfaceFractions","text":"SurfaceFractions{FT<:AbstractFloat} <: AbstractParameters{FT}\n\nParameters for the SurfaceFractions model component.\n\nFields\n\nroof::LocationSpecificSurfaceFractions{FT}: Roof surface fractions.\nground::LocationSpecificSurfaceFractions{FT}: Ground surface fractions.\n\n\n\n\n\n","category":"type"},{"location":"95-reference/#UrbanTethysChloris.ModelComponents.Parameters.LocationSpecificThermalProperties","page":"Reference","title":"UrbanTethysChloris.ModelComponents.Parameters.LocationSpecificThermalProperties","text":"LocationSpecificThermalProperties{FT<:AbstractFloat} <: AbstractHeightDependentParameters{FT}\n\nParameters for the location-specific (wall, roof, ground) thermal properties.\n\nFields\n\nlan_dry::FT: Thermal conductivity dry solid [W/m K]\ncv_s::FT: Volumetric heat capacity solid [J/m^3 K].\n\n\n\n\n\n","category":"type"},{"location":"95-reference/#UrbanTethysChloris.ModelComponents.Parameters.UrbanGeometryParameters","page":"Reference","title":"UrbanTethysChloris.ModelComponents.Parameters.UrbanGeometryParameters","text":"UrbanGeometryParameters{FT<:AbstractFloat} <: AbstractParameters{FT}\n\nParameters for the UrbanGeometry model component.\n\nFields\n\nHeight_canyon::FT: Height of urban canyon [m].\nWidth_canyon::FT: Ground width of urban canyon [m].\nWidth_roof::FT: Roof width of urban canyon [m].\nHeight_tree::FT: Tree height [m].\nRadius_tree::FT: Tree radius (=1/4 fg,tree *Wcan) [m].\nDistance_tree::FT: Distance of wall to tree trunk [m].\nHcan_max::FT: Maximum height of roughness elements (buildings) [m].\nHcan_std::FT: Standard deviation of roughness elements (buildings) [m].\ntrees::Bool: Easy switch to include (=1) and exclude (=0) trees in the urban canyon.\nftree::FT: Tree fraction along canyon axis\nhcanyon::FT: Normalized height of urban canyon [-].\nwcanyon::FT: Normalized ground width of urban canyon [-].\nwroof::FT: Normalized roof width of urban canyon [-]\nhtree::FT: Normalized tree height [-]\nradius_tree::FT: Normalized tree radius [-].\ndistance_tree::FT: Normalized distance of wall to tree trunk [-].\nratio::FT: Height-to-width ratio [-].\nwcanyon_norm::FT: Normalized canyon width overall [-].\nwroof_norm::FT: Normalized roof width overall [-].\n\n\n\n\n\n","category":"type"},{"location":"95-reference/#UrbanTethysChloris.ModelComponents.Parameters.VegetationParameters","page":"Reference","title":"UrbanTethysChloris.ModelComponents.Parameters.VegetationParameters","text":"VegetationParameters{FT<:AbstractFloat} <: AbstractParameters{FT}\n\nContainer for vegetation parameters for different urban surface components.\n\nFields\n\nroof::HeightDependentVegetationParameters{FT}: Vegetation parameters for roof vegetation\nground::HeightDependentVegetationParameters{FT}: Vegetation parameters for ground-level vegetation\ntree::HeightDependentVegetationParameters{FT}: Vegetation parameters for trees\n\n\n\n\n\n","category":"type"},{"location":"95-reference/#UrbanTethysChloris.ModelComponents.Parameters.ThermalProperties","page":"Reference","title":"UrbanTethysChloris.ModelComponents.Parameters.ThermalProperties","text":"ThermalProperties{FT<:AbstractFloat} <: AbstractParameters{FT}\n\nContainer for vegetation parameters for different urban surface components.\n\nFields\n\nroof::LocationSpecificThermalProperties{FT}: Roof thermal properties\nground::LocationSpecificThermalProperties{FT}: Ground thermal properties\nwall::LocationSpecificThermalProperties{FT}: Wall thermal properties\n\n\n\n\n\n","category":"type"},{"location":"95-reference/#UrbanTethysChloris.ModelComponents.Parameters.OpticalProperties","page":"Reference","title":"UrbanTethysChloris.ModelComponents.Parameters.OpticalProperties","text":"OpticalProperties{FT<:AbstractFloat} <: AbstractParameters{FT}\n\nContainer for optical properties for different urban surface components.\n\nFields\n\nroof::VegetatedOpticalProperties{FT}: Roof optical properties\nground::VegetatedOpticalProperties{FT}: Ground optical properties\nwall::SimpleOpticalProperties{FT}: Wall optical properties\ntree::SimpleOpticalProperties{FT}: Tree optical properties\n\n\n\n\n\n","category":"type"},{"location":"95-reference/#UrbanTethysChloris.ModelComponents.Parameters.SimpleOpticalProperties","page":"Reference","title":"UrbanTethysChloris.ModelComponents.Parameters.SimpleOpticalProperties","text":"SimpleOpticalProperties{FT<:AbstractFloat} <: AbstractParameters{FT}\n\nParameters for the location-specific (wall and tree) optical properties, for which albedo and emissivity have been pre-calculated.\n\nFields\n\nalbedo::FT: Surface albedo (-)\nemissivity::FT: Surface emissivity (-)\n\n\n\n\n\n","category":"type"},{"location":"95-reference/#UrbanTethysChloris.ModelComponents.Parameters.VegetatedOpticalProperties","page":"Reference","title":"UrbanTethysChloris.ModelComponents.Parameters.VegetatedOpticalProperties","text":"VegetatedOpticalProperties{FT<:AbstractFloat} <: AbstractParameters{FT}\n\nOptical properties specific to vegetated surfaces (roof and ground)\n\nFields\n\naveg::FT: Vegetation surface albedo (-)\naimp::FT: Impervious surface albedo (-)\nabare::FT: Bare surface albedo (-), only used for ground\neveg::FT: Vegetation surface emissivity (-)\neimp::FT: Impervious surface emissivity (-)\nebare::FT: Bare surface emissivity (-), only used for ground\n\n\n\n\n\n","category":"type"},{"location":"95-reference/#UrbanTethysChloris.ModelComponents.Parameters.BuildingEnergyModelParameters","page":"Reference","title":"UrbanTethysChloris.ModelComponents.Parameters.BuildingEnergyModelParameters","text":"BuildingEnergyModelParameters{FT<:AbstractFloat} <: AbstractParameters{FT}\n\nContainer for all building energy model parameters.\n\nFields\n\nindoor_optical::IndoorOpticalProperties{FT}: Indoor surface optical properties\nthermal::ThermalBuilding{FT}: Building thermal properties\nwindows::WindowParameters{FT}: Window parameters\nhvac::HVACParameters{FT}: HVAC system parameters\n\n\n\n\n\n","category":"type"},{"location":"95-reference/#UrbanTethysChloris.ModelComponents.Parameters.IndoorOpticalProperties","page":"Reference","title":"UrbanTethysChloris.ModelComponents.Parameters.IndoorOpticalProperties","text":"IndoorOpticalProperties{FT<:AbstractFloat} <: AbstractParameters{FT}\n\nOptical properties for indoor building surfaces.\n\nFields\n\nabc::FT: Albedo ceiling (-)\nabw::FT: Albedo wall (-)\nabg::FT: Albedo ground (-)\nabm::FT: Albedo internal mass (-)\nec::FT: Emissivity ceiling (-)\neg::FT: Emissivity ground (-)\new::FT: Emissivity wall (-)\nem::FT: Emissivity internal mass (-)\n\n\n\n\n\n","category":"type"},{"location":"95-reference/#UrbanTethysChloris.ModelComponents.Parameters.ThermalBuilding","page":"Reference","title":"UrbanTethysChloris.ModelComponents.Parameters.ThermalBuilding","text":"ThermalBuilding{FT<:AbstractFloat} <: AbstractParameters{FT}\n\nParameters specifying thermal properties and dimensions of building internal mass.\n\nFields\n\nIntMassOn::Int: Include building internal mass in calculation (0/1)\nFloorHeight::FT: Average floor height in building (m)\ndzFloor::FT: Average thickness of floors in building (m)\ndzWall::FT: Average thickness of walls in building (m)\nlan_ground_floor::FT: Building ground thermal conductivity (W/m K)\ncv_ground_floor::FT: Building ground volumetric heat capacity (J/m³ K)\nlan_floor_IntMass::FT: Internal mass floor thermal conductivity (W/m K)\ncv_floor_IntMass::FT: Internal mass floor volumetric heat capacity (J/m³ K)\nlan_wall_IntMass::FT: Internal mass wall thermal conductivity (W/m K)\ncv_wall_IntMass::FT: Internal mass wall volumetric heat capacity (J/m³ K)\n\n\n\n\n\n","category":"type"},{"location":"95-reference/#UrbanTethysChloris.ModelComponents.Parameters.WindowParameters","page":"Reference","title":"UrbanTethysChloris.ModelComponents.Parameters.WindowParameters","text":"WindowParameters{FT<:AbstractFloat} <: AbstractParameters{FT}\n\nParameters for building windows.\n\nFields\n\nWindowsOn::Int: Include windows in simulation (0/1)\nGlazingRatio::FT: Window-to-wall ratio (-)\nUvalue::FT: U-value of windows (W/m² K)\nlan_windows::FT: Thermal conductivity of windows (W/m K)\ncv_glass::FT: Volumetric heat capacity of glass (J/m³ K)\ndztot::FT: Total thickness of glass layers (m)\nSHGC::FT: Solar heat gain coefficient (-)\nSolarTransmittance::FT: Solar radiation transmittance through windows (-)\nSolarAbsorptivity::FT: Fraction of solar radiation absorbed by window (-)\nSolarAlbedo::FT: Window albedo (-)\n\n\n\n\n\n","category":"type"},{"location":"95-reference/#UrbanTethysChloris.ModelComponents.Parameters.HVACParameters","page":"Reference","title":"UrbanTethysChloris.ModelComponents.Parameters.HVACParameters","text":"HVACParameters{FT<:AbstractFloat} <: AbstractParameters{FT}\n\nParameters for HVAC system.\n\nFields\n\nACon::Int: Enable air conditioning (0/1)\nHeatingon::Int: Enable heating (0/1)\nTsetpointCooling::FT: Cooling setpoint temperature (K)\nTsetpointHeating::FT: Heating setpoint temperature (K)\nRHsetpointCooling::FT: Cooling setpoint relative humidity (%)\nRHsetpointHeating::FT: Heating setpoint relative humidity (%)\nACH::FT: Air changes per hour (1/h)\nCOPAC::FT: Coefficient of performance for AC (-)\nCOPHeat::FT: Coefficient of performance for heating (-)\nf_ACLatentToQ::FT: Fraction of latent heat removed by AC that is condensed (-)\n\n\n\n\n\n","category":"type"},{"location":"95-reference/#UrbanTethysChloris.ModelComponents.Parameters.SoilParameters","page":"Reference","title":"UrbanTethysChloris.ModelComponents.Parameters.SoilParameters","text":"SoilParameters{FT<:AbstractFloat} <: AbstractParameters{FT}\n\nContainer for soil parameters for different urban surface components.\n\nFields\n\nroof::VegetatedSoilParameters{FT}: Roof soil parameters\nground::VegetatedSoilParameters{FT}: Ground soil parameters\nSp_In_T::FT: Specific water retained by a tree (mm m^2 VEG area m^-2 plant area)\n\n\n\n\n\n","category":"type"},{"location":"95-reference/#UrbanTethysChloris.ModelComponents.Parameters.VegetatedSoilParameters","page":"Reference","title":"UrbanTethysChloris.ModelComponents.Parameters.VegetatedSoilParameters","text":"VegetatedSoilParameters{FT<:AbstractFloat} <: AbstractParameters{FT}\n\nSoil parameters specific to vegetated surfaces (roof and ground)\n\nFields\n\nPcla::FT: Fraction of clay in the soil (-)\nPsan::FT: Fraction of sand in the soil (-)\nPorg::FT: Fraction of organic material in the soil (-)\nIn_max_imp::FT: Maxiumum interception capacity of impervious area (mm)\nIn_max_ground::FT: Maxiumum interception capacity of ground under roof vegetation (mm)\nIn_max_underveg::FT: Maxiumum interception capacity of vegetated ground area (mm)\nIn_max_bare::FT: Maxiumum interception capacity of bare ground area (mm)\nSp_In::FT: specific water retained by a vegetated surface (mm m^2 VEG area m^-2 plant area)\nKimp::FT: Hydraulic conductivity of impervious area (mm/h)\nKfc::FT: Conductivity at field capacity (mm/h)\nPhy::FT: Suction at the residual/hygroscopic water content (kPa)\nSPAR::Int: SOIL PARAMETER TYPE 1-VanGenuchten 2-Saxton-Rawls\nKbot::FT: Conductivity at the bedrock layer (mm/h)\n\n\n\n\n\n","category":"type"},{"location":"95-reference/#UrbanTethysChloris.ModelComponents.Parameters.PersonParameters","page":"Reference","title":"UrbanTethysChloris.ModelComponents.Parameters.PersonParameters","text":"PersonParameters{FT<:AbstractFloat} <: AbstractParameters{FT}\n\nParameters for the person in the urban canyon for MRT calculations.\n\nFields\n\nPositionPx::FT: Position within canyon [m]\nPositionPz::FT: Height of centre of person [m]\nPersonWidth::FT: Horizontal radius of ellipse describing person (=hip width / 2) [-]\nPersonHeight::FT: Vertical radius of ellipse describing person (= height / 2) [-]\nHeightWind::FT: Height for wind speed to calculate OTC [m]\n\n\n\n\n\n","category":"type"},{"location":"95-reference/#UrbanTethysChloris.ModelComponents.Parameters.ParameterSet","page":"Reference","title":"UrbanTethysChloris.ModelComponents.Parameters.ParameterSet","text":"ParameterSet{FT<:AbstractFloat} <: AbstractParameterSet{FT}\n\nParameters for the Urban Tethys-Chloris model.\n\nFields\n\nbuilding_energy::BuildingEnergyModelParameters{FT}: Parameters for the building energy model.\nperson::PersonParameters{FT}: Parameters for the person.\nsoil::SoilParameters{FT}: Parameters for the soil.\nsurfacefractions::SurfaceFractions{FT}: Parameters for the surface fractions.\nthermal::ThermalProperties{FT}: Parameters for the thermal properties.\noptical::OpticalProperties{FT}: Parameters for the optical properties.\nurbangeometry::UrbanGeometryParameters{FT}: Parameters for the urban geometry.\nvegetation::VegetationParameters{FT}: Parameters for the vegetation.\n\n\n\n\n\n","category":"type"},{"location":"95-reference/#UrbanTethysChloris.ModelComponents.Parameters.TreeThermalProperties","page":"Reference","title":"UrbanTethysChloris.ModelComponents.Parameters.TreeThermalProperties","text":"TreeThermalProperties{FT<:AbstractFloat} <: AbstractHeightDependentParameters{FT}\n\nParameters for the tree thermal properties.\n\nFields\n\nCthermal_leaf::FT: [J m-2 K-1] Heat capacity per single leaf area\n\n\n\n\n\n","category":"type"},{"location":"95-reference/","page":"Reference","title":"Reference","text":"","category":"page"},{"location":"90-contributing/#contributing","page":"Contributing guidelines","title":"Contributing guidelines","text":"","category":"section"},{"location":"90-contributing/","page":"Contributing guidelines","title":"Contributing guidelines","text":"First of all, thanks for the interest!","category":"page"},{"location":"90-contributing/","page":"Contributing guidelines","title":"Contributing guidelines","text":"We welcome all kinds of contribution, including, but not limited to code, documentation, examples, configuration, issue creating, etc.","category":"page"},{"location":"90-contributing/","page":"Contributing guidelines","title":"Contributing guidelines","text":"Be polite and respectful, and follow the code of conduct.","category":"page"},{"location":"90-contributing/#Bug-reports-and-discussions","page":"Contributing guidelines","title":"Bug reports and discussions","text":"","category":"section"},{"location":"90-contributing/","page":"Contributing guidelines","title":"Contributing guidelines","text":"If you think you found a bug, feel free to open an issue. Focused suggestions and requests can also be opened as issues. Before opening a pull request, start an issue or a discussion on the topic, please.","category":"page"},{"location":"90-contributing/#Working-on-an-issue","page":"Contributing guidelines","title":"Working on an issue","text":"","category":"section"},{"location":"90-contributing/","page":"Contributing guidelines","title":"Contributing guidelines","text":"If you found an issue that interests you, comment on that issue what your plans are. If the solution to the issue is clear, you can immediately create a pull request (see below). Otherwise, say what your proposed solution is and wait for a discussion around it.","category":"page"},{"location":"90-contributing/","page":"Contributing guidelines","title":"Contributing guidelines","text":"tip: Tip\nFeel free to ping us after a few days if there are no responses.","category":"page"},{"location":"90-contributing/","page":"Contributing guidelines","title":"Contributing guidelines","text":"If your solution involves code (or something that requires running the package locally), check the developer documentation. Otherwise, you can use the GitHub interface directly to create your pull request.","category":"page"},{"location":"90-contributing/","page":"Contributing guidelines","title":"Contributing guidelines","text":"","category":"page"},{"location":"","page":"UrbanTethysChloris","title":"UrbanTethysChloris","text":"CurrentModule = UrbanTethysChloris","category":"page"},{"location":"#UrbanTethysChloris","page":"UrbanTethysChloris","title":"UrbanTethysChloris","text":"","category":"section"},{"location":"","page":"UrbanTethysChloris","title":"UrbanTethysChloris","text":"Documentation for UrbanTethysChloris.","category":"page"},{"location":"#Contributors","page":"UrbanTethysChloris","title":"Contributors","text":"","category":"section"},{"location":"","page":"UrbanTethysChloris","title":"UrbanTethysChloris","text":"<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->\n<!-- prettier-ignore-start -->\n<!-- markdownlint-disable -->\n<table>\n  <tbody>\n    <tr>\n      <td align=\"center\" valign=\"top\" width=\"14.28%\"><a href=\"https://github.com/sphamba\"><img src=\"https://avatars.githubusercontent.com/u/17217484?v=4?s=100\" width=\"100px;\" alt=\"Son Pham-Ba\"/><br /><sub><b>Son Pham-Ba</b></sub></a><br /><a href=\"#projectManagement-sphamba\" title=\"Project Management\">📆</a> <a href=\"#review-sphamba\" title=\"Reviewed Pull Requests\">👀</a></td>\n    </tr>\n  </tbody>\n</table>\n\n<!-- markdownlint-restore -->\n<!-- prettier-ignore-end -->\n\n<!-- ALL-CONTRIBUTORS-LIST:END -->","category":"page"},{"location":"","page":"UrbanTethysChloris","title":"UrbanTethysChloris","text":"","category":"page"},{"location":"01-modelcomponents/#model_components","page":"Model Components","title":"Model Components","text":"","category":"section"},{"location":"01-modelcomponents/","page":"Model Components","title":"Model Components","text":"The model is structured using model components as composite types. Here is a general overview of the model components:","category":"page"},{"location":"01-modelcomponents/","page":"Model Components","title":"Model Components","text":"graph TD\n  A[ParameterSet] --> B[SurfaceFractions] --> F[LocationSpecificSurfaceFractions]\n  A --> C[ThermalProperties]\n  A --> D[UrbanGeometryParameters]\n  A --> E[VegetationParameters] --> G[HeightDependentVegetationParameters]\n  A --> H[Optical Properties] --> I[SimpleOpticalProperties]\n  H --> J[VegetatedOpticalProperties]\n  A --> K[BuildingEnergyModelParameters]\n  K --> L[IndoorOpticalProperties]\n  K --> M[ThermalBuilding]\n  K --> N[WindowParameters]\n  K --> O[HVACParameters]\n  A --> P[PersonParameters]\n  A --> Q[SoilParameters]\n  Q --> R[VegetatedSoilParameters]\n  C --> S[TreeThermalProperties]\n  C --> T[LocationSpecificThermalProperties]\n  Q --> U[WallSoilParameters]","category":"page"},{"location":"01-modelcomponents/","page":"Model Components","title":"Model Components","text":"","category":"page"}]
}
