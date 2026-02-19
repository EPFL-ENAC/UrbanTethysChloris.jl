"""
    ModelOptions <: AbstractOptions

Structure that contains the model options

Available options
- `mc_sample_size` and `n_rays`: Number of samples for Monte Carlo ray tracing and number of
rays per sample, respectively. These options control the accuracy and computational cost of
the ray tracing calculations.
- `RESPreCalc` and `fconvPreCalc`: Boolean flags to enable or disable pre-calculation of
resistance and convective heat transfer coefficients, respectively.
- `BEM_on`: A boolean flag to enable or disable the building energy model. When enabled,
the model will include calculations related to building energy balance, which can provide
more detailed insights into the energy dynamics of the urban environment.
- `output_level`: An option to specify the level of detail for the outputs saved during the simulation.
This can be set to different levels (e.g., `no_outputs`, `essential_outputs`,
`extended_energy_climate_outputs`, `extended_outputs`) to control the amount of data
generated and saved, which can be useful for managing storage and focusing on specific
aspects of the simulation results.
"""
Base.@kwdef struct ModelOptions <: AbstractModelOptions
    mc_sample_size::Int = 1000
    n_rays::Int = 200
    RESPreCalc::Bool = true
    fconvPreCalc::Bool = true
    BEM_on::Bool = true
    output_level::UrbanTethysChloris.ModelComponents.AbstractOutputsToSave=plot_outputs
    OPT_Obhukov::AbstractZeroFindingStrategies = SimpleBrentStrategy(
        Float64; abstol=1e-6, maxiters=400
    )
    OPT_SM::AbstractODEOptions = ODEOptions(abstol=0.05)
end
