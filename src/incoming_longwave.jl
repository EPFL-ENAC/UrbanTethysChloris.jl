
"""
    incoming_longwave(Ta::FT, ea::FT, N::FT) where {FT<:AbstractFloat}

Calculate incoming longwave radiation.

# Arguments
- `Ta`: Air temperature [°C]
- `ea`: Vapor pressure [Pa]
- `N`: Cloudiness [0-1] or incoming longwave radiation if N > 1 [W/m²]

# Returns
- `Latm`: Incoming longwave radiation [W/m²]
"""
function incoming_longwave(Ta::FT, ea::FT, N::FT) where {FT<:AbstractFloat}
    if N ≤ one(FT)  # N is cloudiness
        # Convert air temperature to Kelvin
        Ta_k = Ta + 273.15

        # Constants
        σ = 5.6704e-8  # Stefan-Boltzmann constant [W/m² K⁴]

        # Calculate precipitable water [kg/m²] (Prata 1996)
        w = 4.65 * ea / Ta_k

        # Calculate clear sky emissivity (Dilley and O'Brien 1998)
        e_cs = (59.38 + 113.7 * (Ta_k/273.16)^6 + 96.96 * sqrt(w/25)) / (σ * Ta_k^4)

        # Calculate cloud cover attenuation (Unsworth and Monteith 1975)
        K = (1 - 0.84 * N) + (0.84 * N) / e_cs

        # Calculate incoming longwave radiation [W/m²]
        Latm = K * e_cs * σ * Ta_k^4
    else
        # N is the incoming longwave radiation
        Latm = N
    end

    return Latm
end
