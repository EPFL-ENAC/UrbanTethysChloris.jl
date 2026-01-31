"""
    root_soil_conductance(
        Ks::FT,
        Rl::FT,
        rcyl::FT,
        rroot::FT,
        Zr::FT
    ) where {FT<:AbstractFloat}

Calculate soil-to-root conductance using the Holtta/Sperry model.

# Arguments
- `Ks`: Saturated hydraulic conductivity [mm/h]
- `Rl`: Root length index [m root/m² ground]
- `rcyl`: Radius of the cylinder of soil to which root has access [m]
- `rroot`: Root radius [m]
- `Zr`: Rooting depth [m]

# Returns
- `Ksr`: Soil-to-root conductance [mmol H₂O / m² ground s MPa]
"""
function root_soil_conductance(
    Ks::FT, Rl::FT, rcyl::FT, rroot::FT, Zr::FT
) where {FT<:AbstractFloat}
    Zr = Zr / 1000

    row = 1000
    g = FT(9.81)
    rho = 55555
    CF = 1e6 / (row * g)
    Ks = Ks / 3600000

    gsr = Ks * Rl * (2 * π) / (log(rcyl / rroot))
    Ksr = gsr * CF * row * rho

    return Ksr
end
