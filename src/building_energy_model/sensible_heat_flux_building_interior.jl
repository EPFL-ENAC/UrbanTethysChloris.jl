"""
    sensible_heat_flux_building_interior(
        Tbin::FT,
        Tinwallsun::FT,
        Tinwallshd::FT,
        Tceiling::FT,
        Tinground::FT,
        Tintmass::FT,
        Twindow::FT
    ) where {FT<:AbstractFloat}

Calculate sensible heat flux from building surfaces to interior air.

# Arguments
- `Tbin`: Building interior air temperature [K]
- `Tinwallsun`: Interior sunlit wall temperature [K]
- `Tinwallshd`: Interior shaded wall temperature [K]
- `Tceiling`: Ceiling temperature [K]
- `Tinground`: Interior ground temperature [K]
- `Tintmass`: Internal mass temperature [K]
- `Twindow`: Window temperature [K]

# Returns
- `HbinWallSun::FT`: Sensible heat flux from sunlit wall [W/m²]
- `HbinWallshd::FT`: Sensible heat flux from shaded wall [W/m²]
- `HBinRoof::FT`: Sensible heat flux from roof [W/m²]
- `HBinGround::FT`: Sensible heat flux from ground [W/m²]
- `HbinIntMass::FT`: Sensible heat flux from internal mass [W/m²]
- `HbinWindow::FT`: Sensible heat flux from window [W/m²]
"""
function sensible_heat_flux_building_interior(
    Tbin::FT,
    Tinwallsun::FT,
    Tinwallshd::FT,
    Tceiling::FT,
    Tinground::FT,
    Tintmass::FT,
    Twindow::FT,
) where {FT<:AbstractFloat}
    # Convective heat transfer coefficients (Bueno et al. 2012, Oleson et al. 2019)
    hvert = FT(3.076)    # Vertical surfaces (walls) [W/m²K]
    hhorRed = FT(0.948)  # Horizontal with reduced convection [W/m²K]
    hhorEnh = FT(4.040)  # Horizontal with enhanced convection [W/m²K]

    # Calculate sensible heat fluxes
    HbinWallSun = hvert * (Tinwallsun - Tbin)
    HbinWallshd = hvert * (Tinwallshd - Tbin)
    HbinWindow = hvert * (Twindow - Tbin)
    HbinIntMass = hvert * (Tintmass - Tbin)

    # Roof and ground fluxes depend on temperature difference direction
    HBinRoof = ifelse(Tceiling > Tbin, hhorRed, hhorEnh) * (Tceiling - Tbin)
    HBinGround = ifelse(Tinground > Tbin, hhorEnh, hhorRed) * (Tinground - Tbin)

    return HbinWallSun, HbinWallshd, HBinRoof, HBinGround, HbinIntMass, HbinWindow
end
