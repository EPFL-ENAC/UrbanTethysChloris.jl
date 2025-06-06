"""
    urban_roughness(
        hc_H::FT,
        hc_L::FT,
        Csoil::Bool,
        Croad::Bool,
        Croof::Bool
    ) where {FT<:AbstractFloat}

Calculate urban surface roughness parameters.

# Arguments
- `hc_H`: Canopy height of high vegetation layer [m]
- `hc_L`: Canopy height of low vegetation layer [m]
- `Csoil`: Presence of soil
- `Croad`: Presence of road
- `Croof`: Presence of roof

# Returns
- `zom`: Roughness momentum [m]
- `zoh`: Roughness heat [m]
- `zom_ground`: Ground roughness momentum [m]
- `zoh_ground`: Ground roughness heat [m]
- `disp_h`: Maximum displacement height [m]
- `zom_H`: High vegetation roughness momentum [m]
- `zom_L`: Low vegetation roughness momentum [m]
- `zoh_H`: High vegetation roughness heat [m]
- `zoh_L`: Low vegetation roughness heat [m]
- `d_H`: High vegetation displacement height [m]
- `d_L`: Low vegetation displacement height [m]
- `zom_other`: Other urban surfaces roughness momentum [m]
"""
function urban_roughness(
    hc_H::FT, hc_L::FT, Csoil::Bool, Croad::Bool, Croof::Bool
) where {FT<:AbstractFloat}
    # Roughness parameters
    zom_soil = FT(0.003) * FT(Csoil)        # Bare soil roughness momentum [m]
    zom_road = FT(0.003) * FT(Croad)        # Road roughness momentum [m]
    zom_roof = FT(0.01) * FT(Croof)         # Roof roughness momentum [m]
    zom_H = FT(0.123) * hc_H                # High vegetation roughness momentum [m]
    zom_L = FT(0.123) * hc_L                # Low vegetation roughness momentum [m]

    # Other surface roughness
    zom_other = max(zom_soil, zom_road, zom_roof)

    # Heat roughness parameters
    zoh_L = FT(0.1) * zom_L
    zoh_H = FT(0.1) * zom_H
    zoh_other = FT(0.1) * zom_other

    # Patch scale roughness
    zom = max(zom_H, zom_L, zom_other)
    zoh = max(zoh_H, zoh_L, zoh_other)
    zom_ground = max(zom_L, zom_other)
    zoh_ground = max(zoh_L, zoh_other)

    # Displacement heights
    d_L = FT(2/3) * hc_L
    d_H = FT(2/3) * hc_H
    disp_h = max(d_H, d_L)

    return (
        zom,
        zoh,
        zom_ground,
        zoh_ground,
        disp_h,
        zom_H,
        zom_L,
        zoh_H,
        zoh_L,
        d_H,
        d_L,
        zom_other,
    )
end
