"""
    evaporation_layers(Zs, Zdes)

Calculate the evaporation layer fraction for each soil layer based on the desorption depth.

# Arguments
- `Zs::Vector{FT}`: Vector of depth layers in mm, from top to bottom [1....m+1]
- `Zdes::FT`: Depth of desorption in mm

# Returns
- `EvL_Zs::Vector{FT}`: Vector of evaporation layer fractions [1...m], where each element
  represents the fraction of evaporation for the corresponding layer

# Notes
- The function checks if the desorption depth is within the valid range defined by the layers
- If Zdes < first layer depth: returns zeros with error message "ERROR FIRST LAYER TOO DEPTH"
- If Zdes > last layer depth: returns zeros with error message "ERROR LAST LAYER TOO SHALLOW"
- The sum of all fractions equals 1 when Zdes is within valid range
"""
function evaporation_layers(Zs::Vector{FT}, Zdes::FT) where {FT<:AbstractFloat}
    n = length(Zs) - 1
    EvL_Zs = zeros(n)

    if Zdes < Zs[1]
        println("ERROR FIRST LAYER TOO DEEP")
        return EvL_Zs
    end

    if Zdes > Zs[n + 1]
        println("ERROR LAST LAYER TOO SHALLOW")
        return EvL_Zs
    end

    NiL = sum(Zs .< Zdes)  # Number interested Layer
    dz = diff(Zs)  # [mm] Thickness of the Layers
    dz[NiL] = Zdes - Zs[NiL]  # Thickness last layer
    dz = dz[1:NiL]  # Interested Layer
    EvL_Zs[1:NiL] = dz / Zdes

    return EvL_Zs
end
