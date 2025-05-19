"""
    volume_correction(V, EvL_Zs, RfH_Zs, RfL_Zs, EG, T_H, T_L, Lk)

Volume correction for negative values.

# Arguments
- `V`: Volume array
- `EvL_Zs`: Evaporation layers array
- `RfH_Zs`: Root fraction high vegetation layers array
- `RfL_Zs`: Root fraction low vegetation layers array
- `EG`: Ground evaporation
- `T_H`: High vegetation transpiration array
- `T_L`: Low vegetation transpiration array
- `Lk`: Leakage value

# Returns
- `V`: Corrected volume array
- `T_H`: Updated high vegetation transpiration
- `T_L`: Updated low vegetation transpiration
- `EG`: Updated ground evaporation
- `Lk`: Updated leakage
"""
function volume_correction(
    V::Vector{FT},
    EvL_Zs::Vector{FT},
    RfH_Zs::Matrix{FT},
    RfL_Zs::Matrix{FT},
    EG::FT,
    T_H::Vector{FT},
    T_L::Vector{FT},
    Lk::FT,
) where {FT<:AbstractFloat}

    # TODO: refactor functions to avoid code duplication
    V = copy(V)
    T_H = copy(T_H)
    T_L = copy(T_L)

    # Number of active layers
    lay_H = vec(count(>(0), RfH_Zs; dims=2))  # for Transp H
    lay_L = vec(count(>(0), RfL_Zs; dims=2))  # for Transp L
    lay_G = count(>(0), EvL_Zs)  # for Evaporation

    # Compensatory Mechanism on deeper layers
    for i in 1:lay_G
        if V[i] < 0 && i < lay_G
            V[i + 1] += V[i]
            V[i] = 0
        end
    end

    for i in 1:maximum(lay_L)
        if V[i] < 0 && i < maximum(lay_L)
            V[i + 1] += V[i]
            V[i] = 0
        end
    end

    for i in 1:maximum(lay_H)
        if V[i] < 0 && i < maximum(lay_H)
            V[i + 1] += V[i]
            V[i] = 0
        end
    end

    # Compensatory Mechanism on shallow layers
    if any(<(0), V)
        for i in findlast(<(0), V):-1:2
            if V[i] < 0 && V[i - 1] >= 0
                V[i - 1] += V[i]
                V[i] = 0
            end
        end
    end

    # Brutal Correction
    if any(<(0), V)
        # Ground evaporation correction
        for i in 1:lay_G
            if V[i] < 0
                if abs(V[i]) > EG
                    V[i] += EG
                    EG = 0
                else
                    EG += V[i]
                    V[i] = 0
                end
            end
        end

        # Low vegetation transpiration correction
        for i in 1:first(lay_L)
            if V[i] < 0
                if abs(V[i]) > sum(T_L)
                    V[i] += sum(T_L)
                    T_L .*= 0
                else
                    T_L .+= V[i] * T_L / sum(T_L)
                    V[i] = 0
                end
            end
        end

        # High vegetation transpiration correction
        for i in 1:first(lay_H)
            if V[i] < 0
                if abs(V[i]) > sum(T_H)
                    V[i] += sum(T_H)
                    T_H .*= 0
                else
                    T_H .+= V[i] * T_H / sum(T_H)
                    V[i] = 0
                end
            end
        end

        # Leakage correction
        if V[end] < 0
            for i in length(V):-1:2
                if V[i] < 0
                    Lk += V[i]
                    V[i] = 0
                end
            end
        end

        if V[end] < 0
            Lk += V[1]
            V[1] = 0
        end
    end

    # Final correction
    if any(<(0), V)
        EG += sum(V[V .< 0])
        V[V .< 0] .= 0
    end

    return V, T_H, T_L, EG, Lk
end
