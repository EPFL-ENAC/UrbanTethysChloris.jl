"""
    root_fraction_general(
        Zs::Vector{FT},
        CASE_ROOT::Int,
        ZR95_H::Vector{FT},
        ZR50_H::Vector{FT},
        ZR95_L::Vector{FT},
        ZR50_L::Vector{FT},
        ZRmax_H::Vector{FT},
        ZRmax_L::Vector{FT}
    ) where {FT<:AbstractFloat}

Compute root fraction distributions for high and low vegetation using different root profile models.

# Arguments
- `Zs`: Soil layer depths [mm]
- `CASE_ROOT`: Root profile type:
    1. Exponential (Arora and Boer 2005)
    2. Linear Dose Response (Schenk and Jackson 2002)
    3. Constant Profile
    4. Deep (Tap) Root Profile
- `ZR95_H`: 95th percentile root depth for high vegetation [mm]
- `ZR50_H`: 50th percentile root depth for high vegetation [mm]
- `ZR95_L`: 95th percentile root depth for low vegetation [mm]
- `ZR50_L`: 50th percentile root depth for low vegetation [mm]
- `ZRmax_H`: Maximum root depth for high vegetation [mm]
- `ZRmax_L`: Maximum root depth for low vegetation [mm]

# Returns
- `RfH_Zs`: Root fraction distribution for high vegetation [dimensionless]
- `RfL_Zs`: Root fraction distribution for low vegetation [dimensionless]

"""
function root_fraction_general(
    Zs::Vector{FT},
    CASE_ROOT::Int,
    ZR95_H::Vector{FT},
    ZR50_H::Vector{FT},
    ZR95_L::Vector{FT},
    ZR50_L::Vector{FT},
    ZRmax_H::Vector{FT},
    ZRmax_L::Vector{FT},
) where {FT<:AbstractFloat}
    # TODO: split function in two, as it could two calls to a height-dependent function
    n = length(Zs) - 1
    cc = length(ZR95_H)
    RfH_Zs = zeros(FT, cc, n)
    RfL_Zs = zeros(FT, cc, n)

    # Check root depth validity
    for j in 1:cc
        if ZR95_H[j] > Zs[n + 1] ||
            ZR95_L[j] > Zs[n + 1] ||
            ZRmax_H[j] > Zs[n + 1] ||
            ZRmax_L[j] > Zs[n + 1]
            @warn "ERROR: LAST LAYER TOO SHALLOW FOR ACCOMODATING ROOTS"
            return RfH_Zs, RfL_Zs
        end
    end

    if CASE_ROOT == 1  # Exponential Profile
        # Arora and Boer 2005
        eta_H = 3 ./ ZR95_H  # [1/mm] Shape of Root Distribution
        eta_L = 3 ./ ZR95_L

        for j in 1:cc
            i = 1
            if ZR95_H[j] ≠ 0.0
                while i <= n
                    if ZR95_H[j] > Zs[i + 1]
                        RfH_Zs[j, i] = exp(-eta_H[j] * Zs[i]) - exp(-eta_H[j] * Zs[i + 1])
                    else
                        RfH_Zs[j, i] = exp(-eta_H[j] * Zs[i]) - exp(-eta_H[j] * ZR95_H[j])
                        break
                    end
                    i += 1
                end
            end

            i = 1
            if ZR95_L[j] ≠ 0.0
                while i <= n
                    if ZR95_L[j] > Zs[i + 1]
                        RfL_Zs[j, i] = exp(-eta_L[j] * Zs[i]) - exp(-eta_L[j] * Zs[i + 1])
                    else
                        RfL_Zs[j, i] = exp(-eta_L[j] * Zs[i]) - exp(-eta_L[j] * ZR95_L[j])
                        break
                    end
                    i += 1
                end
            end

            # Root proportion normalization
            Rto1 = 0.9502
            RfH_Zs[j, :] ./= Rto1
            RfL_Zs[j, :] ./= Rto1
        end

    elseif CASE_ROOT == 2  # Linear Dose Response
        # Schenk and Jackson 2002, Collins and Bras 2007
        c_H = 2.94 ./ log.(ZR50_H ./ ZR95_H)
        c_L = 2.94 ./ log.(ZR50_L ./ ZR95_L)

        for j in 1:cc
            if ZR95_H[j] != 0
                i = 1
                while i <= n
                    if ZR95_H[j] > Zs[i + 1]
                        RfH_Zs[j, i] =
                            1 / (1 + (Zs[i + 1] / ZR50_H[j])^c_H[j]) -
                            1 / (1 + (Zs[i] / ZR50_H[j])^c_H[j])
                    else
                        RfH_Zs[j, i] =
                            1 / (1 + (ZR95_H[j] / ZR50_H[j])^c_H[j]) -
                            1 / (1 + (Zs[i] / ZR50_H[j])^c_H[j])
                        break
                    end
                    i += 1
                end
            end

            # Similar for low vegetation
            if ZR95_L[j] != 0
                i = 1
                while i <= n
                    if ZR95_L[j] > Zs[i + 1]
                        RfL_Zs[j, i] =
                            1 / (1 + (Zs[i + 1] / ZR50_L[j])^c_L[j]) -
                            1 / (1 + (Zs[i] / ZR50_L[j])^c_L[j])
                    else
                        RfL_Zs[j, i] =
                            1 / (1 + (ZR95_L[j] / ZR50_L[j])^c_L[j]) -
                            1 / (1 + (Zs[i] / ZR50_L[j])^c_L[j])
                        break
                    end
                    i += 1
                end
            end
        end

        Rto1 = 0.9498
        RfH_Zs[j, :] ./= Rto1
        RfL_Zs[j, :] ./= Rto1

    elseif CASE_ROOT == 3  # Constant Profile
        for j in 1:cc
            if ZR95_H[j] != 0
                i = 1
                while i <= n
                    if ZR95_H[j] > Zs[i + 1]
                        RfH_Zs[j, i] = (Zs[i + 1] - Zs[i]) / ZR95_H[j]
                    else
                        RfH_Zs[j, i] = (ZR95_H[j] - Zs[i]) / ZR95_H[j]
                        break
                    end
                    i += 1
                end
            end

            if ZR95_L[j] != 0
                i = 1
                while i <= n
                    if ZR95_L[j] > Zs[i + 1]
                        RfL_Zs[j, i] = (Zs[i + 1] - Zs[i]) / ZR95_L[j]
                    else
                        RfL_Zs[j, i] = (ZR95_L[j] - Zs[i]) / ZR95_L[j]
                        break
                    end
                    i += 1
                end
            end
        end

    elseif CASE_ROOT == 4  # Deep (Tap) Root Profile
        c_H = 2.94 ./ log.(ZR50_H ./ ZR95_H)
        c_L = 2.94 ./ log.(ZR50_L ./ ZR95_L)

        for j in 1:cc
            if ZR95_H[j] != 0
                i = 1
                while i <= n
                    if ZR95_H[j] > Zs[i + 1]
                        RfH_Zs[j, i] =
                            1 / (1 + (Zs[i + 1] / ZR50_H[j])^c_H[j]) -
                            1 / (1 + (Zs[i] / ZR50_H[j])^c_H[j])
                    elseif ZR95_H[j] <= Zs[i + 1] && ZR95_H[j] > Zs[i]
                        RfH_Zs[j, i] =
                            1 / (1 + (ZR95_H[j] / ZR50_H[j])^c_H[j]) -
                            1 / (1 + (Zs[i] / ZR50_H[j])^c_H[j])
                        if ZRmax_H[j] <= Zs[i + 1]
                            RfH_Zs[j, i] +=
                                0.0502 * (ZRmax_H[j] - ZR95_H[j]) / (ZRmax_H[j] - ZR95_H[j])
                            break
                        else
                            RfH_Zs[j, i] +=
                                0.0502 * (Zs[i + 1] - ZR95_H[j]) / (ZRmax_H[j] - ZR95_H[j])
                        end
                    elseif ZRmax_H[j] > Zs[i + 1]
                        RfH_Zs[j, i] =
                            0.0502 * (Zs[i + 1] - Zs[i]) / (ZRmax_H[j] - ZR95_H[j])
                    else
                        RfH_Zs[j, i] =
                            0.0502 * (ZRmax_H[j] - Zs[i]) / (ZRmax_H[j] - ZR95_H[j])
                        break
                    end
                    i += 1
                end
            end

            # Similar logic for low vegetation...
            if ZR95_L[j] != 0
                i = 1
                while i <= n
                    if ZR95_L[j] > Zs[i + 1]
                        RfL_Zs[j, i] =
                            1 / (1 + (Zs[i + 1] / ZR50_L[j])^c_L[j]) -
                            1 / (1 + (Zs[i] / ZR50_L[j])^c_L[j])
                    elseif ZR95_L[j] <= Zs[i + 1] && ZR95_L[j] > Zs[i]
                        RfL_Zs[j, i] =
                            1 / (1 + (ZR95_L[j] / ZR50_L[j])^c_L[j]) -
                            1 / (1 + (Zs[i] / ZR50_L[j])^c_L[j])
                        if ZRmax_L[j] <= Zs[i + 1]
                            RfL_Zs[j, i] +=
                                0.0502 * (ZRmax_L[j] - ZR95_L[j]) / (ZRmax_L[j] - ZR95_L[j])
                            break
                        else
                            RfL_Zs[j, i] +=
                                0.0502 * (Zs[i + 1] - ZR95_L[j]) / (ZRmax_L[j] - ZR95_L[j])
                        end
                    elseif ZRmax_L[j] > Zs[i + 1]
                        RfL_Zs[j, i] =
                            0.0502 * (Zs[i + 1] - Zs[i]) / (ZRmax_L[j] - ZR95_L[j])
                    else
                        RfL_Zs[j, i] =
                            0.0502 * (ZRmax_L[j] - Zs[i]) / (ZRmax_L[j] - ZR95_L[j])
                        break
                    end
                    i += 1
                end
            end
        end
    end

    return RfH_Zs, RfL_Zs
end
