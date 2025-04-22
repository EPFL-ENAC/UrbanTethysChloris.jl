"""
    view_factors_geometry(H::FT, W::FT, a::FT, ht::FT, d::FT, person::PersonParameters{FT},
                         option_surface::Int, mc_sample_size::Int, n_rays::Int) where {FT<:AbstractFloat}

Compute view factors for urban canyon geometry.

# Arguments
- `H`: canyon height [m]
- `W`: canyon width [m]
- `a`: normalized tree radius [-]
- `ht`: normalized tree height [-]
- `d`: normalized tree distance from wall [-]
- `person`: PersonParameters object containing position information
- `option_surface`: specifies emitting surface (1=wall1, 2=wall2, 3=ground, 4=tree1, 5=tree2, 6=sky, 7=point)
- `mc_sample_size`: number of emitting points per surface
- `n_rays`: number of rays emitted per point

# Returns
Tuple of view factors (VG, VW1, VW2, VS, VT1, VT2) for ground, walls, sky and trees
"""
function view_factors_geometry(
    H::FT,
    W::FT,
    a::FT,
    ht::FT,
    d::FT,
    person::ModelComponents.Parameters.PersonParameters{FT},
    option_surface::Int,
    mc_sample_size::Int,
    n_rays::Int,
) where {FT<:AbstractFloat}
    # Normalize dimensions
    h = H/W
    w = W/W
    pz = person.PositionPz/W
    px = person.PositionPx/W

    # Define geometry coordinates
    x2 = FT[1, 1 + w] # Ground
    z2 = FT[0, 0]

    x3 = FT[1, 1] # Wall 1
    z3 = FT[h, 0]

    x4 = FT[1 + w, 1 + w] # Wall 2
    z4 = FT[0, h]

    x5 = FT[1, 1 + w] # Sky
    z5 = FT[h, h]

    # Tree positions
    r = a*w
    if r > 0
        xc = 1 + d*w # Tree 1 center
        yc = ht*w
        xc2 = 1 + w - d*w # Tree 2 center
    else
        xc = xc2 = yc = zero(FT)
    end

    # Monte Carlo parameters
    rand_sz = rand(FT, mc_sample_size)
    delta_rays = LinRange(zero(FT), one(FT), n_rays÷2 + 1)
    angle_dist = asin.(delta_rays)
    ray_angle_q1 = reverse(π/2 .- angle_dist)
    ray_angle_q2 = π/2 .+ angle_dist
    ray_angle = [ray_angle_q1[1:(end - 1)]; ray_angle_q2]

    stc = FT(1e-10)  # Offset from surface

    # Generate emission points and angles based on surface
    if option_surface == 1 # Wall 1
        XSv = fill(1+stc, mc_sample_size)
        YSv = h .* rand_sz
        dthe = ray_angle .- π/2
    elseif option_surface == 2 # Wall 2
        XSv = fill(1+w-stc, mc_sample_size)
        YSv = h .* rand_sz
        dthe = ray_angle .+ π/2
    elseif option_surface == 3 # Ground
        XSv = 1 .+ w .* rand_sz
        YSv = fill(stc, mc_sample_size)
        dthe = ray_angle
    elseif option_surface == 4 # Tree 1
        ang = 2π .* rand_sz
        xt = (r+stc) .* cos.(ang)
        yt = (r+stc) .* sin.(ang)
        XSv = xc .+ xt
        YSv = yc .+ yt

        if r == 0
            XSv = [0.0]
            YSv = [0.0]
        end

        dthe = (ones(mc_sample_size) .* (ray_angle .- π/2)') .+ ang
    elseif option_surface == 5 # Tree 2
        ang = 2π .* rand_sz
        xt = (r+stc) .* cos.(ang)
        yt = (r+stc) .* sin.(ang)
        XSv = xc2 .+ xt
        YSv = yc .+ yt

        if r == 0
            XSv = [0.0]
            YSv = [0.0]
        end

        dthe = (ones(mc_sample_size) .* (ray_angle .- π/2)') .+ ang
    elseif option_surface == 6 # Sky
        XSv = 1 .+ w .* rand_sz
        YSv = fill(h-stc, mc_sample_size)
        dthe = ray_angle .+ π
    else # Point
        rp6 = FT(1/1000)
        ang = 2π .* rand_sz
        xp6 = (rp6+stc) .* cos.(ang)
        yp6 = (rp6+stc) .* sin.(ang)
        XSv = 1+px .+ xp6
        YSv = pz .+ yp6
        dthe = (ones(mc_sample_size) .* (ray_angle .- π/2)') .+ ang
    end

    # Parameters for ray tracing
    dmax = sqrt(h^2 + w^2) + sqrt(h^2 + w^2)/100
    sz = w/1000

    return view_factors_computation(
        XSv, YSv, dmax, sz, dthe, x2, z2, x3, z3, x4, z4, xc, yc, r, xc2, x5, z5
    )
end
