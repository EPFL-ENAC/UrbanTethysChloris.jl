"""
    view_factors_canyon(
        geometry::ModelComponents.Parameters.UrbanGeometryParameters{FT},
        person::ModelComponents.Parameters.PersonParameters{FT},
        mc_sample_size::Int = 1000,
        n_rays::Int = 200
    ) where {FT<:AbstractFloat}

Calculate view factors for an urban street canyon, considering both cases with and without trees.

# Arguments
- `geometry`: Urban geometry parameters including canyon dimensions and tree properties
- `person`: Parameters describing a person's position in the canyon
- `mc_sample_size`: Number of Monte Carlo samples for ray tracing (default: 1000)
- `n_rays`: Number of rays to emit per sample point (default: 200)

# Returns
A tuple containing:
- `vf`: ViewFactor object containing both tree and no-tree view factors
- `vf_point`: Point-specific view factors
"""
function view_factors_canyon(
    geometry::ModelComponents.Parameters.UrbanGeometryParameters{FT},
    person::ModelComponents.Parameters.PersonParameters{FT},
    mc_sample_size::Int=1000,
    n_rays::Int=200,
) where {FT<:AbstractFloat}
    if geometry.trees
        radius_tree = geometry.radius_tree
        htree = geometry.htree
        distance_tree = geometry.distance_tree
    else
        radius_tree = zero(FT), htree = - geometry.Height_canyon / 10000
        distance_tree = zero(FT)
    end

    # Compute view factors with ray tracing and reciprocity
    vf_trees, vf_point, _ = view_factors_ray_tracing_reciprocity(
        geometry.Height_canyon,
        geometry.Width_canyon,
        radius_tree,
        htree,
        distance_tree,
        person,
        mc_sample_size,
        n_rays,
    )

    # Calculate view factors without trees using analytical solution
    vf_no_trees = view_factors_analytical(geometry.Height_canyon, geometry.Width_canyon)

    vf = ViewFactor(vf_no_trees, vf_trees)

    return vf, vf_point
end
