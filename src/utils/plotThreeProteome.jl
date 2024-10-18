function plotThreeProteome(
    df::DataFrame,
    species_to_log2fc::Dict{String, Float32},
    plot_name::String;
    cv_threshold::Real = 20.0f0)
    
    dfc = copy(df)
    filter!(x->x.max_cv<=cv_threshold, dfc)
    n = size(dfc, 1)
    gdf = groupby(dfc, :Species)
    nt = NamedTuple.(keys(gdf))
    p = plot(legend=:outertopright, 
            show = true, 
            title = "$n $plot_name \n at 1% FDR and CV<=20%",
            topmargin=5mm, 
            dpi = 300)
    i = 1
    for (k,v) in pairs(gdf)
        density!(p, gdf[k][:,:log2_diff], label=nothing, color = i, bins = LinRange(-3, 3, 100), xlim = (-3.5, 3.5), ylim = (0, 4.0), show = true, normalize = :probability, alpha = 1.0)
        vline!(p, [log2(species_to_log2fc[nt[i][:Species]])], color = i, label = label="$(nt[i][:Species])", legend = false, lw = 2)
        vline!(p, [median(gdf[k][:,:log2_diff])], color = i, linestyle = :dash, label = label="$(nt[i][:Species])", legend = false)
        i += 1
    end
    return p
end