#=
function getTaskRanges(N::Int, M::Int)
    # Calculate size of each chunk (rounded up to include all elements)
    chunk_size = ceil(Int, N/M)
    
    # Create the ranges
    ranges = [((i-1)*chunk_size + 1):min(i*chunk_size, N) for i in 1:12]
    
    return ranges
end


function analyzeMMCC(
    analyte_df::DataFrame,
    id_col::Symbol,
    quant_col::Symbol;
    n_bootstraps = 1000,
    power = 0.5,
    min_t = 0.5
    )
    sort!(analyte_df, [id_col,:Condition])
    analyte_gdf = groupby(analyte_df, id_col)
    task_ranges = getTaskRanges(length(analyte_gdf), Threads.nthreads())
    lloq_curve_dfs = Vector{Any}(undef, Threads.nthreads())
    @time begin
        Threads.@threads for (idx, task_range) in collect(enumerate(task_ranges))
            lloq_curve_dfs[idx] = BootstrapGDF(analyte_gdf[task_range], 
            id_col,
            quant_col, 
            power = power,
            min_t = min_t,
            n_bootstraps = 1000,
            n_xbins = 100)
        end
    end
    lloq_curve_df = vcat(lloq_curve_dfs...)
    filter!(x->coalesce(x.lloq, typemax(Float32))<100.0, lloq_curve_df)
    return lloq_curve_df
end
=#