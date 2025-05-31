function getSpeciesStats(
    species_group::SubDataFrame,
    expected_log2fc::Real
)
    getSpeciesStats(
        species_group[!,:log2_diff],
        species_group[!,:max_cv],
        species_group[!,:p_adj],
        expected_log2fc
    )
end

function getSpeciesStats(
    log2_diff::AbstractVector{Union{Missing, Float32}},
    max_cv::AbstractVector{Union{Missing, Float32}},
    p_adj::AbstractVector{Float32},
    expected_log2fc::Real)::@NamedTuple{
        total::Int64,
        expected_log2fc::Float32,
        median_log2fc::Float32,
        log2fc_bias::Float32,
        log2fc_std::Float32,
        iqr_outliers_at_cv20::Int64,
        cv20::Int64,
        p_adj_05::Int64,
        p_adj_05_at_fc_threshold::Int64
    }

    p_adj_05 = p_adj.<=0.05
    large_fc = abs.(log2_diff) .> abs(log2(2/3))
    
    log2_diff_filt = log2_diff[max_cv.<=20.0]
    _iqr_ = iqr(log2_diff_filt)
    Q1, Q3 = quantile(log2_diff_filt, 1/4), quantile(log2_diff_filt, 3/4)

    return (
        total = length(log2_diff),
        expected_log2fc = expected_log2fc,
        median_log2fc = median(log2_diff),
        log2fc_bias = median(log2_diff) - expected_log2fc,
        log2fc_std = std(log2_diff),
        iqr_outliers_at_cv20 =sum(log2_diff_filt .< (Q1 - 1.5*_iqr_)) + sum(log2_diff_filt .> (Q3 + 1.5*_iqr_)),
        cv20 = sum(max_cv.<=20.0),
        p_adj_05 = sum(p_adj_05),
        p_adj_05_at_fc_threshold = sum(p_adj_05 .& large_fc)
    )
end


function summarizeExperiment(
    out_dir::String,
    out_prefix::String,
    precursors::Bool,
    experiment::String,
    experiment_df::SubDataFrame,
    groupby_cols::Vector{Symbol},
    experiment_to_species_ratios::Dict{String, Dict{String, Float32}};
    min_n::UInt8 = UInt8(3)
)

    #E5H50Y45_vs_E45H50Y5 => "E5H50Y45", "vs", "E45H50Y5"
    species_ratios = experiment_to_species_ratios[experiment]
    group_a, group_b = split(experiment, "vs")
    group_a = group_a[1:end-1]
    group_b = group_b[2:end]
    CSV.write("/Users/n.t.wamsley/Desktop/experimenttest.csv", experiment_df)
    grouped_experiment_df = groupby(experiment_df, groupby_cols)
    exp_summary = combine(
        exp_group -> summarizeExperimentGroup(exp_group, group_a, group_b),
        grouped_experiment_df
    )
    CSV.write("/Users/n.t.wamsley/Desktop/tempdf.csv", exp_summary)
    filter!(x->coalesce(x.min_n, -1)>=min_n, exp_summary)
    filter!(x->!ismissing(x.p_value), exp_summary)
    exp_summary[!,:p_adj] = adjust(coalesce.(exp_summary[!,:p_value], 1.0f0), BenjaminiHochberg())
    CSV.write(joinpath(out_dir, "tables", experiment*"_"*out_prefix*"_log2diff.csv"), exp_summary)

    plot_name = "Precursors"
    if precursors==false
        plot_name = "Protein Groups"
    end

    p = plotThreeProteome(
        exp_summary,
        species_ratios,
        plot_name;
        cv_threshold = 20.0
    )
    savefig(p, joinpath(out_dir, "plots", experiment*"_"*out_prefix*"_log2diff.pdf"))
    species_grouped = groupby(
        exp_summary,:Species
    )
    
    experiment_stats = combine(species_grouped) do sdf
        expected_log2_fc = log2(species_ratios[sdf.Species[1]])
        getSpeciesStats(
            sdf, 
            expected_log2_fc
            )
    end 
    experiment_stats[!,:Experiment] .= experiment

    return experiment_stats
end

function summarizeExperimentGroup(
    experiment_group::SubDataFrame,
    cond_a::AbstractString, 
    cond_b::AbstractString
)
    if size(experiment_group, 1) > 2
        throw("Experiment group has more than three entries. Indicates error.")
    end

    try 
    summarizeExperimentGroup(
        experiment_group[!,:Condition],
        cond_a, cond_b,
        experiment_group[!,:mean],
        experiment_group[!,:var],
        experiment_group[!,:CV],
        experiment_group[!,:n],
        experiment_group[!,:max_qval]
    )
    catch
        println("typeof( experiment_group[!,:Condition]): ", typeof( experiment_group[!,:Condition]))
        println("typeof( experiment_group[!,:mean]): ", typeof(experiment_group[!,:mean]))
        println("typeof( experiment_group[!,:var]): ", typeof( experiment_group[!,:var]))
        println("typeof( experiment_group[!,:CV]): ", typeof( experiment_group[!,:CV]))
        println("typeof( experiment_group[!,:n]): ", typeof( experiment_group[!,:n]))
        println("typeof( experiment_group[!,:max_qval]): ", typeof( experiment_group[!,:max_qval]))
        println("cond_a $cond_a cond_b $cond_b")

        throw("darn")
    end
end


function summarizeExperimentGroup(
    conditions::AbstractArray{S},
    cond_a::AbstractString,
    cond_b::AbstractString,
    mean_abundance::AbstractVector{Float32},
    group_var::AbstractVector{Float32},
    cv::AbstractVector{Float32},
    n::AbstractVector{UInt8},
    max_qval::AbstractVector{Float32}
)::@NamedTuple{log2_diff::Union{Missing, Float32},
              log2_mean::Union{Missing, Float32},
              max_cv::Union{Missing, Float32},
              min_n::Union{Missing, UInt8},
              max_qval::Union{Missing, Float32},
              p_value::Union{Missing, Float32}} where {S<:AbstractString}
    cond_a = findfirst(x->x==cond_a, conditions)
    cond_b = findfirst(x->x==cond_b, conditions)
    if isnothing(cond_b) | isnothing(cond_a)
        return (
            log2_diff = missing, log2_mean = missing, max_cv = missing, min_n = missing, max_qval = missing, p_value = missing
        )
    end

    min_n = min(n[cond_a], n[cond_b])
    if min_n > 2
    return (
        log2_diff = log2(mean_abundance[cond_b]) - log2(mean_abundance[cond_a]),
        log2_mean = log2((mean_abundance[cond_a] + mean_abundance[cond_b])/2),
        max_cv = max(cv[cond_a], cv[cond_b]),
        min_n = min_n,
        max_qval = max(max_qval[cond_a], max_qval[cond_b]),
        p_value = pvalue(EqualVarianceTTest(
            Int64(n[cond_a]), Int64(n[cond_b]),
            mean_abundance[cond_a], mean_abundance[cond_b],
            group_var[cond_a], group_var[cond_b]
        ))
    )
    else
        return (
            log2_diff = log2(mean_abundance[cond_b]) - log2(mean_abundance[cond_a]),
            log2_mean = log2((mean_abundance[cond_a] + mean_abundance[cond_b])/2),
            max_cv = max(cv[cond_a], cv[cond_b]),
            min_n = min_n,
            max_qval = max(max_qval[cond_a], max_qval[cond_b]),
            p_value = missing
        )
    end 
end