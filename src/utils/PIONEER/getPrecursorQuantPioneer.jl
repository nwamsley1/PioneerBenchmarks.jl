function containsForbiddenMods(
    modified_sequence::String,
    forbidden_mods::Vector{String})
    for mod in forbidden_mods
        if occursin(mod, modified_sequence)
            return true
        end
    end
    return false
end

function getPrecursorQuantPioneer(
    pioneer_results_dir::String,
    run_to_condition::Dict{String, String},
    condition_to_experiment::Dict{String, String};
    forbidden_mods::Vector{String} = ["UNIMOD:35"]
)
    #Load dataframe
    pioneer_pr = DataFrame(Tables.columntable(Arrow.Table(
        joinpath(
                pioneer_results_dir, 
                "precursors_long.arrow"
                 )
    )))

    #Are the expected columns present?
    #Remove unneccessary columns
    column_names = Set(Symbol.(names(pioneer_pr)))
    cols_to_keep = [
    :file_name, :species, :sequence, :structural_mods, 
    :charge, :peak_area
    ]
    checkHasColumns(cols_to_keep, column_names)
    select!(pioneer_pr, cols_to_keep)
    #Standardize sequence names to agree with DIA-NN 
    pioneer_pr[!,:modified_sequence] = getModifiedSequence.(pioneer_pr[!,:sequence], "", pioneer_pr[!,:structural_mods])
    pioneer_pr[!,:modified_sequence]  = pioneer_pr[!,:modified_sequence].*string.(pioneer_pr[!,:charge])

    #Standardize column names 
    select!(pioneer_pr, [:file_name,:species,:modified_sequence,:peak_area])
    rename!(pioneer_pr, [:file_name,:species,:modified_sequence,:peak_area] .=> [:Run, :Species, :PrecursorId,:PrecursorQuantity])
    #Parse file names to get the condition and experimemnt for each precursor id 
    pioneer_pr[!,:Condition] = [run_to_condition[run_name] for run_name in pioneer_pr[!,:Run]]
    pioneer_pr[!,:Experiment] .= ""
    for i in range(1, size(pioneer_pr, 1))
        if haskey(condition_to_experiment, pioneer_pr[i,:Condition])
            exp = condition_to_experiment[pioneer_pr[i,:Condition]]
            pioneer_pr[i,:Experiment] = exp
        end
    end
    filter!(x->x.Experiment.!="", pioneer_pr)
    #pioneer_pr[!,:Experiment] = [condition_to_experiment[condition] for condition in pioneer_pr[!,:Condition]]
    #QValue already filtered at 1% by Pioneer. Add columns to harmonize with diann output 
    pioneer_pr[!,:QValue] .= zero(Float32)
    #filter!(x->x.Run!="E5H50Y45_3",pioneer_pr)
    #Reorder columns
    pioneer_pr = pioneer_pr[!,[:Run,:Condition,:Experiment,:Species,:PrecursorId,:PrecursorQuantity,:QValue]]

    #Remove missing values and change columns eltypes from Union{Missing, Float32} => Float32
    filter!(x->!ismissing(x.PrecursorQuantity), pioneer_pr)
    filter!(x->!ismissing(x.Condition), pioneer_pr)
    pioneer_pr[!,:PrecursorQuantity]= Float32.(pioneer_pr[!,:PrecursorQuantity])
    pioneer_pr[!,:Condition] = string.(pioneer_pr[!,:Condition])
    filter!(x->!containsForbiddenMods(x.PrecursorId, forbidden_mods), pioneer_pr)
    return pioneer_pr
end