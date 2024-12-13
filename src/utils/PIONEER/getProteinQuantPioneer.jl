function getProteinQuantPioneer(
    pioneer_results_dir::String,
    run_to_condition::Dict{String, String},
    condition_to_experiment::Dict{String, String}
)

    #Load dataframe
    pioneer_pg = DataFrame(Tables.columntable(Arrow.Table(
        joinpath(
                pioneer_results_dir, 
                "protein_groups_long.arrow"
                 )
    )))
    #Check for columns and rename to harmonize with DIANN results
    column_names = Set(Symbol.(names(pioneer_pg)))
    checkHasColumns(
        [:file_name,:species,:protein,:abundance,:target,:n_peptides],
        column_names
    )
    cols_to_keep = [:file_name, :species, :protein, :abundance]
    select!(pioneer_pg, cols_to_keep)
    rename!(pioneer_pg, [:file_name,:species,:protein,:abundance] .=> [:Run, :Species, :ProteinGroup, :PGMaxLFQ])
    pioneer_pg[!,:Condition] = [run_to_condition[run_name] for run_name in pioneer_pg[!,:Run]]
    pioneer_pg[!,:Experiment] = [condition_to_experiment[condition] for condition in pioneer_pg[!,:Condition]]
    pioneer_pg[!,:PGQValue] .= zero(Float32)
    #filter!(x->x.Run!="E5H50Y45_3",pioneer_pg)
    return pioneer_pg
end