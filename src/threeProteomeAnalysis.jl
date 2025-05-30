function threeProteomeAnalysis(
    data_path::String,
    fasta_dir::String,
    run_to_condition_path::String,
    key_file_path::String,
    data_type::String,
    out_dir::String)
    if !isdir(out_dir)
        mkdir(out_dir)
    end
    if !isdir(joinpath(out_dir, "plots"))
        mkdir(joinpath(out_dir, "plots"))
    end
    if !isdir(joinpath(out_dir, "tables"))
        mkdir(joinpath(out_dir, "tables"))
    end
    if !isdir(joinpath(out_dir, "stats"))
        mkdir(joinpath(out_dir, "stats"))
    end

    if !any([occursin(data_type, _type_) for _type_ in ["diann","pioneer"]])
        throw("Unrecognized data type $data_type. Valid options are ['diann', 'pioneer']")
    end

    if !isfile(key_file_path)
        throw("Experiment key file path: $key_file_path \n was 
        invalid. File does not exist")
    end

    if !isfile(run_to_condition_path)
        throw("run_to_condition_path: $run_to_condition_path \n was 
        invalid. File does not exist")
    end
    run_to_condition = readRunToConditionTable(
        run_to_condition_path
    )
    #Parse file key
    #run_to_condition, condition_to_experiment, experiment_to_species_ratios = conditionToExperiment(key_file_path)
    condition_to_experiment, experiment_to_species_ratios = conditionToExperiment(key_file_path)
    protein_groups_table, precursors_table = nothing, nothing
    if data_type=="diann"
        if !isfile(data_path)
            throw("DIANN results file not found. Path supplies was: $data_path")
        end
        if !isdir(fasta_dir)
            throw("fasta directory not found. Supplied : $fasta_dir")
        end
        diann_r = LoadDiannResults(
            data_path,
            fasta_dir,
            run_to_condition,
            condition_to_experiment
        )
        protein_groups_table = getProteinQuantDIANN(diann_r)
        precursors_table = getPrecursorQuantDIANN(diann_r)
    elseif data_type=="pioneer"
        precursors_table = getPrecursorQuantPioneer(
            data_path,
            run_to_condition,
            condition_to_experiment
        )
        protein_groups_table = getProteinQuantPioneer(
            data_path,
            run_to_condition,
            condition_to_experiment
        )
    else
        throw("Unrecognized data type $data_type. Valid 
        options are ['diann', 'pioneer']")
    end
        
    ##########
    #Precursor Summary 
    pr_condition_summary = combine(
        condition->summarizePrecursorCondition(condition),    
        groupby(precursors_table,[:Experiment,:Condition,:Species,:PrecursorId])
    );
    species_order = Dict("HUMAN" => 1, "YEAST" => 2, "ECOLI" => 3)
    #println("head pr_condition_summary: ", first(pr_condition_summary, 5))
    #println("unique(pr_condition_summary.Species): ", unique(pr_condition_summary.Species))
    filter!(x->!occursin("CONTAMINANT", x.Species), pr_condition_summary)
    filter!(x->!occursin(";", x.Species), pr_condition_summary)
    filter!(x->!occursin("CONTAMINANT", x.Species), protein_groups_table)
    filter!(x->!occursin(";", x.Species), protein_groups_table)
    println("protein_groups_table", first(protein_groups_table, 5))  
    #println("unique(pr_condition_summary.Species): ", unique(pr_condition_summary.Species))
    #return
    # Sort the DataFrame using the custom order
    sort!(pr_condition_summary, :Species, by = x -> species_order[x])
    CSV.write("/Users/n.t.wamsley/Desktop/pr_condition.csv", pr_condition_summary)
    precursor_group_stats = combine(groupby(pr_condition_summary,:Experiment)) do experiment_df
        println("size(experiment_df) ", size(experiment_df))
        summarizeExperiment(
                out_dir,
                "precursors",
                true,
                experiment_df.Experiment[1], 
                experiment_df,
                [:Species,:PrecursorId],
                experiment_to_species_ratios,
                min_n = UInt8(3)
                )
    end
    sort!(precursor_group_stats,[:Experiment,:Species])
    CSV.write(joinpath(out_dir, "stats", "precursors_stats.csv"), precursor_group_stats)

    ##########
    #Protein Groups Summary
    CSV.write(joinpath(out_dir, "stats", "protein_groups_table.csv"), protein_groups_table)

    pg_condition_summary = combine(
        condition->summarizeProteinGroupCondition(condition),    
        groupby(protein_groups_table,[:Experiment,:Condition,:Species,:ProteinGroup])
    );
    species_order = Dict("HUMAN" => 1, "YEAST" => 2, "ECOLI" => 3)
    # Sort the DataFrame using the custom order
    sort!(pg_condition_summary, :Species, by = x -> species_order[x])
    protein_group_stats = combine(groupby(pg_condition_summary,:Experiment)) do experiment_df
        summarizeExperiment(
                out_dir,
                "protein_groups",
                false,
                experiment_df.Experiment[1], 
                experiment_df,
                [:Species,:ProteinGroup],
                experiment_to_species_ratios,
                min_n = UInt8(3)
                )
    end
    sort!(protein_group_stats,[:Experiment,:Species])
    CSV.write(joinpath(out_dir, "stats", "protein_groups_stats.csv"), protein_group_stats)
    return nothing
end