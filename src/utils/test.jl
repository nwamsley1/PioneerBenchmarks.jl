


test_template = first(getPrecursorQuantDIANN(diann_r), 4)

out_dir = "/Users/n.t.wamsley/Desktop/PIONEER_AltimeterOlsenExploris200ng_M0M1M2"
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


##########
#Format pioneer precursors 
pioneer_pr = DataFrame(Tables.columntable(Arrow.Table("/Users/n.t.wamsley/Desktop/ALTIMETER_PIONEER_101624/OlsenMixedSpeciesExploris500ng_M0M1M2/RESULTS/RESULTS/precursors_long.arrow")))
cols_to_keep = [
    :file_name, :species, :sequence, :structural_mods, 
    :charge, :peak_area
]
pioneer_pr = pioneer_pr[!,cols_to_keep]
pioneer_pr[!,:modified_sequence] = getModifiedSequence.(pioneer_pr[!,:sequence], "", pioneer_pr[!,:structural_mods])
pioneer_pr[!,:modified_sequence]  = pioneer_pr[!,:modified_sequence].*string.(pioneer_pr[!,:charge])
pioneer_pr = pioneer_pr[!,[:file_name,:species,:modified_sequence,:peak_area]]
rename!(pioneer_pr, [:file_name,:species,:modified_sequence,:peak_area] .=> [:Run, :Species, :PrecursorId,:PrecursorQuantity])
pioneer_pr[!,:Condition] = first.(split.(pioneer_pr[!,:Run],'_'))
pioneer_pr[!,:Experiment] = [condition_to_experiment[condition] for condition in pioneer_pr[!,:Condition]]
pioneer_pr[!,:QValue] .= zero(Float32)
pioneer_pr = pioneer_pr[!,[:Run,:Condition,:Experiment,:Species,:PrecursorId,:PrecursorQuantity,:QValue]]
filter!(x->!ismissing(x.PrecursorQuantity), pioneer_pr)
pioneer_pr[!,:PrecursorQuantity]= Float32.(pioneer_pr[!,:PrecursorQuantity])
pioneer_pr[!,:Condition] = string.(pioneer_pr[!,:Condition])
filter!(x->x.Run!="E5H50Y45_3",pioneer_pr)
filter!(x->!occursin("UNIMOD:35", x.PrecursorId),  pioneer_pr)


pr_condition_summary = combine(
    condition->summarizePrecursorCondition(condition),    
    groupby(pioneer_pr,[:Experiment,:Condition,:Species,:PrecursorId])
);
species_order = Dict("HUMAN" => 1, "YEAST" => 2, "ECOLI" => 3)
# Sort the DataFrame using the custom order
sort!(pr_condition_summary, :Species, by = x -> species_order[x])
precursor_group_stats = combine(groupby(pr_condition_summary,:Experiment)) do experiment_df
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
CSV.write(joinpath(out_dir, "stats", "precursors_stats.csv"), protein_group_stats)


##########
#Format pioneer protein groups 

pioneer_pg = DataFrame(Tables.columntable(Arrow.Table("/Users/n.t.wamsley/Desktop/ALTIMETER_PIONEER_101624/OlsenMixedSpeciesExploris500ng_M0M1M2/RESULTS/RESULTS/protein_groups_long.arrow")));
filter!(x->x.target, pioneer_pg)
filter!(x->x.n_peptides>=2, pioneer_pg)
cols_to_keep = [
    :file_name, :species, :protein, :abundance
]
pioneer_pg = pioneer_pg[!,cols_to_keep]
rename!(pioneer_pg, [:file_name,:species,:protein,:abundance] .=> [:Run, :Species, :ProteinGroup, :PGMaxLFQ])
pioneer_pg[!,:Condition] = first.(split.(pioneer_pg[!,:Run],'_'))
pioneer_pg[!,:Experiment] = [condition_to_experiment[condition] for condition in pioneer_pg[!,:Condition]]
pioneer_pg[!,:PGQValue] .= zero(Float32)
filter!(x->x.Run!="E5H50Y45_3",pioneer_pg)


pg_condition_summary = combine(
    condition->summarizeProteinGroupCondition(condition),    
    groupby(pioneer_pg,[:Experiment,:Condition,:Species,:ProteinGroup])
);
pg_condition_summary[!,:Condition] = string.(pg_condition_summary[!,:Condition])
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


#diann_results_path = "/Users/n.t.wamsley/Desktop/DIANN_ASTRAL/OlsenMixedSpeciesAstral200ngNonNormReport.parquet"
diann_results_path = "/Users/n.t.wamsley/Desktop/DIANN_EXPLORIS/OlsenMixedSpeciesExploris500ngReportNoNorm.parquet"
fasta_dir = "/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/SPEC_LIBS/FASTA/"
run_to_condition =  Dict(
    "20230324_OLEP08_200ng_30min_E10H50Y40_180K_2Th3p5ms_01.mzML" => "E10H50Y40",
    "20230324_OLEP08_200ng_30min_E10H50Y40_180K_2Th3p5ms_02.mzML" => "E10H50Y40",
    "20230324_OLEP08_200ng_30min_E10H50Y40_180K_2Th3p5ms_03.mzML" => "E10H50Y40",
    "20230324_OLEP08_200ng_30min_E20H50Y30_180K_2Th3p5ms_01.mzML" => "E20H50Y30",
    "20230324_OLEP08_200ng_30min_E20H50Y30_180K_2Th3p5ms_02.mzML" => "E20H50Y30",
    "20230324_OLEP08_200ng_30min_E20H50Y30_180K_2Th3p5ms_03.mzML" => "E20H50Y30",
    "20230324_OLEP08_200ng_30min_E30H50Y20_180K_2Th3p5ms_01.mzML" => "E30H50Y20",
    "20230324_OLEP08_200ng_30min_E30H50Y20_180K_2Th3p5ms_02.mzML" => "E30H50Y20",
    "20230324_OLEP08_200ng_30min_E30H50Y20_180K_2Th3p5ms_03.mzML" => "E30H50Y20",
    "20230324_OLEP08_200ng_30min_E40H50Y10_180K_2Th3p5ms_01.mzML" => "E40H50Y10",
    "20230324_OLEP08_200ng_30min_E40H50Y10_180K_2Th3p5ms_02.mzML" => "E40H50Y10",
    "20230324_OLEP08_200ng_30min_E40H50Y10_180K_2Th3p5ms_03.mzML" => "E40H50Y10",
    "20230324_OLEP08_200ng_30min_E45H50Y5_180K_2Th3p5ms_01.mzML" => "E45H50Y5",
    "20230324_OLEP08_200ng_30min_E45H50Y5_180K_2Th3p5ms_02.mzML" => "E45H50Y5",
    "20230324_OLEP08_200ng_30min_E45H50Y5_180K_2Th3p5ms_03.mzML" => "E45H50Y5",
    "20230324_OLEP08_200ng_30min_E5H50Y45_180K_2Th3p5ms_01.mzML" => "E5H50Y45",
    "20230324_OLEP08_200ng_30min_E5H50Y45_180K_2Th3p5ms_02.mzML" => "E5H50Y45",
    "20230324_OLEP08_200ng_30min_E5H50Y45_180K_2Th3p5ms_03.mzML" => "E5H50Y45"
)

run_to_condition =  Dict(
    "20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E10H50Y40_30SPD_DIA_1.raw" => "E10H50Y40",
    "20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E10H50Y40_30SPD_DIA_2.raw" => "E10H50Y40",
    "20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E10H50Y40_30SPD_DIA_3.raw" => "E10H50Y40",
    "20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E10H50Y40_30SPD_DIA_4.raw" => "E10H50Y40",
    "20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E20H50Y30_30SPD_DIA_1.raw" => "E20H50Y30",
    "20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E20H50Y30_30SPD_DIA_2.raw" => "E20H50Y30",
    "20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E20H50Y30_30SPD_DIA_3.raw" => "E20H50Y30",
    "20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E20H50Y30_30SPD_DIA_4.raw" => "E20H50Y30",
    "20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E30H50Y20_30SPD_DIA_1.raw" => "E30H50Y20",
    "20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E30H50Y20_30SPD_DIA_2.raw" => "E30H50Y20",
    "20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E30H50Y20_30SPD_DIA_3.raw" => "E30H50Y20",
    "20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E30H50Y20_30SPD_DIA_4.raw" => "E30H50Y20",
    "20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E40H50Y10_30SPD_DIA_1.raw" => "E40H50Y10",
    "20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E40H50Y10_30SPD_DIA_2.raw" => "E40H50Y10",
    "20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E40H50Y10_30SPD_DIA_3.raw" => "E40H50Y10",
    "20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E40H50Y10_30SPD_DIA_4.raw" => "E40H50Y10",
    "20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E45H50Y5_30SPD_DIA_1.raw" => "E45H50Y5",
    "20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E45H50Y5_30SPD_DIA_2.raw" => "E45H50Y5",
    "20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E45H50Y5_30SPD_DIA_3.raw" => "E45H50Y5",
    "20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E45H50Y5_30SPD_DIA_4.raw" => "E45H50Y5",
    "20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E5H50Y45_30SPD_DIA_1.raw" => "E5H50Y45",
    "20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E5H50Y45_30SPD_DIA_2.raw" => "E5H50Y45",
    "20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E5H50Y45_30SPD_DIA_3.raw" => "E5H50Y45",
    "20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E5H50Y45_30SPD_DIA_4.raw" => "E5H50Y45",
)


key_file_path = "/Users/n.t.wamsley/Desktop/astral_test_key_080324.txt"
condition_to_experiment, experiment_to_species_ratios = conditionToExperiment(key_file_path)

diann_r = LoadDiannResults(
    diann_results_path,
    fasta_dir,
    run_to_condition,
    condition_to_experiment
)
filter!(x->x["Run"]!="20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E5H50Y45_30SPD_DIA_3.raw",diann_r)


out_dir = "/Users/n.t.wamsley/Desktop/DIANN_OlsenExploris200ng"
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


pg_condition_summary = combine(
    condition->summarizeProteinGroupCondition(condition),    
    groupby(getProteinQuantDIANN(diann_r),[:Experiment,:Condition,:Species,:ProteinGroup])
);

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
#filter!(x->x.n>2, pg_condition_summary);

pr_condition_summary = combine(
    condition->summarizePrecursorCondition(condition),    
    groupby(getPrecursorQuantDIANN(diann_r),[:Experiment,:Condition,:Species,:PrecursorId])
);
#filter!(x->x.n>2, pr_condition_summary);

#for (key, experiment_df) in pairs(groupby(pr_condition_summary,:Experiment))
precursor_group_stats = combine(groupby(pr_condition_summary,:Experiment)) do experiment_df
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
CSV.write(joinpath(out_dir, "stats", "precursors_stats.csv"), protein_group_stats)

pr_experiments = groupby(getPrecursorQuantDIANN(testdf),:Experiment)


test_exp_summary = summarizeExperiment(
    "E5H50Y45_vs_E45H50Y5", 
    test_experiment,
    min_n = UInt8(3),
    [:Species,:PrecursorId],
    )

    test_exp_summary = summarizeExperiment(
        "E10H50Y40_vs_E40H50Y10", 
        groupby(pg_condition_summary,:Experiment)[1],
        [:Species,:ProteinGroup],
        experiment_to_species_ratios,
        min_n = UInt8(3)
        )


        test_exp_summary = summarizeExperiment(
            "E20H50Y30_vs_E30H50Y20", 
            groupby(pg_condition_summary,:Experiment)[2],
            [:Species,:ProteinGroup],
            experiment_to_species_ratios,
            min_n = UInt8(3)
            )

            