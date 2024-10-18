function conditionToExperiment(
        experiments::Vector{SubString{String}}
)
    #= example
    julia> experiments
    3-element Vector{SubString{String}}:
    "E5H50Y45,E45H50Y5,ECOLI;9.0|HUMAN;1.0|YEAST;0.1111111"
    "E10H50Y40,E40H50Y10,ECOLI;4.0|HUMAN;1.0|YEAST;0.25"
    "E20H50Y30,E30H50Y20,ECOLI;1.5|HUMAN;1.0|YEAST;0.6666666"
    =#
    condition_to_experiment = Dict{String, String}()
    experiment_to_species_ratios = Dict{String, Dict{String, Float32}}()
    for experiment in experiments 
        if length(experiment)==0
            continue
        end
        condition_a, condition_b, species_to_ratio = split(experiment, ',')
        experiment_name = condition_a*"_vs_"*condition_b
        condition_to_experiment[condition_a] = experiment_name
        condition_to_experiment[condition_b] = experiment_name
        experiment_to_species_ratios[experiment_name] = Dict{String, Float32}()
        for spec_ratio in split(species_to_ratio,'|')
            spec, ratio = split(spec_ratio, ';')
            experiment_to_species_ratios[experiment_name][spec] = parse(Float32, ratio)
        end
    end

    #=
    julia> condition_to_experiment
    Dict{String, String} with 6 entries:
    "E45H50Y5"  => "E5H50Y45_vs_E45H50Y5"
    "E10H50Y40" => "E10H50Y40_vs_E40H50Y10"
    "E30H50Y20" => "E20H50Y30_vs_E30H50Y20"
    "E5H50Y45"  => "E5H50Y45_vs_E45H50Y5"
    "E40H50Y10" => "E10H50Y40_vs_E40H50Y10"
    "E20H50Y30" => "E20H50Y30_vs_E30H50Y20"

    julia> experiment_to_species_ratios
    Dict{String, Dict{String, Float32}} with 3 entries:
    "E20H50Y30_vs_E30H50Y20" => Dict("YEAST"=>0.666667, "ECOLI"=>1.5, "HUMAN"=>1.0)
    "E10H50Y40_vs_E40H50Y10" => Dict("YEAST"=>0.25, "ECOLI"=>4.0, "HUMAN"=>1.0)
    "E5H50Y45_vs_E45H50Y5"   => Dict("YEAST"=>0.111111, "ECOLI"=>9.0, "HUMAN"=>1.0)
    =#
    return condition_to_experiment, experiment_to_species_ratios
end

function conditionToExperiment(
    experiment_key_path::String
)

    content_ = read(experiment_key_path, String)
    #=
        julia> print(content_)
        #condition keys
        E10H50Y40,E40H50Y10,E5H50Y45,E45H50Y5,E20H50Y30,E30H50Y20
        #condition pairs
        E5H50Y45,E45H50Y5,ECOLI;9.0|HUMAN;1.0|YEAST;0.1111111
        E10H50Y40,E40H50Y10,ECOLI;4.0|HUMAN;1.0|YEAST;0.25
        E20H50Y30,E30H50Y20,ECOLI;1.5|HUMAN;1.0|YEAST;0.6666666
    =#
    return conditionToExperiment(split(content_, '\n')[4:end])
end

#key_file_path = "/Users/n.t.wamsley/Desktop/astral_test_key_080324.txt"