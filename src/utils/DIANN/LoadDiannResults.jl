function getSpeciesFromProteinNames(
    protein_ids::String,
    acc_to_spec::Dict{String, String}
)::Vector{String}
    species_names = Set{String}()
    for protein_id in split(protein_ids, ';')
        push!(species_names, 
            acc_to_spec[protein_id]
        )
    end
    return collect(species_names)
end

function checkAndFilterDiannColumns(
    diann_df::DataFrame)
    
    required_cols = [
        "Run",
        "Precursor.Id","Modified.Sequence","Stripped.Sequence","Precursor.Charge",
        "Precursor.Quantity","Precursor.Normalised","Q.Value",
        "Protein.Ids","Protein.Group","Protein.Names",
        "PG.Normalised","PG.MaxLFQ","PG.Q.Value"]

    cols_in_df = Set(names(diann_df))
    missing_columns = setdiff(Set(required_cols), cols_in_df)
    if length(missing_columns)>0
        throw("Missing colmns in diann results "*string(collect(missing_columns)))
    end
    diann_df = diann_df[!,collect(required_cols)]
    rename!(diann_df, [Symbol(replace(x, "."=>"")) for x in names(diann_df)])
    return diann_df
end

function getCondition(
    run::String,
    run_to_condition::Dict{String, String})
    return run_to_condition[run]
end

function getExperiment(
    condition::String,
    condition_to_experiment::Dict{String, String})
    if haskey(condition_to_experiment, condition)
        return condition_to_experiment[condition]
    else
        return ""
    end
end


function LoadDiannResults(
    diann_results_path::String,
    fasta_dir::String,
    run_to_condition::Dict{String, String},
    condition_to_experiment::Dict{String, String})

    _, extension = splitext(diann_results_path)
    if extension != ".parquet"
        throw("DIANN results path should have .parquet extention. Path supplied was: $diann_results_path")
    end
    println("Loading DIANN Results DF...")
    @time begin
        ds = Dataset(diann_results_path)
        diann_r = DataFrame(ds; copycols=true);
    end
    println("Checking and filtering columns...")
    diann_r = checkAndFilterDiannColumns(diann_r)
    println("Getting Accesion ID to Species Map From Fastas...")
    @time begin
        acc_to_spec = getAccessionToSpecies(fasta_dir)
    end
    #Get unique species names for each row
    diann_r[!,:SpeciesNames] = map(x->getSpeciesFromProteinNames(x, acc_to_spec), diann_r[!,:ProteinIds])
    #Remove rows for precursors that are not unique to a single species
    filter!(x->length(x.SpeciesNames)==1, diann_r)
    #Get the species to which each entry belongs 
    diann_r[!,:Species] = map(x->first(x), diann_r[!,:SpeciesNames])
    diann_r[!,:Condition] = map(x->getCondition(x, run_to_condition), diann_r[!,:Run])
    #diann_r[!,:Experiment] = map(x->getExperiment(x, condition_to_experiment), diann_r[!,:Condition])


    diann_r[!,:Experiment] .= ""
    for i in range(1, size(diann_r, 1))
            diann_r[i,:Experiment] = getExperiment(diann_r[i,:Condition], condition_to_experiment)
    end
    filter!(x->x.Experiment.!="", diann_r)

    columns_order = [
    :Run,:Condition,:Experiment,:Species,
    :PrecursorId,:ModifiedSequence,:StrippedSequence,:PrecursorCharge,:PrecursorQuantity,:PrecursorNormalised,:QValue,
    :ProteinIds,:ProteinGroup,:ProteinNames,:PGNormalised,:PGMaxLFQ,:PGQValue]

    return diann_r[!,columns_order]
end

function LoadDiannResultsMMCC(
    diann_results_path::String,
    run_to_condition::Dict{String, String})
    
    _, extension = splitext(diann_results_path)
    if extension != ".parquet"
        throw("DIANN results path should have .parquet extention. Path supplied was: $diann_results_path")
    end
    println("Loading DIANN Results DF...")
    @time begin
        ds = Dataset(diann_results_path)
        diann_r = DataFrame(ds; copycols=true);
    end
    println("Checking and filtering columns...")
    diann_r = checkAndFilterDiannColumns(diann_r)
    println("Getting Accesion ID to Species Map From Fastas...")
    diann_r[!,:Condition] = map(x->getCondition(x, run_to_condition), diann_r[!,:Run])

    columns_order = [
    :Run,:Condition,
    :PrecursorId,:ModifiedSequence,:StrippedSequence,:PrecursorCharge,:PrecursorQuantity,:PrecursorNormalised,:QValue,
    :ProteinIds,:ProteinGroup,:ProteinNames,:PGNormalised,:PGMaxLFQ,:PGQValue]

    return diann_r[!,columns_order]
end