function filterRowsOnSequenceLength!(
    protein_group_df::DataFrame,
    min_seq_length::Int64,
    max_seq_length::Int64)
    protein_group_df[!,:SequenceLength] = length.(protein_group_df[!,:StrippedSequence])
    filter!(x->(x.SequenceLength>min_seq_length)&(x.SequenceLength<max_seq_length), protein_group_df)
    return nothing
end

function filterForbiddenMods(
    modified_sequence::String,
    forbidden_mods::Vector{String})
    for mod in forbidden_mods
        if occursin(mod, modified_sequence)
            return true
        end
    end
    return false
end
function getPrecursorQuantDIANN(
    diann_df::DataFrame; 
    min_seq_length::Int64 = 7,
    max_seq_length::Int64 = 30,
    max_precursor_q_val::Real = 0.01,
    forbidden_mods::Vector{String} = ["UniMod:35"])

    #proteinGroupCols = Symbol.(["Condition","Biological","Technical","Species","Protein.Group", "Protein.Ids","Protein.Names","Genes","PG.Normalised","PG.MaxLFQ","PG.Q.Value"])

    precursorCols = [:Run,:Condition,:Experiment,:Species,:PrecursorId,:ModifiedSequence,:StrippedSequence,:PrecursorQuantity,:PrecursorCharge,:QValue]
    diann_pr = copy(diann_df[!,precursorCols])
    ###########
    #Filter Rows 
    filterRowsOnSequenceLength!(
        diann_pr,
        min_seq_length, max_seq_length,
    )
    #Remove length filtering Rows
    select!(diann_pr, Not([:SequenceLength]))
    #Remove rows not passing the q-value threshold or with no quantitation
    filter!(x->x.QValue<=max_precursor_q_val, diann_pr)
    filter!(x->!iszero(x.PrecursorQuantity), diann_pr)

    #Remove rows with forbidden variable modifications
    #(such as MOx which will not have reliable quantitation)
    diann_pr[!,:has_forbidden_mods] = map(x->filterForbiddenMods(x, forbidden_mods), diann_pr[!,:ModifiedSequence])
    filter!(x->x.has_forbidden_mods==false, diann_pr)
    select!(diann_pr, Not([:has_forbidden_mods]))

    return diann_pr
end

function getMMCCPrecursorQuantDIANN(
    diann_df::DataFrame; 
    min_seq_length::Int64 = 7,
    max_seq_length::Int64 = 30,
    max_precursor_q_val::Real = 0.01,
    forbidden_mods::Vector{String} = ["UniMod:35"])

    #proteinGroupCols = Symbol.(["Condition","Biological","Technical","Species","Protein.Group", "Protein.Ids","Protein.Names","Genes","PG.Normalised","PG.MaxLFQ","PG.Q.Value"])

    precursorCols = [:Run,:Condition,:PrecursorId,:ModifiedSequence,:StrippedSequence,:PrecursorQuantity,:PrecursorCharge,:QValue]
    diann_pr = copy(diann_df[!,precursorCols])
    ###########
    #Filter Rows 
    filterRowsOnSequenceLength!(
        diann_pr,
        min_seq_length, max_seq_length,
    )
    #Remove length filtering Rows
    select!(diann_pr, Not([:SequenceLength]))
    #Remove rows not passing the q-value threshold or with no quantitation
    filter!(x->x.QValue<=max_precursor_q_val, diann_pr)
    filter!(x->!iszero(x.PrecursorQuantity), diann_pr)

    #Remove rows with forbidden variable modifications
    #(such as MOx which will not have reliable quantitation)
    diann_pr[!,:has_forbidden_mods] = map(x->filterForbiddenMods(x, forbidden_mods), diann_pr[!,:ModifiedSequence])
    filter!(x->x.has_forbidden_mods==false, diann_pr)
    select!(diann_pr, Not([:has_forbidden_mods]))
    diann_pr[!,:PrecursorQuantity] = Float64.(diann_pr[!,:PrecursorQuantity])
    return diann_pr
end
