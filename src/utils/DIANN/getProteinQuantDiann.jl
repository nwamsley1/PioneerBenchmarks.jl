function filterRowsOnSequenceLength!(
    protein_group_df::DataFrame,
    min_seq_length::Int64,
    max_seq_length::Int64)
    protein_group_df[!,:SequenceLength] = length.(protein_group_df[!,:StrippedSequence])
    filter!(x->(x.SequenceLength>=min_seq_length)&(x.SequenceLength<=max_seq_length), protein_group_df)
    return nothing
end


function getPrecursorsPerProtein!(
    df::DataFrame, 
    group_cols::Vector{Symbol})
    # Group by the two columns
    sort!(df, group_cols)
    df[!,:n_precursors] .= 0
    df[!,:to_keep] .= false
    gdf = groupby(df, group_cols)
    for sdf in gdf
        sdf[1,:n_precursors] = length(unique(sdf[!,:StrippedSequence]))#size(sdf, 1)
        sdf[1,:to_keep] = true
    end
    # Get the first row of each group and count the occurrences
    filter!(x->x.to_keep, df)
    select!(df, Not(:to_keep))
    return nothing
end

function getProteinQuantDIANN(
    diann_df::DataFrame; 
    min_seq_length::Int64 = 7,
    max_seq_length::Int64 = 30,
    max_precursor_q_val::Real = 0.01,
    max_pg_q_val::Real = 0.01,
    min_precursors_per_protein::Integer = 2)

    #proteinGroupCols = Symbol.(["Condition","Biological","Technical","Species","Protein.Group", "Protein.Ids","Protein.Names","Genes","PG.Normalised","PG.MaxLFQ","PG.Q.Value"])

    proteinGroupCols = [:Run,:Condition,:Experiment,:Species,:StrippedSequence,:PrecursorQuantity,:QValue,:ProteinGroup,:ProteinIds,:ProteinNames,:PGNormalised,:PGMaxLFQ,:PGQValue]
    diann_pg = copy(diann_df[!,proteinGroupCols])
    ###########
    #Filter Rows 
    filterRowsOnSequenceLength!(
        diann_pg,
        min_seq_length, max_seq_length,
    )
    #Remove length filtering Rows
    #select!(diann_pg, Not([:StrippedSequence,:SequenceLength]))
    select!(diann_pg, Not([:SequenceLength]))
    #Remove rows not passing the q-value threshold or with no quantitation
    filter!(x->x.QValue<=max_precursor_q_val, diann_pg)
    filter!(x->!iszero(x.PrecursorQuantity), diann_pg)
    select!(diann_pg, Not([:PrecursorQuantity,:QValue]))
    #How many precursors identified for each protein group?
    #Retain only one row per protein group per run 
    getPrecursorsPerProtein!(diann_pg,
    [:Run,:ProteinGroup])

    #Remove protein groups identified based on fewer than a given number of precursors
    filter!(x->x.n_precursors>=min_precursors_per_protein, diann_pg)
    #Remove precursor count column 
    select!(diann_pg, Not([:n_precursors]))

    #Filter based on protein group q value
    filter!(x->x.PGQValue<=max_pg_q_val, diann_pg)
    #println("TEST")
    #Filter based on presence/absence of protein group quantitation
    filter!(x->!iszero(x.PGMaxLFQ), diann_pg)

    return diann_pg
end

function getMMCCProteinQuantDIANN(
    diann_df::DataFrame; 
    min_seq_length::Int64 = 7,
    max_seq_length::Int64 = 30,
    max_precursor_q_val::Real = 0.01,
    max_pg_q_val::Real = 0.01,
    min_precursors_per_protein::Integer = 2)

    #proteinGroupCols = Symbol.(["Condition","Biological","Technical","Species","Protein.Group", "Protein.Ids","Protein.Names","Genes","PG.Normalised","PG.MaxLFQ","PG.Q.Value"])

    proteinGroupCols = [:Run,:Condition,:StrippedSequence,:PrecursorQuantity,:QValue,:ProteinGroup,:ProteinIds,:ProteinNames,:PGNormalised,:PGMaxLFQ,:PGQValue]
    diann_pg = copy(diann_df[!,proteinGroupCols])
    ###########
    #Filter Rows 
    filterRowsOnSequenceLength!(
        diann_pg,
        min_seq_length, max_seq_length,
    )
    #Remove length filtering Rows
    #select!(diann_pg, Not([:StrippedSequence,:SequenceLength]))
    select!(diann_pg, Not([:SequenceLength]))
    #Remove rows not passing the q-value threshold or with no quantitation
    filter!(x->x.QValue<=max_precursor_q_val, diann_pg)
    filter!(x->!iszero(x.PrecursorQuantity), diann_pg)
    select!(diann_pg, Not([:PrecursorQuantity,:QValue]))
    #How many precursors identified for each protein group?
    #Retain only one row per protein group per run 
    getPrecursorsPerProtein!(diann_pg,
    [:Run,:ProteinGroup])

    #Remove protein groups identified based on fewer than a given number of precursors
    filter!(x->x.n_precursors>=min_precursors_per_protein, diann_pg)
    #Remove precursor count column 
    select!(diann_pg, Not([:n_precursors]))
    #Filter based on protein group q value
    #println("TEST")
    filter!(x->x.PGQValue<=max_pg_q_val, diann_pg)

    #Filter based on presence/absence of protein group quantitation
    filter!(x->!iszero(x.PGMaxLFQ), diann_pg)
    diann_pg[!,:PGMaxLFQ] = Float64.(diann_pg[!,:PGMaxLFQ])
    return diann_pg
end
