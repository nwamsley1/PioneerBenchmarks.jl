function mmccAnalysis()
    run_to_condition = Dict(
        "20210906_LRH_MMCC_Static_A_1.raw"=>"A",
        "20210906_LRH_MMCC_Static_A_2.raw"=>"A",
        "20210906_LRH_MMCC_Static_A_3.raw"=>"A",
        "20210906_LRH_MMCC_Static_B_1.raw"=>"B",
        "20210906_LRH_MMCC_Static_B_2.raw"=>"B",
        "20210906_LRH_MMCC_Static_B_3.raw"=>"B",
        "20210906_LRH_MMCC_Static_C_1.raw"=>"C",
        "20210906_LRH_MMCC_Static_C_2.raw"=>"C",
        "20210906_LRH_MMCC_Static_C_3.raw"=>"C",
        "20210906_LRH_MMCC_Static_D_1.raw"=>"D",
        "20210906_LRH_MMCC_Static_D_2.raw"=>"D",
        "20210906_LRH_MMCC_Static_D_3.raw"=>"D",
        "20210906_LRH_MMCC_Static_E_1.raw"=>"E",
        "20210906_LRH_MMCC_Static_E_2.raw"=>"E",
        "20210906_LRH_MMCC_Static_E_3.raw"=>"E",
        "20210906_LRH_MMCC_Static_F_1.raw"=>"F",
        "20210906_LRH_MMCC_Static_F_2.raw"=>"F",
        "20210906_LRH_MMCC_Static_F_3.raw"=>"F",
        "20210906_LRH_MMCC_Static_G_1.raw"=>"G",
        "20210906_LRH_MMCC_Static_G_2.raw"=>"G",
        "20210906_LRH_MMCC_Static_G_3.raw"=>"G",
        "20210906_LRH_MMCC_Static_H_1.raw"=>"H",
        "20210906_LRH_MMCC_Static_H_2.raw"=>"H",
        "20210906_LRH_MMCC_Static_H_3.raw"=>"H",
        "20210906_LRH_MMCC_Static_I_1.raw"=>"I",
        "20210906_LRH_MMCC_Static_I_2.raw"=>"I",
        "20210906_LRH_MMCC_Static_I_3.raw"=>"I",
        "20210906_LRH_MMCC_Static_J_1.raw"=>"J",
        "20210906_LRH_MMCC_Static_J_2.raw"=>"J",
        "20210906_LRH_MMCC_Static_J_3.raw"=>"J",
        "20210906_LRH_MMCC_Static_K_1.raw"=>"K",
        "20210906_LRH_MMCC_Static_K_2.raw"=>"K",
        "20210906_LRH_MMCC_Static_K_3.raw"=>"K"
    )

    conditions = ["A","B","C","D","E","F","G","H","I","J","K"]
    fractions = Float64[0, 0.5, 1, 3, 5, 7, 10, 30, 50, 70, 100]
    cond_to_fraction = Dict(
        zip(
            conditions,fractions
        )
    )


    diann_r = LoadDiannResultsMMCC(
        "/Users/n.t.wamsley/Desktop/DIANN_LUMOS_MMCC/LUMOS_MMCC_HEIL_7o30_report.parquet",
        run_to_condition
    )
    filter!(x->!occursin("-H",x.PrecursorId::String), diann_r)
    filter!(x->occursin("SILAC",x.PrecursorId::String), diann_r)
    filter!(x->x.QValue<=0.01,diann_r)
    filter!(x->!iszero(x.PrecursorQuantity), diann_r)
    precursors_table = getMMCCPrecursorQuantDIANN(diann_r)
    protein_groups_table = getMMCCProteinQuantDIANN(diann_r)
    precursors_table[!,:percent_light] = [Float64(cond_to_fraction[x]) for x in  precursors_table[!,:Condition]]
    protein_groups_table[!,:percent_light] = [Float64(cond_to_fraction[x]) for x in  protein_groups_table[!,:Condition]]

    println("Getting protein group curves...")
    protein_group_curves = analyzeMMCC(
        protein_groups_table,
        :ProteinIds,
        :PGMaxLFQ,
        power = 0.5,
        min_t = 0.5
    )
    println("Getting precursors curves...")
    precursor_curves = analyzeMMCC(
        precursors_table,
        :PrecursorId,
        :PrecursorQuantity,
        power = 0.5, 
        min_t = 0.5
    )
    return protein_group_curves, precursor_curves
end