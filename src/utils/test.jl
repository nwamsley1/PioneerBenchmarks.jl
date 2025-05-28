diann_results_path = "/Users/n.t.wamsley/Desktop/DIANN_EXPLORIS/OlsenMixedSpeciesExploris500ngReportNoNorm.parquet"
fasta_dir = "/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/SPEC_LIBS/FASTA/"
key_file_path = "/Users/n.t.wamsley/Desktop/astral_test_key_080324.txt"
run_to_condition_path = "/Users/n.t.wamsley/Projects/PioneerBenchmarks/data/example_run_to_condition.tsv"
threeProteomeAnalysis(
    diann_results_path,
    fasta_dir,
    run_to_condition_path,
    key_file_path,
    "diann",
    "/Users/n.t.wamsley/Desktop/PioneerBenchmarksTest"
)

pioneer_results_path = "/Users/n.t.wamsley/Desktop/AstralAltimeter_101624_M0M1/"
fasta_dir = "/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/SPEC_LIBS/FASTA/"
run_to_condition_path = "/Users/n.t.wamsley/Projects/PioneerBenchmarks/data/example_run_to_condition_pioneer.csv"
key_file_path = "/Users/n.t.wamsley/Desktop/astral_test_key_080324.txt"
threeProteomeAnalysis(
    pioneer_results_path,
    fasta_dir,
    run_to_condition_path,
    key_file_path,
    "pioneer",
    "/Users/n.t.wamsley/Desktop/PioneerBenchmarksTest2"
)



conditions = ["A","B","C","D","E","F","G","H","I","J","K"]
fractions = Float64[0, 0.5, 1, 3, 5, 7, 10, 30, 50, 70, 100]
cond_to_fraction = Dict(
    zip(
        conditions,fractions
    )
)

function create_12_ranges(N::Int)
    # Calculate size of each chunk (rounded up to include all elements)
    chunk_size = ceil(Int, N/12)
    
    # Create the ranges
    ranges = [((i-1)*chunk_size + 1):min(i*chunk_size, N) for i in 1:12]
    
    return ranges
end

using StatsBase, Distributions, LinearAlgebra, StaticArrays, DataFrames, Polynomials, ProgressBars
protein_groups_table, precursors_table = mmccAnalysis()


precursors_table[!,:PrecursorQuantity] = Float64.( precursors_table[!,:PrecursorQuantity])
precursors_table[!,:percent_light] = [Float64(cond_to_fraction[x]) for x in  precursors_table[!,:Condition]]
sort!( precursors_table, [:PrecursorId,:Condition])
precursors_table_grouped = groupby( precursors_table, :PrecursorId)
include("src/utils/MMCC_UTILS/uniformBspline.jl")
include("src/utils/MMCC_UTILS/bootstrapMMCC.jl")


ranges = create_12_ranges(length(precursors_table_grouped))
test_vecs = Vector{Any}(undef, 12)
@time begin
Threads.@threads for (idx, gdf_range) in collect(enumerate(ranges))
    test_vecs[idx] = BootstrapGDF(precursors_table_grouped[gdf_range], 
    :PrecursorQuantity, 
    power = 0.5,
    n_bootstraps = 1000,
    n_xbins = 100)
end
end
curves = vcat(test_vecs...)
filter!(x->coalesce(x.lloq, typemax(Float32))<99.0, curves)


protein_groups_table[!,:PGMaxLFQ] = Float64.( protein_groups_table[!,:PGMaxLFQ])
protein_groups_table[!,:percent_light] = [Float64(cond_to_fraction[x]) for x in  protein_groups_table[!,:Condition]]
sort!(protein_groups_table, [:ProteinIds,:Condition])
protein_groups_table_grouped = groupby(protein_groups_table, :ProteinIds)

lloq_curves_pg = analyzeMMCC(
    protein_groups_table,
    :ProteinIds,
    :PGMaxLFQ
)


ranges = create_12_ranges(length(protein_groups_table_grouped))
test_vecs = Vector{Any}(undef, 12)
@time begin
Threads.@threads for (idx, gdf_range) in collect(enumerate(ranges))
    test_vecs[idx] = BootstrapGDF(protein_groups_table_grouped[gdf_range], 
    :ProteinIds,
    :PGMaxLFQ, 
    power = 0.5,
    n_bootstraps = 1000,
    n_xbins = 100)
end
end
curves = vcat(test_vecs...)
filter!(x->coalesce(x.lloq, typemax(Float32))<99.0, curves)

histogram(collect(skipmissing(curves[!,:lloq])))
histogram(collect(skipmissing(curves[!,:aic_mean])))
histogram(collect(skipmissing(curves[!,:mean_intercept]))./(collect(skipmissing(curves[!,:max_at_100]))), bins = LinRange(-0.3, 0.3, 25))


protein_groups_table = getProteinQuantDIANN(diann_r)
Threads.@threads for (idx, gdf_range) in collect(enumerate(ranges))
    println("urmom")
end

include("src/utils/MMCC_UTILS/bootstrapMMCC.jl")
test_ = BootstrapGDF(testdf2_grouped[1:100], 
:PrecursorQuantity, 
power = 0.5,
n_bootstraps = 1000,
n_xbins = 100)

cond_to_fraction = Dict(
    zip(
        ["A","B","C","D","E","F","G","H","I","J","K"],
        Float32[0, 0, 0.25, 0.5, 1, 5, 10, 25, 50, 70, 100 ]
    )
)
test_example = copy(testdf2_grouped[100])
test_example[!,:percent_light] = [cond_to_fraction[x] for x in test_example[!,:Condition]]

include("src/utils/MMCC_UTILS/uniformBspline.jl")
X = UniformSplineDMAT(
    test_example[!,:percent_light],
    3,
    3
)
b = X[:,1:5]\test_example[!,:PrecursorQuantity]

y_fit = X[:,1:5]*b

plot( test_example[!,:percent_light], test_example[!,:PrecursorQuantity], seriestype=:scatter)

tspline = UniformSpline(
    test_example[!,:PrecursorQuantity],
    test_example[!,:percent_light],
    3,
    3
)

plot!( LinRange(0, 100, 1000), tspline.(LinRange(0, 100, 1000)))

gdf = testdf2_grouped
quant_col = :PrecursorQuantity
n_bootstraps = 1000
n_xbins = 100
xmin = 0
xmax = 100
power = 0.59
min_t = 0.01


T = Float64
curves = [(precursor_idx = UInt32(0),
           lloq = 0.0,
           aic_mean = 0.0,
           mean_intercept = 0.0,
           intercept_q95 = 0.0,
           max_at_100 = 0.0) for _ in range(1, length(gdf))]

bootstrap_Y = zeros(T, n_xbins, n_bootstraps);
bootstrap_intercept = zeros(T, n_bootstraps);
aic_diffs = zeros(T, n_bootstraps);
x_bins = collect(LinRange(xmin, xmax, n_xbins));
w_bins = diag(getWeights(x_bins, power, min_t));
x_mat_bins_5p = UniformSplineDMAT(x_bins, 3, 3)[:,1:5]
x_mat_bins_lin = hcat(
                    ones(T, length(x_bins)),
                    x_bins
        )

test_example = copy(gdf[100])
test_example[!,:percent_light] = [cond_to_fraction[x] for x in test_example[!,:Condition]]

data_x = Float64.(test_example[!,:percent_light])
#data_y = T.(coalesce.(psms[!,:gmean_area], 0.0))
data_y = Float64.(test_example[!,quant_col])


include("src/utils/MMCC_UTILS/bootstrapMMCC.jl")

i = 1
for (key, psms) in pairs(gdf)
    psms = gdf[12]
    data_x = T.(collect(psms[!,:percent_light]))
    #data_y = T.(coalesce.(psms[!,:gmean_area], 0.0))
    data_y = T.(coalesce.(psms[!,quant_col], 0.0))
BoostrapCalibration!(
            #zero(UInt32),
            zero(UInt32),
            bootstrap_Y,
            bootstrap_intercept,
            aic_diffs,
            data_x,
            data_y,
            x_mat_bins_5p,
            x_mat_bins_lin,
            w_bins,
            x_bins,
            weight_inv_power = power,
            min_t = min_t,
        )





diann_results_path = "/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/DATASETS_ARROW/ASTRAL_MTAC/DIANN/THREE_PROTEOME/report.parquet"
fasta_dir = "/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/DATASETS_ARROW/ASTRAL_MTAC/DIANN/THREE_PROTEOME/fasta"
key_file_path = "/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/DATASETS_ARROW/ASTRAL_MTAC/DIANN/THREE_PROTEOME/key.txt"
run_to_condition_path = "/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/DATASETS_ARROW/ASTRAL_MTAC/DIANN/THREE_PROTEOME/run_to_condition.txt"
threeProteomeAnalysis(
    diann_results_path,
    fasta_dir,
    run_to_condition_path,
    key_file_path,
    "diann",
    "/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/DATASETS_ARROW/ASTRAL_MTAC/DIANN/THREE_PROTEOME"
)






pioneer_results_path = "/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/DATASETS_ARROW/ASTRAL_MTAC/THREE_PROTEOME_5MIN/SEPERATE_TRACES"
fasta_dir = "/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/SPEC_LIBS/FASTA/"
key_file_path = "/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/DATASETS_ARROW/ASTRAL_MTAC/THREE_PROTEOME_5MIN/SEPERATE_TRACES/key.txt"
run_to_condition_path = "/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/DATASETS_ARROW/ASTRAL_MTAC/THREE_PROTEOME_5MIN/SEPERATE_TRACES/run_to_condition.txt"
threeProteomeAnalysis(
    pioneer_results_path,
    fasta_dir,
    run_to_condition_path,
    key_file_path,
    "pioneer",
    "/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/DATASETS_ARROW/ASTRAL_MTAC/THREE_PROTEOME_5MIN/SEPERATE_TRACES/three_proteome"
)



pioneer_results_path = "/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/DATASETS_ARROW/ASTRAL_MTAC/THREE_PROTEOME_5MIN/COMBINE_TRACES/"
fasta_dir = "/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/SPEC_LIBS/FASTA/"
key_file_path = "/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/DATASETS_ARROW/ASTRAL_MTAC/THREE_PROTEOME_5MIN/COMBINE_TRACES/key.txt"
run_to_condition_path = "/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/DATASETS_ARROW/ASTRAL_MTAC/THREE_PROTEOME_5MIN/COMBINE_TRACES/run_to_condition.txt"
threeProteomeAnalysis(
    pioneer_results_path,
    fasta_dir,
    run_to_condition_path,
    key_file_path,
    "pioneer",
    "/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/DATASETS_ARROW/ASTRAL_MTAC/THREE_PROTEOME_5MIN/COMBINE_TRACES/three_proteome"
)



pioneer_results_path = "/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/DATASETS_ARROW/OlsenMixedSpeciesAstral200ng/SEPERATE_TRACES/"
fasta_dir = "/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/SPEC_LIBS/FASTA/"
key_file_path = "/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/DATASETS_ARROW/OlsenMixedSpeciesAstral200ng/SEPERATE_TRACES/key.txt"
run_to_condition_path = "/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/DATASETS_ARROW/OlsenMixedSpeciesAstral200ng/SEPERATE_TRACES/run_to_condition.txt"
threeProteomeAnalysis(
    pioneer_results_path,
    fasta_dir,
    run_to_condition_path,
    key_file_path,
    "pioneer",
    "/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/DATASETS_ARROW/OlsenMixedSpeciesAstral200ng/SEPERATE_TRACES/three_proteome"
)



pioneer_results_path = "/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/PIONEER_ANALYSES/Jan01_2025_SEARCH_NoQuadFit"
fasta_dir = "/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/SPEC_LIBS/FASTA/"
key_file_path = "/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/PIONEER_ANALYSES/Jan01_2025_SEARCH_NoQuadFit/key.txt"
run_to_condition_path = "/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/PIONEER_ANALYSES/Jan01_2025_SEARCH_NoQuadFit/run_to_condition.txt"
threeProteomeAnalysis(
    pioneer_results_path,
    fasta_dir,
    run_to_condition_path,
    key_file_path,
    "pioneer",
    "/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/PIONEER_ANALYSES/Jan01_2025_SEARCH_NoQuadFit/three_proteome"
)
