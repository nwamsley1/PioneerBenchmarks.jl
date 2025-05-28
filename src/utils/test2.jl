psms = load("/Users/n.t.wamsley/Desktop/huber_gpsms.jld2")["psms"]
gpsms = groupby(psms, [:precursor_idx,:scan_idx])
keys = []
for (key, sdf) in pairs(gpsms)
    if ((maximum(sdf[!,:weight]) - minimum(sdf[!,:weight]))/minimum(sdf[!,:weight])) > 0.1
        push!(keys, key)
    end

end

N = 1


test = process_huber_curve(gpsms[N][!,:weight], gpsms[N][!,:huber_δ])
plot(log2.(gpsms[N][!,:huber_δ]), gpsms[N][!,:weight], seriestype=:scatter)
vline!([log2(test[:huber50])])
hline!([test[:w50]])
N += 1


sgpsms = gpsms[keys[N]]
test = process_huber_curve(sgpsms[!,:weight], sgpsms[!,:huber_δ])
plot(log2.(sgpsms[!,:huber_δ]), sgpsms[!,:weight], seriestype=:scatter)
vline!([log2(test[:huber50])])
hline!([test[:w50]])
N += 1


function process_huber_curve(
    weights::AbstractVector{Float32},
    huber_δs::AbstractVector{Float32}
)
    min_w, max_w = minimum(weights), maximum(weights)
    huber50 = missing
    w50 = min_w + (max_w - min_w)/2
    
    if length(weights) > 1
        for i in 1:(length(weights)-1)
            if (w50 >= weights[i]) && (w50 <= weights[i + 1])
                println("i $i ", huber_δs[i])
                huber50 = huber_δs[i] + (huber_δs[i + 1] - huber_δs[i])/2
            elseif (w50 <= weights[i]) && (w50 >= weights[i + 1])
                println("i $i ", huber_δs[i])
                huber50 = huber_δs[i] + (huber_δs[i + 1] - huber_δs[i])/2
            end
        end
    end
    
    return (
        min = min_w,
        max = max_w,
        n = length(weights),
        huber50 = huber50,
        w50 = w50,
        wdiff = (max_w - min_w)/min_w
    )
end




pioneer_results_path = "/Volumes/d.goldfarb/Active/RIS_Goldfarb_Lab/NTW/PIONEER/PIONEER_PAPER/ANALYSES_PIONEER/SCIEX_NSWATH4/Jan01_2025_SEARCH"
fasta_dir = "/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/SPEC_LIBS/FASTA/"
key_file_path = "/Volumes/d.goldfarb/Active/RIS_Goldfarb_Lab/NTW/PIONEER/PIONEER_PAPER/ANALYSES_PIONEER/SCIEX_NSWATH4/Jan01_2025_SEARCH/key_AB.txt"
run_to_condition_path = "/Volumes/d.goldfarb/Active/RIS_Goldfarb_Lab/NTW/PIONEER/PIONEER_PAPER/ANALYSES_PIONEER/SCIEX_NSWATH4/Jan01_2025_SEARCH/run_to_condition.txt"
threeProteomeAnalysis(
    pioneer_results_path,
    fasta_dir,
    run_to_condition_path,
    key_file_path,
    "pioneer",
    "/Volumes/d.goldfarb/Active/RIS_Goldfarb_Lab/NTW/PIONEER/PIONEER_PAPER/ANALYSES_PIONEER/SCIEX_NSWATH4/Jan01_2025_SEARCH/three_proteome"
)

key_file_path = "/Volumes/d.goldfarb/Active/RIS_Goldfarb_Lab/NTW/PIONEER/PIONEER_PAPER/ANALYSES_PIONEER/SCIEX_NSWATH4/Jan01_2025_SEARCH/key_AC.txt"
threeProteomeAnalysis(
    pioneer_results_path,
    fasta_dir,
    run_to_condition_path,
    key_file_path,
    "pioneer",
    "/Volumes/d.goldfarb/Active/RIS_Goldfarb_Lab/NTW/PIONEER/PIONEER_PAPER/ANALYSES_PIONEER/SCIEX_NSWATH4/Jan01_2025_SEARCH/three_proteome"
)

key_file_path = "/Volumes/d.goldfarb/Active/RIS_Goldfarb_Lab/NTW/PIONEER/PIONEER_PAPER/ANALYSES_PIONEER/SCIEX_NSWATH4/Jan01_2025_SEARCH/key_BC.txt"
threeProteomeAnalysis(
    pioneer_results_path,
    fasta_dir,
    run_to_condition_path,
    key_file_path,
    "pioneer",
    "/Volumes/d.goldfarb/Active/RIS_Goldfarb_Lab/NTW/PIONEER/PIONEER_PAPER/ANALYSES_PIONEER/SCIEX_NSWATH4/Jan01_2025_SEARCH/three_proteome"
)
