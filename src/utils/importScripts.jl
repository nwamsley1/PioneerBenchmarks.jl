function importScripts()
    package_root = dirname(dirname(@__DIR__))
    
    [include(joinpath(package_root, "src", "utils", "DIANN", jl_file)) for jl_file in [
                                "LoadDiannResults.jl",
                                "getPrecursorQuantDiann.jl",
                                "getProteinQuantDiann.jl",
                                "summarizeCondition.jl"]];
    
    [include(joinpath(package_root, "src", "utils", "PIONEER", jl_file)) for jl_file in [
        "getPrecursorQuantPioneer.jl",
        "getProteinQuantPioneer.jl",
        "insertMods.jl"]];
    #Utilities
    [include(joinpath(package_root, "src","utils", jl_file)) for jl_file in [
        "checkHasColumns.jl",
        "compareConditions.jl",
        "getCV.jl",
        "importScripts.jl",
        "parseExperimentalDesign.jl",
        "parseFasta.jl",
        "plotThreeProteome.jl",
        "summarizeExperiment.jl",
        "threeProteomeAnalysis.jl"]];
    println("TEST")
end