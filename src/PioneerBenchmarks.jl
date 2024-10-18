module PioneerBenchmarks
    __precompile__(false)
    using Base.Filesystem
    using CSV, DataFrames
    using Parquet2: Dataset
    using Dictionaries
    using HypothesisTests, StatsBase, StatsPlots, Measures, Plots
    using FASTX, CodecZlib
    using MultipleTesting
    # Write your package code here.
    include(joinpath(@__DIR__,"utils","importScripts.jl"))
    importScripts()
    export threeProteomeAnalysis
end
