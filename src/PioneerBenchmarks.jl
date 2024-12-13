module PioneerBenchmarks
    __precompile__(false)
    using Arrow
    using Base.Filesystem
    using CSV, DataFrames
    using StatsBase, Distributions
    using StaticArrays
    using Polynomials
    using LinearAlgebra
    using Parquet2: Dataset
    using Dictionaries
    using HypothesisTests, StatsBase, StatsPlots, Measures, Plots
    using FASTX, CodecZlib
    using MultipleTesting: BenjaminiHochberg
    using MultipleTesting: adjust
    # Write your package code here.
    include(joinpath(@__DIR__,"utils","importScripts.jl"))
    include(joinpath(@__DIR__,"threeProteomeAnalysis.jl"))
    include(joinpath(@__DIR__,"mmcc_analysis.jl"))
    importScripts()
    export threeProteomeAnalysis, mmccAnalysis
end
