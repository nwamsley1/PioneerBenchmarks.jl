function summarizeProteinGroupCondition(
        condition_df::AbstractDataFrame
    )
    summarizeProteinGroupCondition(
        condition_df[!,:PGMaxLFQ],
        condition_df[!,:PGQValue],
    )
end

function summarizeProteinGroupCondition(
        PGMaxLFQ::AbstractVector{Float32},
        PGQValue::AbstractVector{Float32},
    )::@NamedTuple{mean::Float32, CV::Float32, std::Float32, var::Float32, n::UInt8, max_qval::Float32}
    _mean_, _std_, _n_ =  mean(PGMaxLFQ),std(PGMaxLFQ),length(PGMaxLFQ)
    return (
        mean = _mean_,
        CV = getCV(_mean_, _std_, _n_),
        std = _std_,
        var = var(PGMaxLFQ),
        n = _n_,
        max_qval = maximum(PGQValue)
    )
end

function summarizePrecursorCondition(
    condition_df::SubDataFrame
)
    summarizePrecursorCondition(
    condition_df[!,:PrecursorQuantity],
    condition_df[!,:QValue]
)
end

function summarizePrecursorCondition(
    PrecursorQuantity::AbstractVector{Float32},
    QValue::AbstractVector{Float32}
)::@NamedTuple{mean::Float32, CV::Float32, std::Float32, var::Float32, n::UInt8, max_qval::Float32}
    _mean_, _std_, _n_ =  mean(PrecursorQuantity),std(PrecursorQuantity),length(PrecursorQuantity)
    return (
        mean = _mean_,
        CV = getCV(_mean_, _std_, _n_),
        std = _std_,
        var = var(PrecursorQuantity),
        n = _n_,
        max_qval = maximum(QValue)
    )
end
