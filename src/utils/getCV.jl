function getCV(
    _mean_::T,
    _std_::T,
    _n_::I
) where {T<:AbstractFloat,I<:Integer}
    return Float32(100*((1 + 1/(4*_n_))*_std_/_mean_))
end