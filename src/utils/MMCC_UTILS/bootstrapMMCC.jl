function getWeights(t::AbstractVector{T}, power::U, min_t::T) where {T,U<:AbstractFloat}
    #weight w_i should be w_i = 1/σ_i^2. 
    #emperically relationship between t and 
    weights = zeros(T, (length(t), length(t)))
    for i in range(1, length(t))
        adjusted_t = max(t[i], min_t)
        weights[i, i] = 1/(adjusted_t^(power))
    end
    return weights
end

function updateβ!(β::Vector{T}, M::Matrix{T}, y::Vector{T}) where {T<:AbstractFloat}
    for (col_idx, col) in enumerate(eachcol(M))
        β[col_idx] = dot(col, y)
    end
    return nothing
end

function resampleData!(resampled_x::Vector{T},resampled_y::Vector{T}, data_x::AbstractVector{T}, data_y::AbstractVector{T}) where {T<:AbstractFloat}
    #Assums data_x and data_y are sorted in ascending order in by data_x  
    i = 1
    while i < length(data_x)
        xi = data_x[i]
        k = i
        while (data_x[k] - xi) < 1e-5
            k += 1
            if k > length(data_x)
                break
            end
        end
        k = k - 1
        for n in range(i,k)
            rand_idx = rand(i:k)
            resampled_y[n] = data_y[rand_idx]
            resampled_x[n] = data_x[rand_idx]
        end
        i = k + 1
        if i > length(data_x)
            return nothing 
        end
    end
    return nothing
end


function fillWeightedRedisuals!(weighted_r::Vector{T}, 
                                weights::Vector{T},
                                ŷ::Vector{T},
                                y::Vector{T}) where {T<:AbstractFloat}
    for i in range(1, length(weighted_r))
        weighted_r[i] = weights[i]*(ŷ[i] - y[i])
    end
    return nothing
end

function getCorrectedSigma(weighted_r::Vector{T},x::AbstractVector{T},min_x::T) where {T<:AbstractFloat}
    #sum_residuals = zero(T)
    #=
    sum_squared_residuals = zero(T)
    n = 0
    for i in range(1,  length(x))
        if x[i] >= min_x
            #sum_residuals += weighted_r[i]
            sum_squared_residuals += weighted_r[i]^2# - mean
        end
        n += 1
    end
    =#
    #Estimate of standard deviation of residuals
    #mean is assumed zero
    return StatsBase.mad(weighted_r, normalize = true)#sqrt(sum_squared_residuals/(n - 1))
end


function predictXBins(Ŷ::Matrix{T}, Ŷ_col::Int, X::Matrix{T}, β::Vector{T}) where {T<:AbstractFloat}
    for col in range(1, size(X, 2))
        for row in range(1, size(Ŷ, 1))
            Ŷ[Ŷ_col] += β[col]*X[row, col]
        end
    end
end

function solveβevalŷ!(β::Vector{T},
                    ŷ::Vector{T},
                    X::Matrix{T},
                    luA::LU{T, Matrix{T}, Vector{Int64}},
                    F::Vector{T},
                    XtW::Matrix{T},
                    y::Vector{T}) where {T<:AbstractFloat}
    mul!(F, XtW, y)
    ldiv!(β, luA, F)
    mul!(ŷ, X, β);
end

function predictXBins!(Y::Matrix{T}, bootstrap_idx::Int, X::Matrix{T}, β::Vector{T}) where {T<:AbstractFloat}
    for col in range(1, size(X, 2))
        for row in range(1, size(Y, 1))
            Y[row, bootstrap_idx] += β[col]*X[row, col]
        end
    end
end

function predictIntercept!(intercepts::Vector{T}, bootstrap_idx::Int, X::Matrix{T}, β::Vector{T}) where {T<:AbstractFloat}
    for col in range(1, length(β))
        intercepts[bootstrap_idx] += β[col]*X[1, col]
    end
end



function sampleData!(Y::Matrix{T}, bootstrap_idx::Int64, weights::Vector{T}, σ::T) where {T<:AbstractFloat}
            #Simulate data 
        for (row, w) in enumerate(weights)
             Y[row, bootstrap_idx] += max(rand(
                                Normal(0, σ/w)
                            ), zero(T))
        end
end

function getAIC(residuals::Vector{T}, σ::T, k::Int64) where {T<:AbstractFloat}
    aic = 2*k
    err_model = Distributions.Normal(0, σ)
    for r in residuals
        aic -= 2*logpdf(err_model, r)
    end
    n = length(residuals)
    return aic + (2*k^2 + 2*k)/(n - k - 1)
end

function BoostrapCalibration!(precursor_idx::String,
                             bootstrap_Y::Matrix{T},
                             bootstrap_intercept::Vector{T},
                             aic_diffs::Vector{T},
                             x_data::AbstractArray{U}, 
                             y_data::AbstractArray{U},
                             x_mat_bins_5p::Matrix{T},
                             x_mat_bins_lin::Matrix{T},
                             w_bins::Vector{T},
                             x_bins::Vector{T};
                             max_cv::Float64 = 0.2,
                             weight_inv_power::Float64 = 0.5,
                             min_t::Float64 = 0.1)::@NamedTuple{
                                precursor_idx::String,
                                lloq::Float64,
                                aic_mean::Float64,
                                mean_intercept::Float64,
                                intercept_q95::Float64,
                                max_at_100::Float64
                             } where {T,U<:AbstractFloat}
    

    #Design Matrix for 5 parameter model
    #uniform cubic B-spline basis 
    X_5p = Float64.(UniformSplineDMAT(x_data, 3, 3)[:,1:5])

    #Design Matrix for Linear Model
    X_lin = hcat(
                    ones(T, length(x_data)),
                    x_data
        )
    
    #Initialize

    #predictions 
    ŷ_5p = zeros(T, length(x_data));
    ŷ_lin = zeros(T, length(x_data));
    #weighted residuals 
    w_r_5p = zeros(T, length(x_data));
    w_r_lin = zeros(T, length( x_data));
    
    #resampled data placeholder 
    x_r = zeros(eltype(x_data), length(x_data));
    y_r = zeros(eltype(x_data), length(x_data));
    W = getWeights(x_data, weight_inv_power, min_t);
    weights = diag(W);

    ############
    #Given the design matrix we ned to solve an Ax = y problem. 
    #We are resampling y at the same points so A can stay the same.
    #Since we are solving the same system many times with new data, use LU decomposition 

    #Set up Linear systems for  5p 
    luA_5p = LinearAlgebra.lu(X_5p'*W*X_5p)
    XtW_5p = X_5p'*W
    β_5p, F_5p = zeros(T, 5), zeros(T, 5)
    #return "test"
    #Set up Linear systems for  2p 
    luA_lin = LinearAlgebra.lu(X_lin'*W*X_lin)
    XtW_lin = X_lin'*W
    β_lin, F_lin = zeros(T, 2), zeros(T, 2)

    #Begin Bootstrapping
    fill!(bootstrap_Y, zero(T))
    fill!(bootstrap_intercept, zero(T))
    for col_idx in range(1, size(bootstrap_Y, 2))
        #Resample at each point on the calibration curve
        resampleData!(x_r, y_r, x_data, y_data)
        #return β_5p, ŷ_5p, X_5p, XWX_5p, F_5p, XtW_5p, y_r
        solveβevalŷ!(β_5p, ŷ_5p, X_5p, luA_5p, F_5p, XtW_5p, y_r)
        solveβevalŷ!(β_lin, ŷ_lin, X_lin, luA_lin, F_lin, XtW_lin, y_r)
        
        fillWeightedRedisuals!(w_r_5p,weights,ŷ_5p,y_r)
        fillWeightedRedisuals!(w_r_lin,weights,ŷ_lin,y_r)
        #return x_r,ŷ_5p,y_r
        #println("β_5p $β_5p β_lin $β_lin")
        #break
        corrected_σ_5p = getCorrectedSigma(w_r_5p,x_data,15.0)
        corrected_σ_lin = getCorrectedSigma(w_r_lin,x_data,15.0)
        
        AIC_5p, AIC_lin = getAIC(w_r_5p, corrected_σ_5p, 5), getAIC(w_r_lin, corrected_σ_lin, 2)
        #Prefered model has lower AIC. So if the output is 1, the 5-parameter model was a better fit 
        aic_diffs[col_idx] = AIC_5p < AIC_lin ? 1.0 : -1.0
        if AIC_5p < AIC_lin
            #Predict over xbins
            predictXBins!(bootstrap_Y, col_idx, x_mat_bins_5p, β_5p)
            #Simulate data 
            sampleData!(bootstrap_Y, col_idx, w_bins, corrected_σ_5p)
            #Get Intercept 
            predictIntercept!(bootstrap_intercept, col_idx, x_mat_bins_5p, β_5p)

        else
            #Predict over xbins
            predictXBins!(bootstrap_Y, col_idx, x_mat_bins_lin, β_lin)
            #Simulate data 
            sampleData!(bootstrap_Y, col_idx, w_bins, corrected_σ_lin)
            #Get Intercept
            predictIntercept!(bootstrap_intercept, col_idx, x_mat_bins_lin, β_lin)
        end
    end
    
    
    #Get LLOQ
    LLOQ = 0.0
    LLOQ_bin_idx = 1
    intercept_mean = 0.0
    intercept_upper = 0.0
    max_at_100 = 0.0

    intercept_max_median = max(median(bootstrap_intercept), 0.0)#max(median(bootstrap_Y[1,:]), 0.0)
    intercept_q95 = quantile(bootstrap_intercept, 0.95)
    #intercept_mean = mean(bootstrap_Y[1,:])

    #Smallest 'x' where the predictive distribution has a coefficient of variation of less than 'max_cv'
    #Need to subtract the intercept from the mean before calculating CV. 
    for (i,row) in enumerate(eachrow(bootstrap_Y))
        row_mean = mean(row .- intercept_max_median)
        row_std = std((row.-intercept_max_median))
        cv = row_std/row_mean
        if (cv > max_cv) | ((row_mean - row_std)<intercept_q95) | (row_mean <= 0.0)
            LLOQ = x_bins[i]
            LLOQ_bin_idx = i
        end

        if i == size(bootstrap_Y, 1)
            max_at_100 = mean(row)
        end
    end
    
    return (precursor_idx = precursor_idx,
            lloq = max(LLOQ, first(x_data)), 
            aic_mean = mean(aic_diffs), 
            mean_intercept = median(bootstrap_intercept), 
            intercept_q95 = quantile(bootstrap_intercept, 0.95),
            max_at_100 = max_at_100)
end


function BootstrapGDF(gdf::GroupedDataFrame{DataFrame},
                        id_col::Symbol,
                        quant_col::Symbol;
                        n_bootstraps::Int = 1000,
                        n_xbins = 100,
                        xmin = 0,
                        xmax = 100,
                        power = 0.5,
                        min_t = 0.5)
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
    x_mat_bins_5p = Float64.(UniformSplineDMAT(x_bins, 3, 3))[:,1:5]
    x_mat_bins_lin = hcat(
                        ones(T, length(x_bins)),
                        x_bins
            )
    curves = combine(gdf) do psms
            if size(psms, 1)>=12
                try
                    BoostrapCalibration!(
                        psms[1,id_col],
                        bootstrap_Y,
                        bootstrap_intercept,
                        aic_diffs,
                        psms[!,:percent_light],
                        psms[!,quant_col],
                        x_mat_bins_5p,
                        x_mat_bins_lin,
                        w_bins,
                        x_bins,
                        weight_inv_power = power,
                        min_t = min_t,
                    )
                catch e
                    println("size(psms[!,:percent_light]): ", size(psms[!,:percent_light]))
                    println("psms[!,:percent_light]: ", psms[!,:percent_light])
                    rethrow(e)
                end
            else
                (precursor_idx = missing,
                                lloq = missing, 
                                aic_mean = missing, 
                                mean_intercept = missing, 
                                intercept_q95 = missing,
                                max_at_100 = missing)
            end
    end
    return curves
end