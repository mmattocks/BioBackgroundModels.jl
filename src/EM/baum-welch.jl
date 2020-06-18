"""
    ms_mle_step(hmm::AbstractHMM{F}, observations) where F

Perform one step of the EM (Baum-Welch) algorithm.

# Example
```julia
hmm, log_likelihood = ms_mle_step(hmm, observations)
```
"""

function bw_step(hmm::AbstractHMM{F}, observations, obs_lengths) where F
    lls = bw_llhs(hmm, observations)
    log_α = messages_forwards_log(hmm.π0, hmm.π, lls, obs_lengths)
    log_β = messages_backwards_log(hmm.π, lls, obs_lengths)
    log_π = log.(hmm.π)

    K,Tmaxplus1,O = size(lls) #the last T value is the 0 end marker of the longest T

    #transforms to cut down log_ξ, log_γ assignment times
    lls = permutedims(lls, [2,3,1]) # from (K,T,O) to (T,O,K)
    log_α = permutedims(log_α, [2,3,1])
    log_β = permutedims(log_β, [2,3,1])

    # E-step
    log_ξ = fill(-Inf, Tmaxplus1,O,K,K)
    log_γ = fill(-Inf, Tmaxplus1,O,K)
    log_pobs = zeros(O)
 
   @inbounds for o = 1:O
        log_pobs[o] = logsumexp(lps.(log_α[1,o,:], log_β[1,o,:]))
    end

    Threads.@threads for idx in [(i,j,o) for i=1:K, j=1:K, o=1:O]
        i,j,o = idx[1],idx[2],idx[3]
        obsl = obs_lengths[o]
        @inbounds for t = 1:obsl-1 #log_ξ & log_γ calculated to T-1 for each o
           log_ξ[t,o,i,j] = lps(log_α[t,o,i], log_π[i,j], log_β[t+1,o,j], lls[t+1,o,j], -log_pobs[o])
           log_γ[t,o,i] = lps(log_α[t,o,i], log_β[t,o,i], -log_pobs[o])
        end
        t=obsl #log_ξ @ T = 0
        log_ξ[t,o,i,j] = 0
        log_γ[t,o,i] = lps(log_α[t,o,i], log_β[t,o,i], -log_pobs[o])
    end

    ξ = similar(log_ξ)
    ξ .= exp.(log_ξ)
    ∑k_ξ = sum(ξ, dims=[3,4])
    nan_mask = ∑k_ξ .== 0
    ∑k_ξ[nan_mask] .= Inf #prevent NaNs in dummy renorm
    ξ  ./=  ∑k_ξ #dummy renorm across K to keep numerical creep from causing isprobvec to fail on new new_π during hmm creation

    γ = similar(log_γ)
    γ .= exp.(log_γ)
    ∑k_γ = sum(γ, dims=3)
    ∑k_γ[nan_mask[:,:,:]] .= Inf #prevent NaNs in dummy renorm
    γ ./= ∑k_γ #dummy renorm

    # M-step
    new_π = zeros(K,K)
    for i=1:K, j=1:K
        ∑otξ_vec = zeros(O)
        ∑otγ_vec = zeros(O)
        Threads.@threads for o in 1:O
            ∑otξ_vec[o] = sum(ξ[1:obs_lengths[o]-1,o,i,j])
            ∑otγ_vec[o] = sum(γ[1:obs_lengths[o]-1,o,i])
        end
        new_π[i,j] = sum(∑otξ_vec) / sum(∑otγ_vec)
    end
    new_π ./= sum(new_π, dims=2) #dummy renorm
    new_π0 = (sum(γ[1,:,:], dims=1)./sum(sum(γ[1,:,:], dims=1)))[1,:]
    new_π0 ./= sum(new_π0) #dummy renorm

    obs_mask = .!nan_mask
    obs_collection = observations[obs_mask[:,:]]

    D = Distribution{F}[]
    @inbounds for (i, d) in enumerate(hmm.D)
        γ_d = γ[:,:,i]
        push!(D, t_categorical_mle(Categorical, d.support[end], obs_collection, γ_d[obs_mask[:,:]]))
    end

    return typeof(hmm)(new_π0, new_π, D), lps(log_pobs)
end
                function bw_llhs(hmm, observations)
                    lls = zeros(length(hmm.D), size(observations)...)
                    Threads.@threads for d in 1:length(hmm.D)
                        (lls[d,:,:] = logpdf.(hmm.D[d], observations))
                    end
                    return lls
                end
                #Multisequence competent log implementations of forward and backwards algos
                function messages_forwards_log(init_distn, trans_matrix, log_likelihoods, obs_lengths)
                    log_alphas = zeros(size(log_likelihoods))
                    log_trans_matrix = log.(trans_matrix)
                    log_alphas[:,1,:] = log.(init_distn) .+ log_likelihoods[:,1,:]
                    Threads.@threads for o in 1:size(log_likelihoods)[3]
                        @inbounds for t in 2:obs_lengths[o], j in 1:size(log_likelihoods)[1]
                            log_alphas[j,t,o]=logsumexp(log_alphas[:,t-1,o] .+ log_trans_matrix[:,j]) + log_likelihoods[j,t,o]
                        end   
                    end
                    return log_alphas
                end

                function messages_backwards_log(trans_matrix, log_likelihoods, obs_lengths)
                    log_betas = zeros(size(log_likelihoods))
                    log_trans_matrix = log.(trans_matrix)
                    Threads.@threads for o in 1:size(log_likelihoods)[3]
                        @inbounds for t in obs_lengths[o]-1:-1:1
                            tmp = view(log_betas, :, t+1, o) .+ view(log_likelihoods, :, t+1, o)
                            @inbounds for i in 1:size(log_likelihoods)[1]
                                log_betas[i,t,o] = logsumexp(view(log_trans_matrix, i,:) .+ tmp)
                            end
                        end
                    end
                    return log_betas
                end
                #subfunc derived from Distributions.jl categorical.jl fit_mle, threaded
                function t_categorical_mle(::Type{<:Categorical}, k::Integer, x::AbstractArray{T}, w::AbstractArray{F}) where T<:Integer where F<:AbstractFloat
                    Categorical(t_pnormalize!(t_add_categorical_counts!(zeros(k), x, w)))
                end

                t_pnormalize!(v::Vector) = (v ./= sum(v); v)

                function t_add_categorical_counts!(h::Vector{F}, x::AbstractArray{T}, w::AbstractArray{F}) where T<:Integer where F<:AbstractFloat
                    n = length(x)
                    if n != length(w)
                        throw(DimensionMismatch("Inconsistent array lengths."))
                    end
                    hlock = ReentrantLock()
                    Threads.@threads for i=1:n
                        @inbounds xi = x[i]
                        @inbounds wi = w[i]
                        lock(hlock)
                        h[xi] += wi   # cannot use @inbounds, as no guarantee that x[i] is in bound
                        unlock(hlock)
                    end
                    return h
                end
