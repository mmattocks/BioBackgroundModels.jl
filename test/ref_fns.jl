#slow algo written in mimicry of Churbanov & Winters
function old_linear(hmm, observations, obs_lengths)
    O = size(observations)[2]
    a = log.(hmm.A); a = log.(hmm.a)
    N = length(hmm.B); B = length(hmm.B[1].support); b = [log(hmm.B[m].p[γ]) for m in 1:N, γ in 1:B]
    α1oi = zeros(O,N); β1oi = zeros(O,N); Eoi = zeros(O,N,B); Toij = zeros(O,N,N); Aoi = zeros(O,N); log_pobs=zeros(O); γt=0
    
    for o in 1:O
        #INITIALIZATION
        T = obs_lengths[o]; βoi_T = zeros(N) #log betas at T initialised as zeros
        EiT = fill(-Inf,B,N,N); TijT = fill(-Inf,N,N,N) #Ti,j(T,m) = 0 for all m; in logspace
        
        @inbounds for m in 1:N, i in 1:N, γ in 1:B
            observations[T, o] == γ && m == i ? (EiT[γ, i, m] = 0) :
                (EiT[γ, i, m] = -Inf) #log Ei initialisation
        end

        #RECURRENCE
        @inbounds for t in T-1:-1:1
            βoi_t = similar(βoi_T); Tijt = similar(TijT); Eit = similar(EiT)
            Γ = observations[t+1,o]; γt = observations[t,o]
            Γ==0 && println("hooo!! $o $t")
            for m in 1:N
                βoi_t[m] = logsumexp([lps(a[m,j], b[j,Γ], βoi_T[j]) for j in 1:N])
                for i in 1:N
                    for j in 1:N
                        Tijt[i, j, m] = logsumexp([lps(a[m,n], TijT[i, j, n], b[n,Γ]) for n in 1:N])
                        i==m && (Tijt[i, j, m] = logaddexp(Tijt[i, j, m], lps(βoi_T[j], a[m,j], b[j, Γ])))
                    end
                    for γ in 1:B
                        Eit[γ, i, m] = logsumexp([lps(b[n,Γ], a[m,n], EiT[γ, i, n]) for n in 1:N])
                        i==m && γ==γt && (Eit[γ, i, m] = logaddexp(Eit[γ, i, m], βoi_t[m]))
                    end
                end
            end
            βoi_T = βoi_t; TijT = Tijt; EiT = Eit
        end

        #TERMINATION
        Γ = observations[1,o]
        α1oi[o,:] = [lps(a[i], b[i, Γ]) for i in 1:N]
        β1oi[o,:] = βoi_T
        log_pobs[o] = logsumexp(lps.(α1oi[o,:], βoi_T[:]))
        Eoi[o,:,:] = [logsumexp([lps(EiT[γ,i,m], a[m], b[m,γt]) for m in 1:N]) for i in 1:N, γ in 1:B]
        Toij[o,:,:] = [logsumexp([lps(TijT[i,j,m], a[m], b[m,γt]) for m in 1:N]) for i in 1:N, j in 1:N]
    end

    #INTEGRATE ACROSS OBSERVATIONS AND SOLVE FOR NEW BHMM PARAMS
    obs_penalty=log(O)
    a_o=α1oi.+β1oi.-logsumexp.(eachrow(α1oi.+β1oi))
    new_a=logsumexp.(eachcol(a_o)).-obs_penalty

    a_int = Toij.-logsumexp.([Toij[o,i,:] for o in 1:O, i in 1:N])
    new_a = logsumexp.([a_int[:,i,j] for i in 1:N, j in 1:N]).-obs_penalty
    
    e_int=Eoi.-logsumexp.([Eoi[o,j,:] for o in 1:O, j in 1:N])
    new_b=logsumexp.([e_int[:,j,d] for d in 1:B, j in 1:N]).-obs_penalty
    new_D=[Categorical(exp.(new_b[:,i])) for i in 1:N]

    return typeof(hmm)(exp.(new_a), exp.(new_a), new_D), lps(log_pobs)
end


#HMMBase
function mouchet_mle_step(hmm::BHMM, observations) where F
    # NOTE: This function works but there is room for improvement.

    log_likelihoods = mouchet_log_likelihoods(hmm, observations)

    log_α = mouchet_messages_forwards_log(hmm.a, hmm.A, log_likelihoods)
    log_β = mouchet_messages_backwards_log(hmm.A, log_likelihoods)
    log_A = log.(hmm.A)

    normalizer = logsumexp(log_α[1,:] + log_β[1,:])

    # E-step

    T, K = size(log_likelihoods)
    log_ξ = zeros(T, K, K)

    @inbounds for t = 1:T-1, i = 1:K, j = 1:K
        log_ξ[t,i,j] = log_α[t,i] + log_A[i,j] + log_β[t+1,j] + log_likelihoods[t+1,j] - normalizer
    end

    ξ = exp.(log_ξ)
    ξ ./= sum(ξ, dims=[2,3])

    # M-step
    new_A = sum(ξ[1:end-1,:,:], dims=1)[1,:,:] #index fix-MM
    new_A ./= sum(new_A, dims=2)

    new_a = exp.((log_α[1,:] + log_β[1,:]) .- normalizer)
    new_a ./= sum(new_a)

    # TODO: Cleanup/optimize this part
    γ = exp.((log_α .+ log_β) .- normalizer)

    B = Categorical[]
    for (i, d) in enumerate(hmm.B)
        # Super hacky...
        # https://github.com/JuliaStats/Distributions.jl/issues/809
        push!(B, fit_mle(Categorical, permutedims(observations), γ[:,i]))
    end

    typeof(hmm)(new_a, new_A, B), normalizer
end


function mouchet_log_likelihoods(hmm, observations)
    hcat(map(d -> logpdf.(d, observations), hmm.B)...)
end



function mouchet_messages_forwards_log(init_distn, trans_matrix, log_likelihoods)
    # OPTIMIZE
    log_alphas = zeros(size(log_likelihoods))
    log_trans_matrix = log.(trans_matrix)
    log_alphas[1,:] = log.(init_distn) .+ log_likelihoods[1,:]
    @inbounds for t = 2:size(log_alphas)[1]
        for i in 1:size(log_alphas)[2]
            log_alphas[t,i] = logsumexp(log_alphas[t-1,:] .+ log_trans_matrix[:,i]) + log_likelihoods[t,i]
        end
    end
    log_alphas
end

function mouchet_messages_backwards_log(trans_matrix, log_likelihoods)
    # OPTIMIZE
    log_betas = zeros(size(log_likelihoods))
    log_trans_matrix = log.(trans_matrix)
    @inbounds for t = size(log_betas)[1]-1:-1:1
        tmp = view(log_betas, t+1, :) .+ view(log_likelihoods, t+1, :)
        @inbounds for i in 1:size(log_betas)[2]
            log_betas[t,i] = logsumexp(view(log_trans_matrix, i, :) .+ tmp)
        end
    end
    log_betas
end