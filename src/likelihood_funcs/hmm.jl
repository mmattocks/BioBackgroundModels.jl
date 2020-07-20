function obs_lh_given_hmm(observations, hmm; linear=true)
    linear ? (return linear_likelihood(observations, hmm)) : (return bw_likelihood(Matrix(transpose(observations)), hmm))
end

function bw_likelihood(observations, hmm)
    obs_lengths = [findfirst(iszero,observations[:,o])-1 for o in 1:size(observations,2)]

    lls = bw_llhs(hmm, observations)
    log_α = messages_forwards_log(hmm.a, hmm.A, lls, obs_lengths)
    log_β = messages_backwards_log(hmm.A, lls, obs_lengths)

    O = size(lls)[3] #the last T value is the 0 end marker of the longest T
    log_pobs = zeros(O)
    Threads.@threads for o in 1:O
        log_pobs[o] = logsumexp(lps.(log_α[:,1,o], log_β[:,1,o]))
    end

    return sum(log_pobs)
end

function linear_likelihood(observations,hmm)
    O= size(observations,1);
    obs_lengths = [findfirst(iszero,observations[o,:])-1 for o in 1:size(observations,1)] #mask calculations here rather than mle_step to prevent recalculation every iterate

    a = transpose(log.(hmm.a)); A = log.(hmm.A)
    N = length(hmm.B)
    mask=observations.!=0
    #INITIALIZATION
    βoi_T = zeros(O,N); βoi_t = zeros(O,N) #log betas at T initialised as zeros

    #RECURRENCE
    βoi_T = backwards_lh_sweep!(hmm, A, N, βoi_T, βoi_t, observations, mask, obs_lengths)

    #TERMINATION
    lls = c_llhs(hmm,observations[:,1])
    α1om = lls .+ a #first position forward msgs

    return lps([logsumexp(lps.(α1om[o,:], βoi_T[o,:])) for o in 1:O])
end
                #LINEAR_STEP SUBFUNCS
                function backwards_lh_sweep!(hmm, A, N, βoi_T, βoi_t, observations, mask, obs_lengths)
                    for t in maximum(obs_lengths)-1:-1:1
                        last_β=copy(βoi_T)
                        lls = c_llhs(hmm,observations[:,t+1])
                        omask = findall(mask[:,t+1])
                        βoi_T[omask,:] .+= view(lls,omask,:)
                        Threads.@threads for m in 1:N
                            βoi_t[omask,m] = logsumexp.(eachrow(view(βoi_T,omask,:).+transpose(view(A,m,:))))
                        end
                        βoi_T=copy(βoi_t)
                    end
                    return βoi_T
                end

