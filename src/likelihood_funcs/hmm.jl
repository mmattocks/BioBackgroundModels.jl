function obs_lh_given_hmm(observations, hmm; linear=true)
    linear ? (return linear_likelihood(observations, hmm)) : (return bw_likelihood(Matrix(transpose(observations)), hmm))
end

function bw_likelihood(observations, hmm)
    obs_lengths = [findfirst(iszero,observations[:,o])-1 for o in 1:size(observations)[2]]

    lls = bw_llhs(hmm, observations)
    log_α = messages_forwards_log(hmm.π0, hmm.π, lls, obs_lengths)
    log_β = messages_backwards_log(hmm.π, lls, obs_lengths)

    O = size(lls)[3] #the last T value is the 0 end marker of the longest T
    log_pobs = zeros(O)
    Threads.@threads for o in 1:O
        log_pobs[o] = logsumexp(lps.(log_α[:,1,o], log_β[:,1,o]))
    end

    return sum(log_pobs)
end

function linear_likelihood(observations, hmm)
    O = size(observations)[1]; obs_lengths = [findfirst(iszero,observations[o,:])-1 for o in 1:O]
    a = log.(hmm.π); π0 = log.(hmm.π0)
    N = length(hmm.D); D = length(hmm.D[1].support); b = [log(hmm.D[m].p[γ]) for m in 1:N, γ in 1:D]
    α1oi = zeros(O,N); β1oi = zeros(O,N); Eoi = zeros(O,D,N); Toij = zeros(O,N,N); πoi = zeros(O,N); log_pobs=zeros(O); γt=0
    
    Threads.@threads for o in 1:O
        #INITIALIZATION
        T = obs_lengths[o]; βT = zeros(N) #log betas at T initialised as zeros
        #RECURRENCE
        for t in T-1:-1:1
            βt = similar(βT); Γ = observations[o,t+1]
            for m in 1:N
                βt[m] = logsumexp([lps(a[m,j], b[j,Γ], βT[j]) for j in 1:N])
            end
            βT = βt
        end
        #TERMINATION
        Γ = observations[o,1]
        α1oi[o,:] = [lps(π0[i], b[i, Γ]) for i in 1:N]
        log_pobs[o] = logsumexp(lps.(α1oi[o,:], βT[:]))
        end

    return logsumexp(log_pobs)
end