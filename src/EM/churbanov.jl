function linear_step(hmm, observations, obs_lengths)
    O,T = size(observations);
    a = transpose(log.(hmm.a)); A = log.(hmm.A) #less expensive to transpose transmatrix in the backwards_sweep than to transpose here and take a ranged view ie view(A,m:m,:) is more expensive than transpose(view(A,m,:))
    N = length(hmm.B); Γ = length(hmm.B[1].support);
    mask=observations.!=0
    #INITIALIZATION
    βoi_T = zeros(O,N); βoi_t = zeros(O,N) #log betas at T initialised as zeros
    Eoγim_T = fill(-Inf,O,Γ,N,N); Eoγim_t = fill(-Inf,O,Γ,N,N)
    for m in 1:N, i in 1:N, γ in 1:Γ, o in 1:O
        observations[o, obs_lengths[o]] == γ && m == i && (Eoγim_T[o, γ, i, m] = 0)
    end
    Tijm_T = fill(-Inf,O,N,N,N); Tijm_t = fill(-Inf,O,N,N,N) #Ti,j(T,m) = 0 for all m; in logspace
        
    #RECURRENCE
    βoi_T,Tijm_T,Eoγim_T=backwards_sweep!(hmm,A,N,Γ,βoi_T,βoi_t,Tijm_T,Tijm_t,Eoγim_T,Eoγim_t,observations,mask,obs_lengths)
        
    #TERMINATION
    lls = c_llhs(hmm,observations[:,1])
    α1om = lps.(lls,a) #first position forward msgs
    log_pobs = [c_lse(lps.(α1om[o,:], βoi_T[o,:])) for o in 1:O]

    Toij = [c_lse([lps(Tijm_T[o,i,j,m], α1om[o,m], -log_pobs[o]) for m in 1:N]) for o in 1:O, i in 1:N, j in 1:N] #terminate Tijs with forward messages
    Eoiγ = [c_lse([lps(Eoγim_T[o,γ,i,m], α1om[o,m], -log_pobs[o]) for m in 1:N]) for o in 1:O, i in 1:N, γ in 1:Γ] #terminate Eids with forward messages

    #INTEGRATE ACROSS OBSERVATIONS AND SOLVE FOR NEW BHMM PARAMS
    #INITIAL STATE DIST
    a_o=α1om.+βoi_T.-c_lse.(eachrow(α1om.+βoi_T)) #estimate a for each o
    obs_penalty=log(O) #broadcast subtraction to normalise log prob vals by obs number
    new_a=c_lse.(eachcol(a_o)).-obs_penalty #sum over obs and normalise by number of obs
    #TRANSITION MATRIX
    new_A = [c_lse(view(Toij,:,i,j)) for i in 1:N, j in 1:N].-c_lse.(eachcol(c_lse.([view(Toij,o,i,:) for o in 1:O, i in 1:N])))
    #EMISSION MATRIX
    new_b=[c_lse(view(Eoiγ,:,i,γ)) for i in 1:N, γ in 1:Γ].-c_lse.(eachcol(c_lse.([view(Eoiγ,o,i,:) for o in 1:O, i in 1:N])))
    new_B=[Categorical(exp.(new_b[i,:])) for i in 1:N]

    return BHMM(exp.(new_a), exp.(new_A), new_B, hmm.partition), lps(log_pobs)
end
                #LINEAR_STEP SUBFUNCS
                function backwards_sweep!(hmm, A, N, Γ, βoi_T, βoi_t, Tijm_T, Tijm_t, Eoγim_T, Eoγim_t, observations, mask, obs_lengths)
                    for t in maximum(obs_lengths)-1:-1:1
                        last_β=copy(βoi_T)
                        lls = c_llhs(hmm,observations[:,t+1])
                        omask = findall(mask[:,t+1])
                        βoi_T[omask,:] .+= view(lls,omask,:)
                        Threads.@threads for m in 1:N
                            βoi_t[omask,m] = c_lse.(eachrow(view(βoi_T,omask,:).+transpose(view(A,m,:))))
                            for j in 1:N, i in 1:N
                                Tijm_t[omask, i, j, m] .= c_lse.(eachrow(lps.(view(Tijm_T,omask,i,j,:), view(lls,omask,:), transpose(view(A,m,:)))))
                                i==m && (Tijm_t[omask, i, j, m] .= logaddexp.(Tijm_t[omask, i, j, m], lps.(last_β[omask,j], A[m,j],lls[omask,j])))
                            end
                            for i in 1:N, γ in 1:Γ
                                Eoγim_t[omask, γ, i, m] .= c_lse.(eachrow(lps.(view(Eoγim_T,omask,γ,i,:),view(lls,omask,:),transpose(view(A,m,:)))))
                                if i==m
                                    symmask = findall(observations[:,t].==γ)
                                    Eoγim_t[symmask, γ, i, m] .= logaddexp.(Eoγim_t[symmask, γ, i, m], βoi_t[symmask,m])
                                end
                            end
                        end
                        βoi_T=copy(βoi_t); Tijm_T=copy(Tijm_t); Eoγim_T = copy(Eoγim_t);
                    end
                    return βoi_T, Tijm_T, Eoγim_T
                end

                #logsumexp with single dispatch for speed
                function c_lse(X::AbstractArray{T}; dims=:) where {T<:Real}
                    u=maximum(X)
                    u isa AbstractArray || isfinite(u) || return float(u)
                    return u + log(sum(x -> exp(x-u), X))
                end

                #log likelihoods
                function c_llhs(hmm, observation)
                    lls = zeros(length(observation),length(hmm.B))
                    Threads.@threads for d in 1:length(hmm.B)
                        lls[:,d] = logpdf.(hmm.B[d], observation)
                    end
                    return lls
                end

