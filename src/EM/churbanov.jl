function linear_step(hmm, observations, obs_lengths)
    O,T = size(observations);
    a = log.(hmm.π); π0 = transpose(log.(hmm.π0))
    N = length(hmm.D); Γ = length(hmm.D[1].support);
    mask=observations.!=0
    #INITIALIZATION
    βoi_T = zeros(O,N); βoi_t = zeros(O,N) #log betas at T initialised as zeros
    Eoγim_T = fill(-Inf,O,Γ,N,N); Eoγim_t = fill(-Inf,O,Γ,N,N)
    for m in 1:N, i in 1:N, γ in 1:Γ, o in 1:O
        observations[o, obs_lengths[o]] == γ && m == i && (Eoγim_T[o, γ, i, m] = 0)
    end
    Tijm_T = fill(-Inf,O,N,N,N); Tijm_t = fill(-Inf,O,N,N,N) #Ti,j(T,m) = 0 for all m; in logspace
        
    #RECURRENCE
    βoi_T,Tijm_T,Eoγim_T=backwards_sweep!(hmm,a,N,Γ,βoi_T,βoi_t,Tijm_T,Tijm_t,Eoγim_T,Eoγim_t,observations,mask,obs_lengths)
        
    #TERMINATION
    lls = c_llhs(hmm,observations[:,1])
    α1om = lls .+ π0 #first position forward msgs

    Toij = [c_lse([lps(Tijm_T[o,i,j,m], α1om[o,m]) for m in 1:N]) for o in 1:O, i in 1:N, j in 1:N] #terminate Tijs with forward messages

    Eoiγ = [c_lse([lps(Eoγim_T[o,γ,i,m], α1om[o,m]) for m in 1:N]) for o in 1:O, i in 1:N, γ in 1:Γ] #terminate Eids with forward messages

    #INTEGRATE ACROSS OBSERVATIONS AND SOLVE FOR NEW HMM PARAMS
    obs_penalty=log(O) #broadcast subtraction to normalise log prob vals by obs number
    #INITIAL STATE DIST
    π0_o=α1om.+βoi_T.-c_lse.(eachrow(α1om.+βoi_T)) #estimate π0 for each o
    new_π0=c_lse.(eachcol(π0_o)).-obs_penalty #sum over obs and normalise by number of obs
    #TRANSITION MATRIX
    a_int = Toij.-c_lse.([view(Toij,o,i,:) for o in 1:O, i in 1:N])
    new_a = c_lse.([a_int[:,i,j] for i in 1:N, j in 1:N]).-obs_penalty
    #EMISSION MATRIX
    e_int=Eoiγ.-c_lse.([view(Eoiγ,o,j,:) for o in 1:O, j in 1:N])
    new_b=c_lse.([view(e_int,:,j,γ) for γ in 1:Γ, j in 1:N]).-obs_penalty
    new_D=[Categorical(exp.(new_b[:,i])) for i in 1:N]

    return typeof(hmm)(exp.(new_π0), exp.(new_a), new_D), lps([c_lse(lps.(α1om[o,:], βoi_T[o,:])) for o in 1:O])
end
                #LINEAR_STEP SUBFUNCS
                function backwards_sweep!(hmm, a, N, Γ, βoi_T, βoi_t, Tijm_T, Tijm_t, Eoγim_T, Eoγim_t, observations, mask, obs_lengths)
                    for t in maximum(obs_lengths)-1:-1:1
                        last_β=copy(βoi_T)
                        lls = c_llhs(hmm,observations[:,t+1])
                        omask = findall(mask[:,t+1])
                        βoi_T[omask,:] .+= view(lls,omask,:)
                        Threads.@threads for m in 1:N
                            βoi_t[omask,m] = c_lse.(eachrow(view(βoi_T,omask,:).+transpose(view(a,m,:))))
                            for j in 1:N, i in 1:N
                                Tijm_t[omask, i, j, m] .= c_lse.(eachrow(lps.(view(Tijm_T,omask,i,j,:), view(lls,omask,:), transpose(view(a,m,:)))))
                                i==m && (Tijm_t[omask, i, j, m] .= logaddexp.(Tijm_t[omask, i, j, m], (last_β[omask,j].+ a[m,j].+ lls[omask,j])))
                            end
                            for i in 1:N, γ in 1:Γ
                                Eoγim_t[omask, γ, i, m] .= c_lse.(eachrow(lps.(view(Eoγim_T,omask,γ,i,:),view(lls,omask,:),transpose(view(a,m,:)))))
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

                #logsumexp
                function c_lse(X::AbstractArray{T}; dims=:) where {T<:Real}
                    u=maximum(X)
                    u isa AbstractArray || isfinite(u) || return float(u)
                    u + log(sum(x -> exp(x-u), X))
                end

                #log likelihoods
                function c_llhs(hmm, observation)
                    lls = zeros(length(observation),length(hmm.D))
                    Threads.@threads for d in 1:length(hmm.D)
                        lls[:,d] = logpdf.(hmm.D[d], observation)
                    end
                    return lls
                end

