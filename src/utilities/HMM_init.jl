function autotransition_init(K::Integer, order::Integer)
        π0 = rand(Dirichlet(ones(K)/K)) #uninformative prior on initial state probabilities
        π = strong_autotrans_matrix(K)
        no_emission_symbols = Int(4^(order+1)) #alphabet size for the order
        emission_dists = [generate_emission_dist(no_emission_symbols) for i in 1:K]
        #generate the HMM with the appropriate transition matrix and emissions distributions
        hmm = HMM(π0, π, emission_dists)
end
            #function to construct HMM transition matrix with strong priors on auto-transition
            function strong_autotrans_matrix(states::Integer, prior_dope::AbstractFloat=(states*250.0), prior_background::AbstractFloat=.1)
                transition_matrix=zeros(states,states)
                for k in 1:states
                    dirichlet_params = fill(prior_background, states)
                    dirichlet_params[k] = prior_dope
                    transition_matrix[k,:] = rand(Dirichlet(dirichlet_params))
                end
                return transition_matrix
            end

            #function to construct HMM state emission distribution from uninformative dirichlet over the alphabet size
            function generate_emission_dist(no_emission_symbols, prior=Dirichlet(ones(no_emission_symbols)/no_emission_symbols))
                return Categorical(rand(prior))
            end

