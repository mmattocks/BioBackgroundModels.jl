"""
Utility functions for learning and using background genomic hidden markov models
"""

#function to split random sample dataframe into training and test sets (divide total sequence length by half)
function split_obs_sets(sample_dfs::Dict{String,DataFrame})
    training_sets = Dict{String,Vector{BioSequence{DNAAlphabet{4}}}}()
    test_sets = Dict{String,Vector{BioSequence{DNAAlphabet{4}}}}()

    for (partition_id, partition) in sample_dfs
        partition.sampleLength = (partition.SampleEnd - partition.SampleStart) .+ 1
        midway = sum(partition.sampleLength)÷2
        split_index = 0
        counter = 0
        while split_index == 0
            counter += 1
            length_sum = sum(partition.sampleLength[1:counter])
            if length_sum > midway
                split_index = counter
            end
        end

        training_sets[partition_id]  = partition.SampleSequence[1:split_index-1]
        test_sets[partition_id] = partition.SampleSequence[split_index:end]
    end
    return training_sets, test_sets
end

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


    #function to determine required jobids for global search from hmms learnt in survey
    function HMM_global_search_params(hmm_results_dict::Dict)
        params_dict=Dict{String,Tuple{Integer,Integer}}()
        for ((partition, K, order, rep), chain) in hmm_results_dict
            if !haskey(params_dict, partition)
                params_dict[partition]=(K, order)
            else
                params_dict[partition]!=(K,order) && throw(ArgumentError, "More than one state number or order for partition")
            end
        end
        return params_dict
    end

    #function to add additional replicates to HMMs learnt in survey for estimation of global optimum
    function HMM_global_search_setup!(hmm_results_dict::Dict, params_dict::Dict{String,Tuple{Integer,Integer}}, search_replicates::Integer, search_thresh::AbstractFloat,  input_hmms::RemoteChannel, training_sets::Dict{String,Vector{BioSequence{DNAAlphabet{4}}}}, base_alphabet_size::Integer)
        no_input_hmms=search_replicates * length(params_dict)
        code_dict = Dict{String, Vector}()

        @showprogress 1 "Encoding observations..." for (partition_id, partition) in training_sets #build the appropriate sample sets once
            K,order_no=params_dict[partition_id]
            order_seqs = get_order_n_seqs(partition,order_no) #get the kmer sequences at the appropriate order
            coded_seqs = code_seqs(order_seqs, sorted=true) #numerically code the sequences in trainable format
            code_dict[partition_id] = coded_seqs
        end

        @showprogress 1 "Setting up HMMs..." for partition_id in keys(params_dict)
            K,order=params_dict[partition_id]
            for replicate in 1:search_replicates
                jobid=(partition_id, K, order, replicate)
                if haskey(hmm_results_dict, jobid) && length(hmm_results_dict[jobid]) > 0 #true if resuming from incomplete chain
                    chain=hmm_results_dict[jobid]
                    job_convergence=chain[end][5]
                    if job_convergence && chain[end][4] > search_thresh #if marked as converged from earlier search, remark for deeper search
                        hmm_results_dict[jobid][end][5]=0
                    end
                    
                    if !job_convergence #push the last hmm iterate for nonconverged chains to the input channel
                        iterate = chain[end][1]
                        hmm =  chain[end][2]
                        last_norm = chain[end][3]
                        put!(input_hmms, (jobid, iterate, hmm, last_norm, code_dict[partition_id]))
                    else #skip any jobs that have converged from previous runs
                        no_input_hmms -= 1
                    end

                else #initialise first HMM in chain
                    π0 = rand(Dirichlet(ones(K)/K)) #uninformative prior on initial state probabilities
                    π = generate_transition_matrix(K)
                    no_emission_symbols = Int(base_alphabet_size^(order+1)) #alphabet size for the order
                    emission_dists = [generate_emission_dist(no_emission_symbols) for i in 1:K]
                    #generate the HMM with the appropriate transition matrix and emissions distributions
                    hmm = CLHMM.HMM(π0, π, emission_dists)
                    hmm_results_dict[jobid] = [] #initialise the relevant results array
                    put!(input_hmms, (jobid, 1, hmm, 0.0, code_dict[partition_id]))
                end
            end
        end
        
        return no_input_hmms
    end
    