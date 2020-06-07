#FUNCTIONS TO ANALYSE LEARNT HMMS
#returns likelihood of appropriately coded observation set given hmm
function test_hmm(hmm::HMM, test_set, order)
    order_seqs = get_order_n_seqs(test_set,order) #get the kmer sequences at the appropriate order
    coded_seqs = code_seqs(order_seqs) #numerically code the sequences in trainable format
    return lin_obs_set_lh(hmm, coded_seqs) #function in CLHMM - linear.jl
end

function get_diagonal_array(hmm::HMM)
    k = length(hmm.π0)
    if k > 1
        diagonal = zeros(k)
        for i in 1:k
            diagonal[i] = hmm.π[i,i]
        end
        return diagonal
    else
        return [NaN]
    end
end

function chain_diagonal_stability_matrix(chain::Vector{Any})
    diagonals=zeros(length(chain),length(chain[1][2].π0))
    for (s,step) in enumerate(chain)
        diagonals[s,:]=get_diagonal_array(step[2])        
    end
    return diagonals
end

#function to simulate run lengths for vector of diagonal values
function sim_run_lengths(diagonal_value::AbstractArray, samples::Integer)
    if isnan(diagonal_value[1])
        return [NaN]
    else
        mean_run_lengths = zeros(length(diagonal_value))
        for (i, value) in enumerate(diagonal_value)
            runlengths = zeros(Integer, samples)
            for s in 1:samples
                run = true
                runlength = 0
                while run
                    runlength += 1
                    if rand(1)[1] > value
                        run = false
                    end
                end
                runlengths[s] = runlength
            end
            mean_run_lengths[i] = mean(runlengths)
        end
        return mean_run_lengths
    end
end

#function to produce matrix of hmm chain parameter coordinate evolution- selects 3 matched emission parameters from 3+ state HMMs, one per state for 3 states. Matching is performed by minimizing euclidean distance between the state emission probability vectors D for the states to be matched (as EM-optimised replicates will have similar states organized in different orders)

function chain_3devo_coords(chains::Vector{Vector{Any}})
    length(chains)<=0 && throw(ArgumentError, "Argument must be vector of more than one hmm chain")
    coords=[Vector{Tuple{AbstractFloat,AbstractFloat,AbstractFloat}}() for i in 1:length(chains)]

    for step in chains[1]
        hmm=step[2]
        length(hmm.D)<3 && throw(ArgumentError, "3- or greater state hmms required")
        push!(coords[1], (hmm.D[1].p[1],hmm.D[2].p[1],hmm.D[3].p[1]))
    end

    for (c, chain) in enumerate(chains[2:end])
        ref_idxs=Vector{Integer}()
        ref_vecs=[chains[1][end][2].D[1].p,chains[1][end][2].D[2].p,chains[1][end][2].D[3].p]
        end_hmm=chain[end][2]
        length(end_hmm.D)<3 && throw(ArgumentError, "3- or greater state hmms required")

        for vec in ref_vecs
            euclideans=Vector{AbstractFloat}()
            for d in 1:length(end_hmm.D)
                push!(euclideans, euclidean(vec, end_hmm.D[d].p))
            end
            min_euc_idx = findmin(euclideans)[2]
            while min_euc_idx in ref_idxs
                euclideans[min_euc_idx]=1.
                min_euc_idx = findmin(euclideans)[2]
            end
            push!(ref_idxs,findmin(euclideans)[2])
        end

        for step in chain
            hmm=step[2]
            push!(coords[c+1], (hmm.D[ref_idxs[1]].p[1],hmm.D[ref_idxs[2]].p[1],hmm.D[ref_idxs[3]].p[1]))
        end
    end

    return coords
end