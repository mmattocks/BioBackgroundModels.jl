

#function to produce matrix of hmm chain parameter coordinate evolution- selects 3 matched emission parameters from 3+ state HMMs, one per state for 3 states. Matching is performed by minimizing euclidean distance between the state emission probability vectors B for the states to be matched (as EM-optimised replicates will have similar states organized in different orders)

function chain_3devo_coords(chains::Vector{Vector{Any}})
    length(chains)<=0 && throw(ArgumentError, "Argument must be vector of more than one hmm chain")
    coords=[Vector{Tuple{AbstractFloat,AbstractFloat,AbstractFloat}}() for i in 1:length(chains)]

    for step in chains[1]
        hmm=step[2]
        length(hmm.B)<3 && throw(ArgumentError, "3- or greater state hmms required")
        push!(coords[1], (hmm.B[1].p[1],hmm.B[2].p[1],hmm.B[3].p[1]))
    end

    for (c, chain) in enumerate(chains[2:end])
        ref_idxs=Vector{Integer}()
        ref_vecs=[chains[1][end][2].B[1].p,chains[1][end][2].B[2].p,chains[1][end][2].B[3].p]
        end_hmm=chain[end][2]
        length(end_hmm.B)<3 && throw(ArgumentError, "3- or greater state hmms required")

        for vec in ref_vecs
            euclideans=Vector{AbstractFloat}()
            for d in 1:length(end_hmm.B)
                push!(euclideans, euclidean(vec, end_hmm.B[d].p))
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
            push!(coords[c+1], (hmm.B[ref_idxs[1]].p[1],hmm.B[ref_idxs[2]].p[1],hmm.B[ref_idxs[3]].p[1]))
        end
    end

    return coords
end