#function to obtain positional likelihoods for a sequence under a given HMM.
function get_BGHMM_symbol_lh(seq::AbstractMatrix, hmm::AbstractHMM)
    @assert size(seq)[1] == 1
    (seq=Array(transpose(seq))) # one sequence at a time only
    symbol_lhs = zeros(length(seq))
    length_mask = [length(seq)-1]
    
    lls = bw_llhs(hmm, seq) #obtain log likelihoods for sequences and states
    log_α = messages_forwards_log(hmm.π0, hmm.π, lls, length_mask) #get forward messages
    log_β = messages_backwards_log(hmm.π, lls, length_mask) # get backwards messages

    #calculate observation probability and γ weights
    K,Tmaxplus1,Strand = size(lls) #the last T value is the 0 end marker of the longest T

    #transforms to cut down log_ξ, log_γ assignment times
    lls = permutedims(lls, [2,1,3]) # from (K,T) to (T,K)
    log_α = permutedims(log_α, [2,1,3])
    log_β = permutedims(log_β, [2,1,3])

    log_γ = fill(-Inf, Tmaxplus1,K)
    log_pobs = logsumexp(lps.(log_α[1,:], log_β[1,:]))

    @inbounds for i = 1:K, t = 1:length_mask[1]
            log_γ[t,i] = lps(log_α[t,i],log_β[t,i],-log_pobs)
    end

    for t in 1:length_mask[1]
        symbol_lh::AbstractFloat = -Inf #ie log(p=0)
        for k = 1:K #iterate over states
                state_symbol_lh::AbstractFloat = lps(log_γ[t,k], log(hmm.D[k].p[seq[t]])) #state symbol likelihood is the γ weight * the state symbol probability (log implementation)
                symbol_lh = logaddexp(symbol_lh, state_symbol_lh) #sum the probabilities over states
        end
        symbol_lhs[t] = symbol_lh
    end

    return symbol_lhs[1:end-1] #remove trailing index position
end

#function to calculate BGHMM from an observation set and a dict of BGHMMs
function BGHMM_likelihood_calc(observations::DataFrame, BGHMM_dict::Dict, code_partition_dict = get_partition_code_dict(false); symbol=:PadSeq)
    lh_matrix_size = ((findmax(length.(collect(values(observations[!, symbol]))))[1]), length(observations[!, symbol]))
    BGHMM_lh_matrix = zeros(lh_matrix_size) #T, Strand, O

    BGHMM_fragments = fragment_observations_by_BGHMM(observations[!, symbol], observations.MaskMatrix)

    @showprogress 1 "Writing frags to matrix.." for (jobid, frag) in BGHMM_fragments
        (frag_start, o, partition, strand) = jobid

        partition_BGHMM::HMM = BGHMM_dict[code_partition_dict[partition]][1]
        no_symbols = length(partition_BGHMM.D[1].p)
        order = Int(log(4,no_symbols) - 1)

        order_seq = get_order_n_seqs([frag], order)
        coded_seq = code_seqs(order_seq)
        
        subseq_symbol_lh = get_BGHMM_symbol_lh(coded_seq, partition_BGHMM)

        if strand == -1
            subseq_symbol_lh = reverse(subseq_symbol_lh)
        end #positive, unstranded frags  are inserted as-is
        BGHMM_lh_matrix[frag_start:frag_start+length(frag)-1,o] = subseq_symbol_lh
    end

    return BGHMM_lh_matrix
end

function fragment_observations_by_BGHMM(seqs::AbstractVector, masks::AbstractVector)
    likelihood_jobs = Vector{Tuple{Tuple,BioSequence{DNAAlphabet{4}}}}()
    @showprogress 1 "Fragmenting observations by partition..." for (o, obs_seq) in enumerate(seqs)
        mask = masks[o]
        frags = Vector{BioSequence{DNAAlphabet{4}}}() #container for all subsequences in observation
        frag = BioSequence{DNAAlphabet{4}}()

        frag_end=0
        frag_start = 1

        while frag_start < length(obs_seq) # while we're not at the sequence end
            curr_partition = mask[frag_start,1] #get the partition code of the frag start
            curr_strand = mask[frag_start,2] #get the strand of the frag start

            #JOBID COMPOSED HERE
            jobid = (frag_start, o, curr_partition, curr_strand) #compose an identifying index for this frag

            findnext(!isequal(curr_partition),mask[:,1],frag_start) != nothing ? frag_end = findnext(!isequal(curr_partition),mask[:,1],frag_start) -1 : frag_end = length(obs_seq) #find the next position in the frag that has a different partition mask value from hte current one and set that position-1 to frag end, alternately frag end is end of the overall  sequence 
            frag = obs_seq[frag_start:frag_end] #get the frag bases
            if curr_strand == -1 #if the fragment is reverse stranded
                frag = reverse_complement(frag) #use the reverse complement sequence
            end

            push!(likelihood_jobs,(jobid, frag)) #put the frag in the jobs vec
            frag_start = frag_end + 1 #move on
        end        
    end
    return likelihood_jobs
end