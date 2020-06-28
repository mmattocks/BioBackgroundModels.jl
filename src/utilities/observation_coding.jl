#KMER ORDER/SEQUENCE INTEGER CODING UTILITIES
#higher order DNA alphabet
struct CompoundAlphabet
  symbols::Dict{Mer{DNAAlphabet{2}},Integer}
  integers::Dict{Integer,Mer{DNAAlphabet{2}}}
  #build a CompoundAlphabet for DNA of some order_no
    function CompoundAlphabet(alphabet::Tuple, order_no::Integer)
        symbols = Dict{Mer{DNAAlphabet{2}},Integer}()
        integers = Dict{Integer,Mer{DNAAlphabet{2}}}()
    
        tuples = Array{Tuple}
        if order_no > 0
            tuples = collect(Iterators.product([alphabet for order in 0:order_no]...))
        else #0th order compound product of an alphabet is that alphabet
            tuples = collect(Iterators.product(alphabet))
        end
    
        @inbounds for index in eachindex(tuples)
            tuple_seq = Mer{DNAAlphabet{2}}(collect(tuples[index]))
            integers[index] = tuple_seq
            symbols[tuple_seq] = index
        end
    
        new(symbols,integers)
    end
end


struct N_Order_ntSequence
    alphabet::CompoundAlphabet
    seq_lengths::Vector{Integer}
    order_kmers::Vector{Vector{Mer{DNAAlphabet{2}}}}
  end
  

function code_job_obs(job_ids::Vector{Chain_ID}, obs_sets::Dict{String,Vector{LongSequence{DNAAlphabet{2}}}})
    code_jobs=Vector{Tuple{String,Integer}}()
    for id in job_ids #assemble a vector of observation encoding jobs
        code_job=(id.obs_id, id.order)
        !(code_job in code_jobs) && push!(code_jobs, code_job)
    end

    code_dict = Dict{Tuple{String,Integer}, AbstractArray}()
    @showprogress 1 "Encoding observations..." for (obs_id, order) in code_jobs #build the appropriate sample sets once
        order_seqs = get_order_n_seqs(obs_sets[obs_id],order) #get the kmer sequences at the appropriate order
        coded_seqs = code_seqs(order_seqs, sorted=true) #numerically code the sequences in trainable format
        code_dict[(obs_id, order)] = coded_seqs
    end
    return code_dict
end

#from a vector of LongSequences, get
function get_order_n_seqs(seqs::Vector{LongSequence{DNAAlphabet{2}}}, order_no::Integer, base_tuple::Tuple=ACGT)
    kmer_vecs = Vector{Vector{Mer{DNAAlphabet{2}}}}()
    length_vec = Vector{Integer}()
    window = order_no + 1

    for seq in seqs
        kmer_vec = Vector{Mer{DNAAlphabet{2}}}()
        @inbounds for (i, kmer) in collect(each(Mer{DNAAlphabet{2},window},seq))
            push!(kmer_vec, kmer)
        end

        push!(kmer_vecs, kmer_vec)
        push!(length_vec, length(seq))
    end

    return nordseqs = N_Order_ntSequence(CompoundAlphabet(base_tuple, order_no), length_vec, kmer_vecs)
end

#convert tuple kmers to symbol codes
function code_seqs(input::N_Order_ntSequence, offsets::Vector=[0 for i in 1:length(input.order_kmers)]; sorted::Bool=false)
    symbol_no=length(input.alphabet.symbols)
    if symbol_no <= typemax(UInt8)
        integer_type = UInt8
    elseif typemax(UInt8) < symbol_no < typemax(UInt16)
        integer_type = UInt16
    else
        integer_type = UInt32
    end

    alphabet = input.alphabet
    output = zeros(integer_type,  length(input.order_kmers), (maximum([length(seq) for seq in input.order_kmers])+1)) #leave 1 missing value after the longest sequence forindexing sequence length in CLHMM messages
    sorted && (sort_idxs = sortperm(input.seq_lengths,rev=true))
    sorted ? input_seqs = input.order_kmers[sort_idxs] : input_seqs = input.order_kmers

    for (i, seq) in enumerate(input_seqs)
        for t in 1:length(seq)
            curr_kmer = seq[t]
            curr_code = alphabet.symbols[curr_kmer]
            output[i,t+offsets[i]]=curr_code
        end
    end
    return output
end
