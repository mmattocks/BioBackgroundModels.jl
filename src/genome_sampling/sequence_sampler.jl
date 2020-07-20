####SAMPLING FUNCTIONS####
#function for a Distributed worker to produce a set of samples of given parameters from genomic sequences
function get_sample_set(input_sample_jobs::RemoteChannel, completed_sample_jobs::RemoteChannel, progress_updates::RemoteChannel)
    while isready(input_sample_jobs)
        genome_path, genome_index_path, partition_df, partitionid, sample_set_length, sample_window_min, sample_window_max, deterministic = take!(input_sample_jobs)

        stranded::Bool = get_strand_dict()[partitionid]
        scaffold_sequence_record_dict::Dict{String,LongSequence} = build_scaffold_seq_dict(genome_path, genome_index_path)

        sample_df = DataFrame(SampleScaffold=String[],SampleStart=Integer[],SampleEnd=Integer[],SampleSequence=LongSequence[],Strand=Char[])
        metacoordinate_bitarray = trues(partition_df.MetaEnd[end])
        sample_set_counter = 0

        while sample_set_counter < sample_set_length #while we don't yet have enough sample sequence
            sample_scaffold::String, sample_Start::Integer, sample_End::Integer, sample_metaStart::Integer, sample_metaEnd::Integer, sample_sequence::LongSequence, strand::Char = get_sample(metacoordinate_bitarray, sample_window_min, sample_window_max,  partition_df, scaffold_sequence_record_dict; stranded=stranded, deterministic=deterministic)
            push!(sample_df,[sample_scaffold, sample_Start, sample_End, sample_sequence, strand]) #push the sample to the df
            sample_length = sample_End - sample_Start + 1
            sample_set_counter += sample_length #increase the counter by the length of the sampled sequence
            metacoordinate_bitarray[sample_metaStart:sample_metaEnd] = [false for base in 1:sample_length] #mark these residues as sampled
            put!(progress_updates, (partitionid, min(sample_set_counter,sample_set_length)))
        end
        put!(completed_sample_jobs,(partitionid, sample_df))
    end
end
                #get_sample_set() SUBFUNCTIONS
                #function defining whether partitions respect stranding upon fetching sequence (ie is the sequence fetched in the feature strand orientation, or are we agnostic about the strand we sample?)
                function get_strand_dict()
                    return Dict("exon"=>true,"periexonic"=>true,"intergenic"=>false)
                end
                #function to obtain a dict of scaffold sequences from a FASTA reader
                function build_scaffold_seq_dict(genome_fa, genome_index)
                    genome_reader = open(FASTA.Reader, genome_fa, index=genome_index)
                    seq_dict::Dict{String, LongSequence} = Dict{String,FASTA.Record}()
                    @inbounds for record in genome_reader
                        id = FASTA.identifier(record)
                        seq_dict[rectify_identifier(id)]=FASTA.sequence(record)
                    end
                    close(genome_reader)
                    return seq_dict
                end

                # function to convert scaffold ID from that observed by the masked .fna to the more legible one observed by the GRCz11 GFF3
                function rectify_identifier(scaffold_id::String)
                    if length(scaffold_id) >= 4 && scaffold_id[1:4] == "CM00" #marks chromosome scaffold
                        chr_code = scaffold_id[5:10]
                        chr_no = "$(Int((parse(Float64,chr_code)) - 2884.2))"
                        return chr_no
                    else
                        return scaffold_id
                    end
                end

                #function to produce a single sample from a metacoordinate set and the feature df
                function get_sample(metacoordinate_bitarray::BitArray, sample_window_min::Integer, sample_window_max::Integer, partition_df::DataFrame,  scaffold_seq_dict::Dict{String,LongSequence}; stranded::Bool=false, deterministic::Bool=false)
                    proposal_acceptance = false
                    sample_metaStart = 0
                    sample_metaEnd = 0
                    sample_Start = 0
                    sample_End = 0
                    sample_sequence = LongSequence{DNAAlphabet{2}}("")
                    sample_scaffold = ""
                    strand = nothing
                    
                    while proposal_acceptance == false
                        available_indices = findall(metacoordinate_bitarray) #find all unsampled indices
                        window = "FAIL"
                        start_index,feature_metaStart,feature_metaEnd, feature_length = 0,0,0,0
                        strand = '0'
                        while window == "FAIL"
                            start_index = rand(available_indices) #randomly choose from the unsampled indices
                            feature_metaStart, feature_metaEnd, strand = get_feature_params_from_metacoord(start_index, partition_df, stranded)
                            feature_length = length(feature_metaStart:feature_metaEnd) #find the metaboundaries of the feature the index occurs in
                            if feature_length >= sample_window_min #don't bother finding windows on features smaller than min
                                window = determine_sample_window(feature_metaStart, feature_metaEnd, start_index, metacoordinate_bitarray, sample_window_min, sample_window_max) #get an appropriate sampling window around the selected index, given the feature boundaries and params
                            end
                        end

                        sample_scaffoldid, sample_scaffold_start, sample_scaffold_end = meta_to_feature_coord(window[1],window[2],partition_df)

                        proposal_sequence = fetch_sequence(sample_scaffoldid, scaffold_seq_dict, sample_scaffold_start, sample_scaffold_end, strand; deterministic=deterministic) #get the sequence associated with the sample window

                        if mask_check(proposal_sequence)
                            proposal_acceptance = true #if the sequence passes the mask check, accept the proposed sample
                            sample_scaffold = sample_scaffoldid
                            sample_Start = sample_scaffold_start
                            sample_End = sample_scaffold_end
                            sample_metaStart = window[1]
                            sample_metaEnd = window[2]
                            sample_sequence=proposal_sequence
                        end
                    end
                    return sample_scaffold, sample_Start, sample_End, sample_metaStart, sample_metaEnd, sample_sequence, strand
                end
                
                #function to find a valid sampling window
                function determine_sample_window(feature_metaStart::Integer, feature_metaEnd::Integer, metacoord::Integer, metacoord_bitarray::BitArray, sample_window_min::Integer, sample_window_max::Integer)
                    window_start = 0
                    window_end = 0
                    feature_size = feature_metaEnd - feature_metaStart + 1
                    feature_sampled = any(!eval, metacoord_bitarray[feature_metaStart:feature_metaEnd])
                    #attempt to construct the sample window by taking the whole feature
                    if !feature_sampled && feature_size < sample_window_max && feature_size > sample_window_min
                        window_start = feature_metaStart
                        window_end = feature_metaEnd
                        return window_start, window_end
                    else #if this fails, build the biggest window we can from the sampling point, within the feature
                        featurepos = metacoord - feature_metaStart + 1
                        next_sampled_index = findnext(!eval, metacoord_bitarray[feature_metaStart:feature_metaEnd], featurepos)
                        prev_sampled_index = findprev(!eval, metacoord_bitarray[feature_metaStart:feature_metaEnd], featurepos)
                        next_sampled_index === nothing ? (window_end = feature_metaEnd) : (window_end = next_sampled_index + feature_metaStart - 1)
                        prev_sampled_index === nothing ? (window_start = feature_metaStart) : (window_start = prev_sampled_index + feature_metaStart - 1)
                        windowsize = window_end - window_start + 1
                        #check to see if this window is bigger than the min, if not, return a failure code
                        windowsize < sample_window_min && return "FAIL"
                        #check to see if this window is bigger than the max, if so, trim it before returning, removing as evenly as possible around the metacoordinate
                        if windowsize > sample_window_max
                            bases_to_trim = windowsize - sample_window_max
                            clearance_5P = metacoord - window_start + 1
                            clearance_3P = window_end - metacoord + 1
                            trimmed=0
                            while trimmed < bases_to_trim
                                if clearance_5P >= clearance_3P && clearance_5P > 0
                                    clearance_5P -= 1
                                elseif clearance_3P > clearance_5P && clearance_3P > 0
                                    clearance_3P -= 1
                                end
                                trimmed +=1
                            end 
                            window_start = metacoord - clearance_5P + 1
                            window_end = metacoord + clearance_3P - 1
                        end
                        return window_start, window_end
                    end
                end


                # function to check for repetitive stretches or degenerate bases in proposal sequence
                function mask_check(proposal_sequence::LongSequence)
                    proposal_acceptance = true
                    if hasambiguity(proposal_sequence) || isrepetitive(proposal_sequence, (length(proposal_sequence) ÷ 10))
                        proposal_acceptance = false
                    end
                    return proposal_acceptance
                end

####SHARED SEQUENCE FETCHER####
# function to get proposal sequence from dict of scaffold sequences, given coords and scaffold id
function fetch_sequence(scaffold_id::String, scaffold_seq_dict::Dict{String, LongSequence}, proposal_start::Integer, proposal_end::Integer, strand::Char; deterministic=false)
    if strand == '0' #unstranded samples may be returned with no preference in either orientation
        deterministic ? strand = '+' : (rand(1)[1] <=.5 ? strand = '+' : strand = '-')
    end
    if strand == '+'
        proposal_sequence = scaffold_seq_dict[scaffold_id][proposal_start:proposal_end]
    elseif strand == '-'
        proposal_sequence = reverse_complement(scaffold_seq_dict[scaffold_id][proposal_start:proposal_end])
    else
        throw(ArgumentError("Invalid sample code! Must be '+', '-', or '0' (random strand)"))
    end

    return proposal_sequence
end

####SHARED BASIC COORDINATE SUBFUNCTIONS####
# function to push scaffold ID, start, and end points of given featuretype to supplied dataframe
function build_feature_df!(GFF3_path::String, feature_type::String, scaffold_exclusion::String, feature_df::DataFrame)
    reader = open(GFF3.Reader, GFF3_path) # access the GFF3
    @inbounds for record in reader # iterate over Gff3 records
        if GFF3.isfeature(record) # if the record is a feature
            if GFF3.featuretype(record) == feature_type #if the features is of the requested type, get the following info
                seqID = GFF3.seqid(record)
                if seqID != scaffold_exclusion
                    seq_start = GFF3.seqstart(record)
                    seq_end = GFF3.seqend(record)
                    if feature_type == "CDS" || feature_type == "gene"
                        seq_strand = convert(Char, GFF3.strand(record))
                        push!(feature_df, [seqID, seq_start, seq_end, seq_strand])
                    else
                        push!(feature_df, [seqID, seq_start, seq_end]) # push relevant info to the df
                    end
                end
            end
        end
    end
    close(reader)
end

#function to assemble dataframe of scaffold coords + metacoords given gff3
function build_scaffold_df(gff3_path)
    scaffold_df = DataFrame(SeqID = String[], Start = Integer[], End = Integer[])
    build_feature_df!(gff3_path, "supercontig", "MT", scaffold_df)
    build_feature_df!(gff3_path, "chromosome", "MT", scaffold_df)
    add_metacoordinates!(scaffold_df)
    return scaffold_df
end

#function to add pad to either side of some featureset
function add_pad_to_coordinates!(feature_df::DataFrame, scaffold_df::DataFrame, pad_size::Integer; col_symbols::Array{Symbol}=[:Start, :End])
    pad_start_array = zeros(Integer,size(feature_df,1))
    pad_end_array = zeros(Integer,size(feature_df,1))
    feature_df[!, col_symbols[1]] = [max(feature_df.Start[i]-pad_size,1) for i in 1:size(feature_df,1)] #truncate pads at beginning and end of scaffolds
    feature_df[!, col_symbols[2]] = [min(feature_df.End[i]+pad_size,scaffold_df.End[findfirst(isequal(feature_df.SeqID[i]),scaffold_df.SeqID)]) for i in 1:size(feature_df,1)]
end

####SHARED METACOORDINATE FUNCTIONS####
# function to add a metacoordinate column to a dataframe of scaffold positions, allowing sampling across all scaffolds
function add_metacoordinates!(feature_df::DataFrame)
    meta_start_array = Integer[]
    meta_end_array = Integer[]
    metaposition = 1
    @inbounds for feature in eachrow(feature_df)
        push!(meta_start_array, (metaposition))
        metaposition += feature.End - feature.Start
        push!(meta_end_array, (metaposition))
        metaposition += 1
    end
    feature_df.MetaStart = meta_start_array # metacoordinate contains 'start' metaposition across all genomic material
    feature_df.MetaEnd = meta_end_array # metacoordinate contains 'end' metaposition across all genomic material
end

# function to convert metacoordinate set to scaffold-relative coordinates
function meta_to_feature_coord(meta_start::Integer, meta_end::Integer, feature_df::DataFrame)
    feature_row = get_feature_row_index(feature_df, meta_start)
    seqid = feature_df.SeqID[feature_row]
    scaffold_start = feature_df.Start[feature_row] + (meta_start - feature_df.MetaStart[feature_row])
    scaffold_end = feature_df.End[feature_row] - (feature_df.MetaEnd[feature_row] - meta_end)
    return seqid, scaffold_start, scaffold_end
end

#function to obtain the feature boundaries and strand of the feature that a metacoordinate falls within
function get_feature_params_from_metacoord(metacoordinate::Integer, feature_df::DataFrame, stranded::Bool)
    feature_row = get_feature_row_index(feature_df, metacoordinate)
    feature_metaStart = feature_df.MetaStart[feature_row]
    feature_metaEnd = feature_df.MetaEnd[feature_row]
    stranded ? feature_strand = feature_df.Strand[feature_row] : feature_strand = '0'
    return feature_metaStart, feature_metaEnd, feature_strand
end

#function obtain the index of a feature given its metacoordinate
function get_feature_row_index(feature_df::DataFrame, metacoordinate::Integer)
    total_feature_bases = feature_df.MetaEnd[end]
    if metacoordinate == total_feature_bases #edge case of metacoord at end of range
        return size(feature_df,1) #last index in df
    else
        index = findfirst(end_coord->end_coord>metacoordinate,feature_df.MetaEnd) #find the index of the feature whose end metacoord is > the query metacoord
        if index > 1 && metacoordinate == feature_df.MetaEnd[index-1] #if the metacoordinate is the last base of the previous feature
            if metacoordinate >= feature_df.MetaStart[index-1] && metacoordinate <= feature_df.MetaEnd[index-1] #confirm that the metacoord is in the previous feature
                return index-1 #return the previous feature index
            end
        elseif metacoordinate >= feature_df.MetaStart[index] && metacoordinate <= feature_df.MetaEnd[index] #else confirm that the metacoordinate is in the found feature
            return index #return the found feature index
        else
            throw(DomainError("Unexpected metacoordinate $metacoordinate in partition of $total_feature_bases bases, with feature start $(feature_df.MetaStart[index]), end $(feature_df.MetaEnd[index])"))
        end
    end
end