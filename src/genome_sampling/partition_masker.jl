function get_partition_code_dict(dict_forward::Bool=true)
    if dict_forward == true
        return partition_code_dict = Dict("intergenic"=>1,"periexonic"=>2,"exon"=>3)
    else
        return code_partition_dict = Dict(1=>"intergenic", 2=>"periexonic", 3=>"exon")
    end
end

function make_padded_df(position_fasta::String, gff3_path::String, genome_path::String, genome_index_path::String, pad::Integer)
    position_reader = FASTA.Reader(open((position_fasta),"r"))
    genome_reader = open(FASTA.Reader, genome_path, index=genome_index_path)
    scaffold_df = build_scaffold_df(gff3_path)
    position_df = DataFrame(SeqID = String[], Start=Int[], End=Int[], PadSeq = LongSequence[], PadStart=Int[], RelStart=Int[], SeqOffset=Int[])
    scaffold_seq_dict = build_scaffold_seq_dict(genome_path, genome_index_path)

    for entry in position_reader
        scaffold = FASTA.identifier(entry)

        if scaffold != "MT"
            desc_array = split(FASTA.description(entry))
            pos_start = parse(Int, desc_array[2])
            pos_end = parse(Int, desc_array[4])
            scaffold_end = scaffold_df.End[findfirst(isequal(scaffold), scaffold_df.SeqID)]

            pad_start=max(1,pos_start-pad)
            pad_length= pos_start - pad_start
            seq_offset = pad - pad_length
            padded_seq = fetch_sequence(scaffold, scaffold_seq_dict, pad_start, pos_end, '+')

            if !hasambiguity(padded_seq)
                push!(position_df, [scaffold, pos_start, pos_end, padded_seq, pad_start, pad_length, seq_offset])
            end
        end
    end
    
    close(position_reader)
    close(genome_reader)
    return position_df
end

function add_partition_masks!(position_df::DataFrame, gff3_path::String, perigenic_pad::Integer=500, columns::Tuple{Symbol,Symbol,Symbol}=(:SeqID, :PadSeq, :PadStart))
    partitions=["exon", "periexonic", "intergenic"]
    partition_coords_dict = partition_genome_coordinates(gff3_path, perigenic_pad)
    partitioned_scaffolds = divide_partitions_by_scaffold(partition_coords_dict)
    maskcol = Vector{Matrix{Integer}}()

    @showprogress 1 "Masking..." for entry in eachrow(position_df)
        scaffold = entry[columns[1]]
        maskLength = length(entry[columns[2]])
        seqStart = entry[columns[3]]

        scaffold_coords_dict = Dict{String,DataFrame}()
        
        for ((partition, part_scaffold), df) in partitioned_scaffolds
            if scaffold == part_scaffold
                scaffold_coords_dict[partition] = df
            end
        end

        push!(maskcol, mask_sequence_by_partition(maskLength, seqStart, scaffold_coords_dict))
    end

    position_df.MaskMatrix=maskcol
end
                #add_partition_masks!() SUBFUNCTIONS
                function divide_partitions_by_scaffold(partition_coords_dict::Dict{String, DataFrame})
                    scaffold_coords_dict = Dict{Tuple{String, String}, DataFrame}()
                    for (partition_id, partition_df) in partition_coords_dict
                        for scaffold_subframe in groupby(partition_df, :SeqID)
                            scaffold_id = scaffold_subframe.SeqID[1]
                            scaffold_df = copy(scaffold_subframe)
                            scaffold_coords_dict[(partition_id, scaffold_id)] = scaffold_df
                        end
                    end
                    return scaffold_coords_dict
                end

                function mask_sequence_by_partition(maskLength::Integer, seqStart::Integer, scaffold_coords_dict::Dict{String, DataFrame})
                    partition_code_dict = get_partition_code_dict()
                    seqMask = zeros(Integer, (maskLength, 2))
                    position = seqStart
                    while position <= seqStart+maskLength
                        position_partition, partition_extent, position_strand = find_position_partition(position, scaffold_coords_dict)
                
                        partition_code = partition_code_dict[position_partition]
                        mask_position = position - seqStart + 1
                        seqMask[mask_position:min(maskLength,mask_position + partition_extent),1] .= partition_code
                        if position_strand == '+'
                            seqMask[mask_position:min(maskLength,mask_position + partition_extent),2] .= 1
                        elseif position_strand == '-'
                            seqMask[mask_position:min(maskLength,mask_position + partition_extent),2] .= -1
                        else
                            seqMask[mask_position:min(maskLength,mask_position + partition_extent),2] .= 0
                        end
                
                        position += partition_extent + 1
                    end
                
                    if any(iszero, seqMask[:,1])
                        @error "Malformed seqmask!!"
                    else
                        return seqMask
                    end
                end

                function find_position_partition(position::Integer, partition_dict::Dict{String, DataFrame})
                    foundPos = false
                    position_partition_id = ""
                    three_prime_extent = 0
                    sample_strand = 0
                    for (partition_id, partition) in partition_dict
                        hitindex = findfirst(x->x>=position, partition.End)
                        if hitindex != nothing
                            if position >= partition.Start[hitindex]
                                foundPos=true
                                position_partition_id = partition_id
                                three_prime_extent = partition.End[hitindex] - position
                                if partition_id == "exon" || partition_id == "periexonic"
                                    sample_strand = partition.Strand[hitindex]
                                end
                            end
                        end
                    end
                
                    if foundPos == false
                        @error "Position $position not found among partition coordinates!"
                    else
                        return position_partition_id, three_prime_extent, sample_strand
                    end
                end