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
    maskcol = [zeros(Int64,0,0) for i in 1:size(position_df,1)]
    position_df.MaskMatrix=maskcol

    @Threads.threads for entry in eachrow(position_df)
        scaffold = entry[columns[1]]
        maskLength = length(entry[columns[2]])
        seqStart = entry[columns[3]]

        scaffold_coords_dict = Dict{String,DataFrame}()
        
        for ((partition, part_scaffold), df) in partitioned_scaffolds
            if scaffold == part_scaffold
                scaffold_coords_dict[partition] = df
            end
        end

        entry.MaskMatrix=mask_sequence_by_partition(maskLength, seqStart, scaffold_coords_dict)
    end

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
                
                    return seqMask
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
                        throw(DomainError("Position $position not found among partition coordinates!"))
                    else
                        return position_partition_id, three_prime_extent, sample_strand
                    end
                end

#function to partition a genome into coordinate sets of:
#merged exons
#"periexonic" sequences (genes with 5' and 3' boundaries projected -/+perigenic_pad bp, minus exons) - includes promoter elements, introns, 3' elements
#intergenic sequences (everything else)
#given a valid gff3
function partition_genome_coordinates(gff3_path::String, perigenic_pad::Integer=500)
    # construct dataframes of scaffolds and metacoordinates
    scaffold_df = build_scaffold_df(gff3_path)

    #partition genome into intragenic, periexonic, and exonic coordinate sets
    #assemble exonic featureset
    exon_df = DataFrame(SeqID = String[], Start = Integer[], End = Integer[], Strand=Char[])
    build_feature_df!(gff3_path, "CDS", "MT", exon_df)

    #project exon coordinates onto the scaffold bitwise, merging overlapping features and returning the merged dataframe
    merged_exon_df = DataFrame(SeqID = String[], Start = Integer[], End = Integer[], Strand=Char[])
    @showprogress 1 "Partitioning exons..." for scaffold_subframe in DataFrames.groupby(exon_df, :SeqID) # for each scaffold subframe that has exon features
        scaffold_id  = scaffold_subframe.SeqID[1] #get the scaffold id
        scaffold_bitarray = init_scaffold_bitarray(scaffold_id, scaffold_df, false, stranded=true) #init a stranded bitarray of scaffold length
        project_features_to_bitarray!(scaffold_subframe,scaffold_bitarray) #project features as Trues on bitarray of falses
        merged_subframe = get_feature_df_from_bitarray(scaffold_id, scaffold_bitarray) #get a feature df from the projected bitarray
        append!(merged_exon_df,merged_subframe) #append the merged scaffold df to the overall merged df
    end

    #assemble gene featureset
    gene_df =  DataFrame(SeqID = String[], Start = Integer[], End = Integer[], Strand = Char[])
    build_feature_df!(gff3_path, "gene", "MT", gene_df)

    perigenic_pad > 0 && add_pad_to_coordinates!(gene_df, scaffold_df, perigenic_pad) #if a perigenic pad is specified (to capture promoter/downstream elements etc in the periexonic set), apply it to the gene coords

    #build intergenic coordinate set by subtracting gene features from scaffold bitarrays
    intergenic_df = DataFrame(SeqID = String[], Start = Integer[], End = Integer[])
    @showprogress 1 "Partitioning intergenic regions..." for scaffold_subframe in DataFrames.groupby(scaffold_df, :SeqID) # for each scaffold subframe that has exon features
        scaffold_id  = scaffold_subframe.SeqID[1] #get the scaffold id
        scaffold_bitarray = init_scaffold_bitarray(scaffold_id, scaffold_df, true) #init a bitarray of scaffold length
        if any(isequal(scaffold_id),gene_df.SeqID) #if any genes on the scafold
            scaffold_genes=gene_df[findall(isequal(scaffold_id), gene_df.SeqID), :] #get the gene rows by finding the scaffold_id
            subtract_features_from_bitarray!(scaffold_genes,scaffold_bitarray) #subtract the gene positions from the scaffold bitarray
        end
        intragenic_subframe = get_feature_df_from_bitarray(scaffold_id, scaffold_bitarray) #get a feature df from the projected bitarray
        append!(intergenic_df,intragenic_subframe) #append the merged scaffold df to the overall merged df
    end

    #build periexonic set by projecting gene coordinates onto the scaffold bitwise, then subtracting the exons
    periexonic_df = DataFrame(SeqID = String[], Start = Integer[], End = Integer[], Strand=Char[])
    @showprogress 1 "Partitioning periexonic regions..." for gene_subframe in DataFrames.groupby(gene_df, :SeqID) # for each scaffold subframe that has exon features
        scaffold_id  = gene_subframe.SeqID[1] #get the scaffold id
        scaffold_bitarray = init_scaffold_bitarray(scaffold_id, scaffold_df, false, stranded=true) #init a bitarray of scaffold length
        project_features_to_bitarray!(gene_subframe,scaffold_bitarray) #project features as Trues on bitarray of falses
        if any(isequal(scaffold_id),merged_exon_df.SeqID)
            scaffold_exons=merged_exon_df[findall(isequal(scaffold_id), merged_exon_df.SeqID), :]
            subtract_features_from_bitarray!(scaffold_exons,scaffold_bitarray)
        end
        periexonic_subframe = get_feature_df_from_bitarray(scaffold_id, scaffold_bitarray) #get a feature df from the projected bitarray
        append!(periexonic_df,periexonic_subframe) #append the merged scaffold df to the overall merged df
    end

    return Dict("exon"=>merged_exon_df, "periexonic"=>periexonic_df, "intergenic"=>intergenic_df)
end

            #partition_genome_coordinates() SUBFUNCTIONS
            #BITARRAY SCAFFOLD REPRESENTATION SUBFUNCTIONS
            function init_scaffold_bitarray(scaffold_id::String, scaffold_df, value::Bool; stranded::Bool=false)
                scaffold_length = scaffold_df.End[findfirst(isequal(scaffold_id), scaffold_df.SeqID)]
                if stranded
                    value ? (return rscaffold_bitarray = trues(scaffold_length,2)) : (return rscaffold_bitarray = falses(scaffold_length,2))
                else
                    value ? (return rscaffold_bitarray = trues(scaffold_length,1)) : (return rscaffold_bitarray = falses(scaffold_length,1))
                end
            end

            function project_features_to_bitarray!(scaffold_feature_sf::SubDataFrame, scaffold_bitarray::BitArray)
                @inbounds for item in eachrow(scaffold_feature_sf)
                    scaffold_bitarray[item.Start:item.End,1] = [true for base in item.Start:item.End]
                    if size(scaffold_bitarray)[2] == 2 #if the bitarray is stranded
                        if item.Strand == '+'
                            scaffold_bitarray[item.Start:item.End,2] = [true for base in item.Start:item.End]
                        end
                    end
                end
            end

            function subtract_features_from_bitarray!(scaffold_feature_sf::DataFrame, scaffold_bitarray::BitArray)
                @inbounds for item in eachrow(scaffold_feature_sf)
                    scaffold_bitarray[item.Start:item.End,1] = [false for base in item.Start:item.End]
                end
            end

            function get_feature_df_from_bitarray(scaffold_id::String, scaffold_bitarray::BitArray)
                size(scaffold_bitarray)[2] == 2 ? scaffold_feature_df = DataFrame(SeqID = String[], Start = Integer[], End = Integer[], Strand=Char[]) : scaffold_feature_df = DataFrame(SeqID = String[], Start = Integer[], End = Integer[])

                new_feature_start = findnext(view(scaffold_bitarray,:,1), 1)
                while new_feature_start != nothing # while new features are still found on the bitarray
                        if size(scaffold_bitarray)[2] == 2 #if stranded, get strand info
                            scaffold_bitarray[new_feature_start,2] == 1 ? new_feature_strand = '+' : new_feature_strand = '-'
                        end
                        if findnext(!eval,view(scaffold_bitarray,:,1),new_feature_start) != nothing
                            new_feature_end = findnext(!eval,view(scaffold_bitarray,:,1),new_feature_start)-1 #find next false after feature start and subtract 1 for feature end
                        else
                            new_feature_end = size(scaffold_bitarray)[1] #if none is found, the end of the feature is the end of hte scaffold
                        end
                        size(scaffold_bitarray)[2] == 2 ? push!(scaffold_feature_df,[scaffold_id, new_feature_start, new_feature_end, new_feature_strand]) : push!(scaffold_feature_df,[scaffold_id, new_feature_start, new_feature_end]) #push stranded feature info as appropriate
                        new_feature_start = findnext(view(scaffold_bitarray,:,1), new_feature_end+1)
                end
                return scaffold_feature_df
            end
