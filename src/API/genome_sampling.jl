"""
    setup_sample_jobs(genome_path, genome_index_path, gff3_path, sample_set_length, sample_window_min, sample_window_max, perigenic_pad; deterministic)

Return sample job channels for 'execute_sample_jobs', given 'genome_path' to a properly formatted FASTA genome, 'genome_index_path' to the associated .FAI, 'gff3_path' to the genome's GFF3 feature database, a total number of bases to sample, 'sample_set_length', as well as 'sample_window_min' and 'sample_window_max' lengths for individual samples. 'perigenic_pad' specifies the number of bases up- and down-stream of exonic sequences to be considered 'perigenic'; if 0, the perigenic partition will be entirely intronic. Optional boolean 'deterministic' may be set to 'true' to always return '+' strand sequences from unstranded samples; otherwise unstranded samples may be returned with either orientation.
"""

#function to partition genome and set up Distributed RemoteChannels so partitions can be sampled simultaneously
function setup_sample_jobs(genome_path::String, genome_index_path::String, gff3_path::String, sample_set_length::Integer, sample_window_min::Integer, sample_window_max::Integer, perigenic_pad::Integer; deterministic::Bool=false)
    #argument checking
    !ispath(genome_path) && throw(ArgumentError("Bad genome path!"))
    !ispath(genome_index_path) && throw(ArgumentError("Bad genome index path!"))
    !ispath(gff3_path) && throw(ArgumentError("Bad gff3 path!"))
    sample_set_length < 1 && throw(ArgumentError("Sample set length must be a positive integer!"))
    sample_window_min < 1 || sample_window_max < 1 && throw(ArgumentError("Sample window minimum and maximum bounds must be positive integers!"))
    sample_window_min >= sample_window_max && throw(ArgumentError("Sample window minimum size must be smaller than maximum size"))
    perigenic_pad < 0 && throw(ArgumentError("Perigenic pad must be 0 or positive!"))

    coordinate_partitions = partition_genome_coordinates(gff3_path, perigenic_pad)
    sample_sets = DataFrame[]
    input_sample_jobs = RemoteChannel(()->Channel{Tuple}(length(coordinate_partitions))) #channel to hold sampling jobs
    completed_sample_jobs = RemoteChannel(()->Channel{Tuple}(length(coordinate_partitions))) #channel to hold completed sample dfs
    for (partitionid, partition) in coordinate_partitions
        add_metacoordinates!(partition)
        put!(input_sample_jobs, (genome_path, genome_index_path, partition, partitionid, sample_set_length, sample_window_min, sample_window_max, deterministic))
    end
    progress_channel = RemoteChannel(()->Channel{Tuple}(20))
    return (input_sample_jobs, completed_sample_jobs, progress_channel, sample_set_length)
end

"""
    execute_sample_jobs!(channels, worker_pool)

Return sample dataframes for each genomic partition, given appropriately set up job channels (from setup_sample_jobs!) and a vector of worker ids ('worker_pool', [1] for master).
"""

function execute_sample_jobs(channels::Tuple{RemoteChannel,RemoteChannel,RemoteChannel,Integer}, worker_pool::Vector{Int64}; partitions::Integer=3)
    input_sample_channel, completed_sample_channel, progress_channel, sample_set_length = channels

    #argument checking
    length(worker_pool) < 1 && throw(ArgumentError("Worker pool must contain one or more worker IDs!"))
    !isready(input_sample_channel) && throw(ArgumentError("Input sample channel not ready, likely channel set from setup_sample_jobs passed incorrectly"))
    sample_set_length < 1 && throw(ArgumentError("Sample set length must be a positive integer, likely channel set from setup_sample_jobs passed incorrectly"))

    #send sampling jobs to workers
    if isready(input_sample_channel) > 0
        @info "Sampling.."
        #WORKERS SAMPLE
        for (n, worker) in enumerate(worker_pool)
            if n <= partitions #no more workers than partitions to be used
                remote_do(get_sample_set, worker, input_sample_channel, completed_sample_channel, progress_channel)
            end
        end
    else
        @error "No sampling jobs!"
    end

    #progress meters for sampling
    sampling_meters=Dict{String, Progress}()
    overall_sampling_meter=Progress(partitions,"Overall sampling progress:")
    completed_counter = 0
    ProgressMeter.update!(overall_sampling_meter, completed_counter)
    sampling_offset = ones(Bool, partitions)

    #collect progress updates while waiting on completion of sampling jobs
    while completed_counter < partitions
        wait(progress_channel)
        partition_id, progress = take!(progress_channel)
        if haskey(sampling_meters, partition_id)
            ProgressMeter.update!(sampling_meters[partition_id], progress)
        else
            offset = findfirst(sampling_offset)[1]
            sampling_meters[partition_id] = Progress(sample_set_length, "Sampling partition $partition_id:", offset)
            ProgressMeter.update!(sampling_meters[partition_id], progress)
            sampling_offset[offset] = false
        end
        if progress == sample_set_length
            completed_counter += 1
            ProgressMeter.update!(overall_sampling_meter, completed_counter)
        end
    end

    #collect sample dfs by partition id when ready
    sample_record_dfs = Dict{String,DataFrame}()
    collected_counter = 0
    while collected_counter < partitions
        wait(completed_sample_channel)
        partition_id, sample_df = take!(completed_sample_channel)
        sample_record_dfs[partition_id] = sample_df
        collected_counter += 1
        @info "Partition $partition_id completed sampling..."
    end

    return sample_record_dfs
end
