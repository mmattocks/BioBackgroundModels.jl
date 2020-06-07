#function to partition genome and set up Distributed RemoteChannels so partitions can be sampled simultaneously
function setup_sample_jobs(genome_path::String, genome_index_path::String, gff3_path::String, sample_set_length::Integer, sample_window_min::Integer, sample_window_max::Integer, perigenic_pad::Integer; deterministic::Bool=false)
    coordinate_partitions = partition_genome_coordinates(gff3_path, perigenic_pad)
    sample_sets = DataFrame[]
    input_sample_jobs = RemoteChannel(()->Channel{Tuple}(length(coordinate_partitions))) #channel to hold sampling jobs
    completed_sample_jobs = RemoteChannel(()->Channel{Tuple}(length(coordinate_partitions))) #channel to hold completed sample dfs
    for (partitionid, partition) in coordinate_partitions
        add_metacoordinates!(partition)
        put!(input_sample_jobs, (genome_path, genome_index_path, partition, partitionid, sample_set_length, sample_window_min, sample_window_max, deterministic))
    end
    progress_channel = RemoteChannel(()->Channel{Tuple}(20))
    return (input_sample_jobs, completed_sample_jobs, progress_channel)
end

function execute_sample_jobs(channels::Tuple{RemoteChannel,RemoteChannel,RemoteChannel}, worker_pool::Vector{Int64}, sample_set_length::Integer; partitions::Integer=3)
    input_sample_channel, completed_sample_channel, progress_channel = channels
    #send sampling jobs to workers
    if isready(input_sample_channel) > 0
        @info "Sampling.."
        #WORKERS SAMPLE
        for worker in worker_pool[1:partitions]
            remote_do(get_sample_set, worker, input_sample_channel, completed_sample_channel, progress_channel)
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
