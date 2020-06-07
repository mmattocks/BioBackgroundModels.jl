function EM_converge!(hmm_jobs::RemoteChannel, output_hmms::RemoteChannel; EM_func::Function=linear_step, delta_thresh=1e-3, max_iterates=5000, verbose=false)
    while isready(hmm_jobs)
        workerid = myid()
        jobid::Tuple{String, Integer, Integer, Integer}, start_iterate::Integer, hmm::HMM, job_norm::AbstractFloat, observations::Matrix = take!(hmm_jobs)
        jobid == 0 && break #no valid job for this worker according to load_table entry

        start_iterate > max_iterates - 1 && throw(ArgumentError("HMM chain $jobid is already longer ($start_iterate iterates) than specified max_iterates!"))
        curr_iterate = start_iterate

        #mask calculations here rather than ms_mle_step to prevent recalculation every iterate
        #build array of observation lengths
        obs_lengths = [findfirst(iszero,observations[o,:])-1 for o in 1:size(observations)[1]] #mask calculations here rather than mle_step to prevent recalculation every iterate

        start_iterate == 1 && put!(output_hmms, (workerid, jobid, curr_iterate, hmm, 0, 0, false)); #on the first iterate return the initial HMM for the chain right away
        verbose && @info "Fitting HMM, start iterate $start_iterate, job ID $jobid with $(size(hmm.π)[1]) states and $(length(hmm.D[1].support)) symbols..."

        curr_iterate += 1
        if curr_iterate == 2 #no delta value is available
            new_hmm, last_norm = EM_func(hmm, observations, obs_lengths)
            put!(output_hmms, (workerid, jobid, curr_iterate, new_hmm, last_norm, 0, false))
        else #get the delta value from the channel-supplied job value to resume a chain properly
            new_hmm, last_norm = EM_func(hmm, observations, obs_lengths)
            delta = abs(lps(job_norm, -last_norm))
            put!(output_hmms, (workerid, jobid, curr_iterate, new_hmm, last_norm, delta, false))
        end

        for i in curr_iterate:max_iterates
            new_hmm, norm = EM_func(new_hmm, observations, obs_lengths)
            curr_iterate += 1
            delta = abs(lps(norm, -last_norm))
            if delta < delta_thresh
                put!(output_hmms, (workerid, jobid, curr_iterate, new_hmm, norm, delta, true))
                verbose && @info "Job ID $jobid converged after $(curr_iterate-1) EM steps"
                break
            else
                put!(output_hmms, (workerid, jobid, curr_iterate, new_hmm, norm, delta, false))
                verbose && @info "Job ID $jobid EM step $(curr_iterate-1) delta $delta"
                last_norm = norm
            end
        end
    end
end