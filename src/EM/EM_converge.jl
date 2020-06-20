function EM_converge!(hmm_jobs::RemoteChannel, output_hmms::RemoteChannel, no_models::Integer; load_config::LoadConfig=LoadConfig(1:typemax(Int64)-1, 0:typemax(Int64)-1, [""]),  EM_func::Function=linear_step, delta_thresh=1e-3, max_iterates=5000, verbose=false)
    while isready(hmm_jobs)
        workerid = myid()
        jobid, start_iterate, hmm, job_norm, observations = load_balancer(no_models, hmm_jobs, load_config)
        jobid == 0 && break #no valid job for this worker according to load_table entry

        start_iterate > max_iterates - 1 && throw(ArgumentError("HMM chain $jobid is already longer ($start_iterate iterates) than specified max_iterates!"))
        curr_iterate = start_iterate

        #mask calculations here rather than ms_mle_step to prevent recalculation every iterate
        #build array of observation lengths

        obs_lengths = [findfirst(iszero,observations[o,:])-1 for o in 1:size(observations)[1]] #mask calculations here rather than mle_step to prevent recalculation every iterate
        EM_func==bw_step && (observations=transpose(observations))

        start_iterate == 1 && put!(output_hmms, (workerid, jobid, curr_iterate, hmm, 0.0, 0.0, false, 0.0)); #on the first iterate return the initial HMM for the chain right away
        verbose && @info "Fitting HMM on Wk $workerid, start iterate $start_iterate, $jobid with $(size(hmm)[1]) states and $(size(hmm)[2]) symbols..."

        curr_iterate += 1
        if curr_iterate == 2 #no delta value is available
            start=time()
            new_hmm, last_norm = EM_func(hmm, observations, obs_lengths)
            put!(output_hmms, (workerid, jobid, curr_iterate, new_hmm, last_norm, 0.0, false, time()-start))
        else #get the delta value from the channel-supplied job value to resume a chain properly
            start=time()
            new_hmm, last_norm = EM_func(hmm, observations, obs_lengths)
            delta = abs(lps(job_norm, -last_norm))
            put!(output_hmms, (workerid, jobid, curr_iterate, new_hmm, last_norm, delta, false, time()-start))
        end

        for i in curr_iterate:max_iterates
            start=time()
            new_hmm, norm = EM_func(new_hmm, observations, obs_lengths)
            curr_iterate += 1
            delta = abs(lps(norm, -last_norm))
            if delta < delta_thresh
                put!(output_hmms, (workerid, jobid, curr_iterate, new_hmm, norm, delta, true, time()-start))
                verbose && @info "$jobid converged after $(curr_iterate-1) EM steps"
                break
            else
                put!(output_hmms, (workerid, jobid, curr_iterate, new_hmm, norm, delta, false, time()-start))
                verbose && @info "$jobid EM step $(curr_iterate-1) delta $delta"
                last_norm = norm
            end
        end
    end
end