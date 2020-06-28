struct LoadConfig
    k_range::UnitRange
    o_range::UnitRange
    blacklist::Vector{Chain_ID}
    whitelist::Vector
    function LoadConfig(k_range, o_range; blacklist=Vector{Chain_ID}(), whitelist=Vector{Chain_ID}())
        assert_loadconfig(k_range, o_range) && new(k_range, o_range, blacklist, whitelist)
    end
end

function assert_loadconfig(k_range,o_range)
    k_range[1] < 1 && throw(ArgumentError("Minimum value of LoadConfig k_range is 1!"))
    o_range[1] < 0 && throw(ArgumentError("Minimum value of LoadConfig o_range is 0!"))
    return true
end

#subfunc to handle balancing memory load on dissimilar machines in cluster
function load_balancer(no_models::Integer, hmm_jobs::RemoteChannel, config::LoadConfig)
    lb_job_counter = 1

    jobid::Chain_ID, start_iterate::Integer, hmm::HMM, job_norm::AbstractFloat, observations::Matrix = take!(hmm_jobs)

    while (jobid.K < config.k_range[1] || jobid.K > config.k_range[end] || jobid.order < config.o_range[1] || jobid.order > config.o_range[end] || jobid in config.blacklist || (length(config.whitelist)>0 && !(jobid in config.whitelist))) && lb_job_counter <= no_models #while a job prohibited by load table, keep putting the job back and drawing a new one
        put!(hmm_jobs, (jobid, start_iterate, hmm, job_norm, observations))
        jobid, start_iterate, hmm, job_norm, observations = take!(hmm_jobs)
        lb_job_counter += 1
    end
    lb_job_counter > no_models ? (return 0,0,0,0,0) : (return jobid, start_iterate, hmm, job_norm, observations)
end