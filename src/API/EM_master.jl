"""
    setup_EM_jobs!(job_ids, obs_sets; delta_thresh, chains, init_function)

Given a vector of 'Chain_IDs', and the appropriate obs_sets dictionary, return the appropriate tuple of channels and chains for the 'execute_EM_jobs!' function. Given an existing chains dict, resumes any existing non-converged chains (as assessed by the provided 'delta_thresh' value), otherwise initialises hmms for new chains with optionally user-specified 'init_function' (default 'autotransition_init').

"""
function setup_EM_jobs!(job_ids::Vector{Chain_ID}, obs_sets::Dict{String,Vector{LongSequence{DNAAlphabet{2}}}}; delta_thresh::Float64=1e-3,  chains::Dict{Chain_ID,Vector{EM_step}}=Dict{Chain_ID,Vector{EM_step}}(), init_function::Function=autotransition_init)
    #argument checking
    length(job_ids) < 1 && throw(ArgumentError("Empty job id vector!"))
    length(obs_sets) < 1 && throw(ArgumentError("Empty observation sets!"))

    no_input_hmms = length(job_ids)
    input_channel= RemoteChannel(()->Channel{Tuple}(no_input_hmms*3)) #channel to hold BHMM learning jobs
    output_channel= RemoteChannel(()->Channel{Tuple}(Inf)) #channel to take EM iterates off of

    code_dict = code_job_obs(job_ids, obs_sets)

    @showprogress 1 "Setting up HMMs..." for id in job_ids #for each jobid, add an initial BHMM to input_channel for EM_workers
        if haskey(chains, id) && length(chains[id]) > 0 #true if resuming from incomplete chain
            chain_end=chains[id][end]
            if !chain_end.converged || (chain_end.converged && chain_end.delta > delta_thresh)#push the last hmm iterate for nonconverged chains to the input channel with coded observations and values for chain resumption
                put!(input_channel, (id, chain_end.iterate, chain_end.hmm, chain_end.log_p, code_dict[(id.obs_id, id.order)]))
            else #skip any jobs that have converged from previous runs
                no_input_hmms -= 1
            end

        else 
            hmm = init_function(id.K, id.order, id.obs_id) #initialise first BHMM in chain
            chains[id] = Vector{EM_step}() #initialise the relevant chain
            obs=code_dict[(id.obs_id,id.order)]
            lh=obs_lh_given_hmm(obs,hmm,linear=false)
            put!(input_channel, (id, 1, hmm, lh, obs)) 
        end
    end
    return no_input_hmms, chains, input_channel, output_channel
end

"""
    execute_EM_jobs!(worker_pool, no_input_hmms, chains, input_channel, output_channel, chains_path; load_dict, EM_func, delta_thresh, max_iterates, verbose)

Given: a vector of worker IDs ([1] to use the master), the splatted output of 'setup_EM_jobs' (no_input_hmms, chains, input_channel, output_channel), and a valid path to save chains ('chains_path'); assign chains to workers using 'load_dict' to iteratively execute 'EM_func' (default linear_step) until 'delta_thresh' convergence criterion is obtained or 'max_iterates' have been calculated. Specifying 'verbose=true' gives some additional debug output from 'EM_converge!'.
"""
function execute_EM_jobs!(worker_pool::Vector{Int64}, no_input_hmms::Integer, chains::Dict{Chain_ID,Vector{EM_step}},  input_channel::RemoteChannel, output_channel::RemoteChannel, chains_path::String; load_dict=Dict{Int64,LoadConfig}(), EM_func::Function=linear_step, delta_thresh=1e-3, max_iterates=5000, verbose=false)
    #argument checking
    length(worker_pool) < 1 && throw(ArgumentError("Worker pool must contain one or more worker IDs!"))
    no_input_hmms < 1 && throw(ArgumentError("Zero input HMMs reported, likely continuing from chains already converged beyond default delta_thresh for setup_EM_jobs"))
    length(chains) < 1 && throw(ArgumentError("No chains supplied, likely job set from setup_EM_jobs passed incorrectly"))
    !isready(input_channel) && throw(ArgumentError("BHMM input channel has no contents, likely job set from setup_EM_jobs already executed"))

    #SEND BHMM FIT JOBS TO WORKERS
    if isready(input_channel) > 0
        @info "Fitting HMMs.."
        #WORKERS FIT HMMS
        for worker in worker_pool
            if worker in keys(load_dict)
                remote_do(EM_converge!, worker, input_channel, output_channel, no_input_hmms, load_config=load_dict[worker], EM_func=EM_func, delta_thresh=delta_thresh, max_iterates=max_iterates, verbose=verbose)
            else
                remote_do(EM_converge!, worker, input_channel, output_channel, no_input_hmms, EM_func=EM_func, delta_thresh=delta_thresh, max_iterates=max_iterates, verbose=verbose)
            end
        end
    else
        @warn "No input HMMs (all already converged?), skipping fitting.."
    end

    #GET LEARNT HMMS OFF REMOTECHANNEL, SERIALISE AT EVERY ITERATION, UPDATE PROGRESS METERS
    job_counter=no_input_hmms
    learning_meters=Dict{Chain_ID, ProgressHMM}()
    overall_meter=Progress(no_input_hmms,"Overall batch progress:")

    while job_counter > 0
        wait(output_channel)
        workerid, jobid, iterate, hmm, log_p, delta, converged, steptime = take!(output_channel)
        #either update an existing ProgressHMM meter or create a new one for the job
        if haskey(learning_meters, jobid) && iterate > 2
                update!(learning_meters[jobid], delta, steptime)
        else
            offset = workerid - 1
            if iterate <=2
                learning_meters[jobid] = ProgressHMM(delta_thresh, "$jobid on Wk $workerid:", offset, 2)
                update!(learning_meters[jobid], delta, steptime)
            else
                learning_meters[jobid] = ProgressHMM(delta_thresh, "$jobid on Wk $workerid:", offset, iterate)
                update!(learning_meters[jobid], delta, steptime)
            end
        end
        #push the hmm and related params to the results_dict
        push!(chains[jobid], EM_step(iterate, hmm, log_p, delta, converged))
        #try to serialize the results; catch interrupts and other errors to prevent corruption
        try
            serialize(chains_path, chains)
        catch e
            @warn "Serializing failed!"
            println(e)
        end
       
        #decrement the job counter, update overall progress meter, and save the current results dict on convergence or max iterate
        if converged || iterate == max_iterates
            job_counter -= 1
            ProgressMeter.update!(overall_meter, (no_input_hmms-job_counter))
            if !isready(input_channel) #if there are no more jobs to be learnt, retire the worker
                workerid!=1 && rmprocs(workerid)
            end
        end
    end

    #count converged & unconverged jobs, report results
    converged_counter = 0
    unconverged_counter = 0
    for (id, chain) in chains
        chain[end].converged == true ? (converged_counter += 1) : (unconverged_counter += 1)
    end

    @info "Background HMM batch learning task complete, $converged_counter converged jobs, $unconverged_counter jobs failed to converge in $max_iterates iterates since job start."
end
