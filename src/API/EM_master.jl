#function to setup an HMM chains dictionary and RemoteChannel for learning jobs, given a vector of state #s, order_nos, replicates to train, the dictionary to fill, the RemoteChannel and the training sequences
#resumes any existing non-converged chains, otherwise initialises hmms for new chains given provided constants
function setup_EM_jobs!(job_ids::Vector{Chain_ID}, obs_sets::Dict{String,Vector{LongSequence{DNAAlphabet{2}}}};  chains::Dict{Chain_ID,Vector{EM_step}}=Dict{Chain_ID,Vector{EM_step}}(), init_function::Function=autotransition_init)
    no_input_hmms = length(job_ids)
    input_channel= RemoteChannel(()->Channel{Tuple}(no_input_hmms*3)) #channel to hold HMM learning jobs
    output_channel= RemoteChannel(()->Channel{Tuple}(Inf)) #channel to take EM iterates off of

    code_dict = code_job_obs(job_ids, obs_sets)

    @showprogress 1 "Setting up HMMs..." for id in job_ids #for each jobid, add an initial HMM to input_channel for EM_workers
        if haskey(chains, id) && length(chains[id]) > 0 #true if resuming from incomplete chain
            chain_end=chains[id][end]
            if !chain_end.converged #push the last hmm iterate for nonconverged chains to the input channel with coded observations and values for chain resumption
                put!(input_channel, (id, chain_end.iterate, chain_end.hmm, chain_end.log_p, code_dict[(id.obs_id, id.order)]))
            else #skip any jobs that have converged from previous runs
                no_input_hmms -= 1
            end

        else 
            hmm = init_function(id.K, id.order) #initialise first HMM in chain
            chains[id] = Vector{EM_step}() #initialise the relevant chain
            put!(input_channel, (id, 1, hmm, 0.0, code_dict[(id.obs_id, id.order)])) 
        end
    end

    return no_input_hmms, chains, input_channel, output_channel
end

function execute_EM_jobs!(worker_pool::Vector{Int64}, no_input_hmms::Integer, chains::Dict{Chain_ID,Vector{EM_step}},  input_channel::RemoteChannel, output_channel::RemoteChannel, chains_path::String,; EM_func::Function=linear_step, delta_thresh=1e-3, max_iterates=5000, verbose=false)
    #SEND HMM FIT JOBS TO WORKERS
    if isready(input_channel) > 0
        @info "Fitting HMMs.."
        #WORKERS FIT HMMS
        for worker in worker_pool
            remote_do(EM_converge!, worker, input_channel, output_channel, EM_func=EM_func, delta_thresh=delta_thresh, max_iterates=max_iterates, verbose=verbose)
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
        workerid, jobid, iterate, hmm, log_p, delta, converged = take!(output_channel)
        #either update an existing ProgressHMM meter or create a new one for the job
        if haskey(learning_meters, jobid) && iterate > 2
                update!(learning_meters[jobid], delta)
        else
            offset = workerid - 1
            if iterate <=2
                learning_meters[jobid] = ProgressHMM(delta_thresh, "Jobid $jobid on worker $workerid converging:", offset, 2)
            else
                learning_meters[jobid] = ProgressHMM(delta_thresh, "Jobid $jobid on worker $workerid converging:", offset, iterate)
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
                rmprocs(workerid)
            end
        end
    end

    #count converged & unconverged jobs, report results
    converged_counter = 0
    unconverged_counter = 0
    for (id, chain) in chains
        chain.converged == true ? (converged_counter += 1) : (unconverged_counter += 1)
    end

    @info "Background HMM batch learning task complete, $converged_counter converged jobs, $unconverged_counter jobs failed to converge in $max_iterates iterates since job start."
end
