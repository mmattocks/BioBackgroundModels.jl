struct Replicate_Report
    ids::Vector{Chain_ID} #replicates 
    state_vecs::Dict{Chain_ID,Matrix{Float64}} #vectors of autotransition probabilities evolving by iterate, concatanated to matrices, indexed by Chain_ID
    emission_arrays::Dict{Chain_ID,Array{Float64}} #iterate*symbol*state_vecs
    sorted_id1_states::Vector{Integer}
    sort_dicts::Dict{Chain_ID,Dict{Integer,Integer}}#emission channels sorted on the basis of euclidean closeness to ids[1]
    sorted_symbols::Dict{Integer,Vector{Integer}}


    Replicate_Report(ids,state_vecs,emission_arrays, sorted_id1_states,sort_dicts,sorted_symbols) = assert_repreport(ids) && new(ids,state_vecs,emission_arrays,sorted_id1_states,sort_dicts,sorted_symbols)
end

function assert_repreport(ids::Vector{Chain_ID})
    Ks=Vector{Int64}()
    orders=Vector{Int64}()
    replicates=Vector{Int64}()
    obsids=Vector{String}()

    for id in ids
        push!(Ks,id.K)
        push!(orders,id.order)
        push!(obsids,id.obs_id)
        push!(replicates, id.replicate)
    end

    length(unique(Ks)) > 1 && throw(ArgumentError("Replicate set includes chains with different K state numbers! All chains must have HMMs with the same number of states."))
    length(unique(orders)) > 1 && throw(ArgumentError("Replicate set includes chains with different order coding numbers! All chains must have HMMs with the same order coding."))
    !allunique(replicates) && throw(ArgumentError("Replicate set includes chains with the same replicate #. Replicates should be unique!"))

    return true
end

function report_replicates(repset::Vector{Chain_ID},chains::Dict{Chain_ID,Vector{EM_step}})
    assert_repreport(repset) #check that the repset will pass its constructor before doing the work
    emission_arrays=Dict{Chain_ID,AbstractArray{AbstractFloat}}()
    state_arrays=Dict{Chain_ID,AbstractArray{AbstractFloat}}()
    
    for id in repset
        emission_array=zeros(length(chains[id]),4^(id.order+1),id.K)
        state_array=zeros(length(chains[id]),id.K)
        for (it, step) in enumerate(chains[id])
            for k in 1:id.K
                emission_array[it,:,k]=step.hmm.B[k].p
                state_array[it,k]=step.hmm.A[k,k]
            end
        end
        emission_arrays[id]=emission_array
        state_arrays[id]=state_array
    end

    println(repset)
    println(emission_arrays)
    println(state_arrays)
    sort_dicts,sorted_id1_states,sorted_symbols=sort_emitters_by_distance!(repset,emission_arrays)

    return Replicate_Report(repset,state_arrays,emission_arrays,sorted_id1_states, sort_dicts,sorted_symbols)
end

function sort_emitters_by_distance!(repset,emission_arrays)
    id1=repset[1]
    length(repset)==1 && (return Dict{Chain_ID,Dict{Integer,Integer}}(),[1:id1.K...],[1:4^id1.order+1])

    final_emissions = zeros(length(repset),size(emission_arrays[id1])[2:3]...)
    for (n,id) in enumerate(repset)
        final_emissions[n,:,:] = emission_arrays[id][end,:,:]
    end

    distances=zeros(length(repset)-1,id1.K,id1.K)
    for rep in 1:length(repset)-1
        for k in 1:id1.K, j in 1:id1.K
            distances[rep,k,j]=euclidean(final_emissions[1,:,k],final_emissions[1+rep,:,j])
        end
    end

    sortdicts=Dict{Chain_ID,Dict{Integer,Integer}}()
    sorted_id1_states=Vector{Integer}()
    sorted_symbols=Dict{Integer,Vector{Integer}}()
    for rep in 1:length(repset)-1
        sortdicts[repset[rep+1]]=Dict{Integer,Integer}()
    end

    for i in 1:id1.K
        id1_mindist_K=findmin([sum([minimum(distances[rep,k,:]) for rep in 1:length(repset)-1]) for k in 1:id1.K])[2]#find the state from the first replicate that has minimum cumulative distance to the closest state in the other replicates, excluding already-chosen states by spiking their values
        push!(sorted_id1_states,id1_mindist_K)
        for rep in 1:length(repset)-1
            rep_k=findmin(distances[rep,id1_mindist_K,:])[2]
            sortdicts[repset[rep+1]][id1_mindist_K]=rep_k
            distances[rep,:,rep_k].=1.0 #mask
            distances[1,id1_mindist_K,:].=1.0 #mask
        end

        symbol_distances=sum(hcat([euclidean.(final_emissions[1,:,id1_mindist_K],final_emissions[n+1,:,sortdicts[id][id1_mindist_K]]) for (n,id) in enumerate(repset[2:end])]...),dims=2)
        sorted_symbols[id1_mindist_K]=sortperm(symbol_distances[:,1])
    end

    return sortdicts,sorted_id1_states,sorted_symbols
end


function Base.show(io::IO, report::Replicate_Report; top_states::Integer=report.ids[1].K, top_symbols::Integer=2)
    printstyled(io, "REPLICATE CONVERGENCE REPORT\n", color=:green)
    println(" ̶̶ ̶ ̶ ̶ ̶̶ ̶ ̶ ̶ ̶̶ ̶ ̶ ̶ ̶̶ ̶ ̶ ̶ ̶̶ ̶ ̶ ̶ ̶̶ ̶ ̶ ̶ ̶̶ ̶ ̶ ̶ ̶̶ ̶ ̶ ̶ ̶̶ ̶ ̶ ̶ ̶̶ ̶ ̶ ̶ ̶̶ ̶ ̶ ̶ ̶̶ ̶ ̶ ̶ ̶̶ ̶ ̶ ̶ ̶̶ ̶ ̶ ̶ ̶̶ ̶ ̶ ̶")
    for i in 1:top_states
        id=report.ids[1]
        id1_state=report.sorted_id1_states[i]
        itmax=maximum([size(report.state_vecs[id],1) for id in report.ids])

        autoplt=lineplot(report.state_vecs[id][:,id1_state],
                        title="Autotransition prob. convergence",
                        name="$(report.ids[1]) K$id1_state",
                        xlim=(0,ceil(itmax, sigdigits=2)),
                        ylim=(0,1),
                        xlabel="iterate",
                        ylabel="p")

        symplot=lineplot(report.emission_arrays[id][:,report.sorted_symbols[id1_state][1],id1_state],
                        report.emission_arrays[id][:,report.sorted_symbols[id1_state][2],id1_state], 
                        title="Symbol prob. convergence",
                        name="$(report.ids[1]) K$id1_state",
                        xlim=(0,1),
                        ylim=(0,1),
                        xlabel="Symbol 1 p",
                        ylabel="S2 p")

        for (n,id) in enumerate(report.ids)
            if n > 1
                lineplot!(autoplt, report.state_vecs[id][:,report.sort_dicts[id][id1_state]],
                            name="$(report.ids[n]) K$(report.sort_dicts[id][id1_state])")
                lineplot!(symplot,
                            report.emission_arrays[id][:,report.sorted_symbols[id1_state][1],report.sort_dicts[id][id1_state]],
                            report.emission_arrays[id][:,report.sorted_symbols[id1_state][2],report.sort_dicts[id][id1_state]],
                            name="$(report.ids[n]) K$(report.sort_dicts[id][id1_state])")
            end
        end

        show(autoplt)
        println("\n")
        show(symplot)
        println("\n")
    end
end