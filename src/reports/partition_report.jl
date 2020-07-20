struct Order_Report
    converged_K::Vector{Float64}
    converged_lh::Vector{Float64}
    failed_K::Vector{Float64}
    failed_lh::Vector{Float64}
end

struct Partition_Report
    partition_id::String
    naive_lh::Float64
    orddict::Dict{Integer,Order_Report}
    best_model::Tuple{Chain_ID,BHMM}
    best_repset::Vector{Chain_ID}
end

function Base.show(io::IO, report::Partition_Report)
    printstyled(io, "BACKGROUND HMM TRAINING REPORT\n", bold=true)
    printstyled(io, "Genome partition id: $(report.partition_id)\n", bold=true, color=:magenta)

    for (order, ordreport) in report.orddict
        if length(ordreport.failed_lh)>0
            ordplot=scatterplot(ordreport.converged_K, ordreport.converged_lh; name="Converged", title="Order $order HMMs", xlabel="K# States", ylabel="P(O|θ)", color=:green, xlim=(1,maximum(ordreport.converged_K)), ylim=(floor(min(minimum(ordreport.converged_lh),minimum(ordreport.failed_lh),report.naive_lh),sigdigits=2),ceil(max(maximum(ordreport.converged_lh),maximum(ordreport.failed_lh),report.naive_lh),sigdigits=2)))
            length(ordreport.failed_K)>=1 && scatterplot!(ordplot, ordreport.failed_K, ordreport.failed_lh; name="Unconverged",color=:red)
            lineplot!(ordplot, [report.naive_lh for i in 1:maximum(ordreport.converged_K)], name="Naive",color=:magenta)
            show(ordplot)
            println()
        else
            ordplot=scatterplot(ordreport.converged_K, ordreport.converged_lh; name="Converged", title="Order $order HMMs", xlabel="K# States", ylabel="P(O|θ)", color=:green, xlim=(1,maximum(ordreport.converged_K)), ylim=(floor(min(minimum(ordreport.converged_lh),report.naive_lh),sigdigits=2),ceil(max(maximum(ordreport.converged_lh),report.naive_lh),sigdigits=2)))
            lineplot!(ordplot, [report.naive_lh for i in 1:maximum(ordreport.converged_K)], name="Naive",color=:magenta)
            show(ordplot)
            println()
        end
    end
end

function report_partitions(chain_reports::Dict{Chain_ID,Chain_Report})
    partitions=Vector{String}()
    orders=Vector{Integer}()
    Ks=Vector{Integer}()
    reps=Vector{Integer}()
    for id in keys(chain_reports)
        !in(id.obs_id, partitions) && push!(partitions, id.obs_id)
        !in(id.order, orders) && push!(orders, id.order)
        !in(id.K, Ks) && push!(Ks, id.K)
        !in(id.replicate, reps) && push!(reps, id.replicate)
    end

    reports=Vector{Partition_Report}()
    for partition in partitions
        best_lh,best_rep=-Inf,Chain_ID("empty",1,0,1)
        naive_lh=1.
        orddict=Dict{Integer,Order_Report}()
        for order in orders
            conv_statevec=Vector{Float64}()
            fail_statevec=Vector{Float64}()
            conv_lh_vec=Vector{Float64}()
            fail_lh_vec=Vector{Float64}()
            for (id,report) in chain_reports
                if id.obs_id==partition && id.order==order
                    report.test_lh > best_lh && (best_lh=report.test_lh;best_rep=id)
                    naive_lh==1. && (naive_lh=chain_reports[id].naive_lh)
                    chain_reports[id].converged ? (push!(conv_statevec,id.K);
                    push!(conv_lh_vec,chain_reports[id].test_lh)) : (push!(fail_statevec,id.K);
                    push!(fail_lh_vec,chain_reports[id].test_lh))
                end
            end
            orddict[order]=Order_Report(conv_statevec,conv_lh_vec,fail_statevec,fail_lh_vec)
        end
        best_model=(best_rep,chain_reports[best_rep].final_hmm)
        best_repset=[Chain_ID(partition, best_rep.K, best_rep.order, rep) for rep in reps]
        push!(reports,Partition_Report(partition,naive_lh,orddict,best_model,best_repset))
    end
    return reports
end

