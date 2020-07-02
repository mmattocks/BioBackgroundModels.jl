struct Report_Folder
    partition_id::String
    partition_report::Partition_Report
    replicate_report::Replicate_Report
    chain_reports::Dict{Chain_ID,Chain_Report}
end

function generate_reports(chains::Dict{Chain_ID,Vector{EM_step}}, test_sets::Dict{String,Vector{LongSequence{DNAAlphabet{2}}}})
    #check chains dict for problems

    length(chains)==0 && throw(ArgumentError("Empty chains dict!")) 
    any(chain->length(chain)<2, values(chains)) && throw(ArgumentError("Some chains are too short (<2 EM steps)! Probably not all chains have been operated on by EM workers yet. Try EM_converge first!"))
    unconv=sum(!chain[end].converged for chain in values(chains))
    unconv > 0 && @warn "Not all chains are converged to the selected step delta."
    length(test_sets)==0 && throw(ArgumentError("Empty test_sets dict!")) 

    chain_reports = report_chains(chains, test_sets)
    partition_reports = report_partitions(chain_reports)

    report_folders=Dict{String, Report_Folder}()

    for part_report in partition_reports
        rep_report=report_replicates(part_report.best_repset, chains)
        chain_subset=Dict{Chain_ID,Chain_Report}()
        for id in keys(chain_reports)
            id.obs_id==part_report.partition_id && (chain_subset[id]=chain_reports[id])
        end
        report_folders[part_report.partition_id]=Report_Folder(part_report.partition_id, part_report, rep_report, chain_subset)
    end

    return report_folders
end

