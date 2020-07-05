module BioBackgroundModels

import Base:copy,size
import Distances: euclidean
import Distributions:Distribution,VariateForm,Univariate,Dirichlet,Categorical,logpdf,isprobvec
import Distributed: RemoteChannel, myid, remote_do, rmprocs
import MCMCChains: Chains, ChainDataFrame, heideldiag
import Printf: @sprintf
import Random: rand
import Serialization: serialize
import StatsFuns: logsumexp, logaddexp
import Statistics: mean
import UnicodePlots: lineplot,lineplot!, scatterplot,scatterplot!
using BioSequences, DataFrames, FASTX, GFF3, ProgressMeter

include("HMM/HMM.jl")
export HMM
include("HMM/chain.jl")
export Chain_ID, EM_step
include("API/genome_sampling.jl")
export setup_sample_jobs, execute_sample_jobs
include("API/EM_master.jl")
export setup_EM_jobs!, execute_EM_jobs!

include("EM/baum-welch.jl")
include("EM/churbanov.jl")
include("utilities/load_balancer.jl")
export LoadConfig
include("EM/EM_converge.jl")
include("genome_sampling/partition_masker.jl")
include("genome_sampling/sequence_sampler.jl")
include("likelihood_funcs/bg_lh_matrix.jl")
include("likelihood_funcs/hmm.jl")
export obs_lh_given_hmm
include("reports/chain_report.jl")
export generate_reports
include("utilities/observation_coding.jl")
include("utilities/BBG_analysis.jl")
include("utilities/BBG_progressmeter.jl")
include("utilities/HMM_init.jl")
include("utilities/model_display.jl")
include("utilities/utilities.jl")
export split_obs_sets
include("utilities/log_prob_sum.jl")
export lps
include("reports/partition_report.jl")
include("reports/replicate_convergence.jl")
include("dev.jl")

end # module
