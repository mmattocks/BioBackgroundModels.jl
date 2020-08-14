module BioBackgroundModels

import Base:copy,size
import Distances: euclidean
import Distributions:Univariate,Dirichlet,Categorical,logpdf,isprobvec
import Distributed: RemoteChannel, myid, remote_do, rmprocs
import HMMBase: AbstractHMM, assert_hmm, istransmat
import MCMCChains: Chains, ChainDataFrame, heideldiag
import Printf: @sprintf
import Random: rand, shuffle
import Serialization: serialize
import StatsFuns: logsumexp, logaddexp
import Statistics: mean
import UnicodePlots: lineplot,lineplot!, scatterplot,scatterplot!
using BioSequences, DataFrames, FASTX, GFF3, ProgressMeter

include("BHMM/BHMM.jl")
export BHMM
include("EM/chain.jl")
export Chain_ID, EM_step
include("API/genome_sampling.jl")
export setup_sample_jobs, execute_sample_jobs
include("API/EM_master.jl")
export setup_EM_jobs!, execute_EM_jobs!

include("reports/chain_report.jl")
include("reports/partition_report.jl")
include("reports/replicate_convergence.jl")
include("API/reports.jl")
export generate_reports

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
include("utilities/observation_coding.jl")
include("utilities/BBG_analysis.jl")
include("utilities/BBG_progressmeter.jl")
include("utilities/HMM_init.jl")
include("utilities/model_display.jl")
include("utilities/utilities.jl")
export split_obs_sets
include("utilities/log_prob_sum.jl")
export lps
end # module
