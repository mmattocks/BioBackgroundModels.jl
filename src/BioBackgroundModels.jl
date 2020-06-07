module BioBackgroundModels

import Base:copy,size
import Distances: euclidean
import Distributions:Distribution,VariateForm,Univariate,Dirichlet,Categorical,logpdf,isprobvec
import Distributed: RemoteChannel, myid, remote_do, rmprocs
import Printf: @sprintf
import Random: rand
import Serialization: serialize
import StatsFuns: logsumexp, logaddexp
import Statistics: mean
using BioSequences, BioSequences.FASTA, DataFrames, GenomicFeatures, ProgressMeter

include("HMM/HMM.jl")
export HMM
include("HMM/chain.jl")
export Chain_ID
include("API/genome_sampling.jl")
export setup_sample_jobs, execute_sample_jobs
include("API/EM_worker.jl")
include("API/EM_master.jl")
export setup_EM_jobs!, execute_EM_jobs!

include("EM/baum-welch.jl")
include("EM/churbanov.jl")
include("genome_sampling/partition_masker.jl")
include("genome_sampling/sequence_sampler.jl")
include("likelihood_funcs/bghmm_lh_matrix.jl")
export linear_likelihood, obs_set_likelihood
include("likelihood_funcs/hmm.jl")
include("utilities/observation_coding.jl")
include("utilities/BBG_analysis.jl")
include("utilities/BBG_progressmeter.jl")
include("utilities/utilities.jl")
export split_obs_sets
include("utilities/log_prob_sum.jl")
export lps

end # module
