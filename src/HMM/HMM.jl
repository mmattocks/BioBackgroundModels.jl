#AbstractHMM and HMM forked from HMMBase to maintain compatibility with MS_HMMBase

"""
    AbstractHMM{F<:VariateForm}
An HMM type must at-least implement the following interface:
```julia
struct CustomHMM{F,T} <: AbstractHMM{F}
    π0::AbstractVector{T}              # Initial state distribution
    π::AbstractMatrix{T}               # Transition matrix
    D::AbstractVector{Distribution{F}} # Observations distributions
    # Custom fields ....
end
```
"""
abstract type AbstractHMM{F<:VariateForm} end

"""
    HMM([π0::AbstractVector{T}, ]π::AbstractMatrix{T}, D::AbstractVector{<:Distribution{F}}) where F where T
Build an HMM with transition matrix `π` and observations distributions `D`.  
If the initial state distribution `π0` is not specified, a uniform distribution is assumed. 
Observations distributions can be of different types (for example `Normal` and `Exponential`).  
However they must be of the same dimension (all scalars or all multivariates).
# Example
```julia
hmm = HMM([0.9 0.1; 0.1 0.9], [Normal(0,1), Normal(10,1)])
```
"""
struct HMM{F,T} <: AbstractHMM{F}
    π0::Vector{T}
    π::Matrix{T}
    D::Vector{Distribution{F}}
    HMM{F,T}(π0, π, D) where {F,T} = assert_hmm(π0, π, D) && new(π0, π, D) 
end

HMM(π0::AbstractVector{T}, π::AbstractMatrix{T}, D::AbstractVector{<:Distribution{F}}) where {F,T} = HMM{F,T}(π0, π, D)
HMM(π::AbstractMatrix{T}, D::AbstractVector{<:Distribution{F}}) where {F,T} = HMM{F,T}(ones(size(π)[1])/size(π)[1], π, D)

"""
    assert_hmm(π0, π, D)
Throw an `ArgumentError` if the initial state distribution and the transition matrix rows does not sum to 1,
and if the observations distributions does not have the same dimensions.
"""
function assert_hmm(π0::AbstractVector, 
                    π::AbstractMatrix, 
                    D::AbstractVector{<:Distribution})
    !isprobvec(π0) && throw(ArgumentError("Initial state vector π0 is not a valid probability vector!")) 
    !istransmat(π) && throw(ArgumentError("Transition matrix π is not valid!")) 

    !all([length(d.p) for d in D] .== length(D[1].p)) && throw(ArgumentError("All distributions must have the same dimensions"))
    !(length(π0) == size(π,1) == length(D)) && throw(ArgumentError("Length of initial state vector π0, dimension of transition matrix π, and number of distributions D are not the same!"))
    return true
end

issquare(A::AbstractMatrix) = size(A,1) == size(A,2)
istransmat(A::AbstractMatrix) = issquare(A) && all([isprobvec(A[i,:]) for i in 1:size(A,1)])

size(hmm::AbstractHMM, dim = :) = (length(hmm.D), length(hmm.D[1].p))[dim]

copy(hmm::HMM) = HMM(copy(hmm.π0), copy(hmm.π), copy(hmm.D))