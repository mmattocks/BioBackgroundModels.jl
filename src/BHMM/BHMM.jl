"""
    BHMM([a::AbstractVector{T}, ]A::AbstractMatrix{T}, B::AbstractVector{<:Categorical}) where F where T
Build an BHMM with transition matrix `A` and observations distributions `B`.  
If the initial state distribution `a` is not specified, a uniform distribution is assumed. 
Observations distributions can be of different types (for example `Normal` and `Exponential`).  
However they must be of the same dimension (all scalars or all multivariates).
# Example
```julia
hmm = BHMM([0.9 0.1; 0.1 0.9], [Normal(0,1), Normal(10,1)])
```
"""
struct BHMM{T} <: AbstractHMM{Univariate}
    a::AbstractVector{T}
    A::AbstractMatrix{T}
    B::AbstractVector{<:Categorical}
    partition::String
    BHMM{T}(a, A, B, partition="") where {T} = assert_BHMM(a, A, B) && new(a, A, B, partition) 
end

BHMM(a::AbstractVector{T}, A::AbstractMatrix{T}, B::AbstractVector{<:Categorical}, partition ...) where {T} = BHMM{T}(a, A, B, partition ...)
BHMM(A::AbstractMatrix{T}, B::AbstractVector{<:Categorical}, partition...) where {T} = BHMM{T}(ones(size(A)[1])/size(A)[1], A, B, partition...)

copy(hmm::BHMM) = BHMM(copy(hmm.a), copy(hmm.A), copy(hmm.B), hmm.partition)
size(hmm::BHMM) = (length(hmm.B), length(hmm.B[1].p))


"""
    assert_BHMM(a, A, B)
Throw an `ArgumentError` if the initial state distribution and the transition matrix rows does not sum to 1,
and if the observations distributions does not have the same dimensions.
"""
function assert_BHMM(a::AbstractVector, 
                    A::AbstractMatrix, 
                    B::AbstractVector{<:Categorical})
    !isprobvec(a) && throw(ArgumentError("Initial state vector a is not a valid probability vector!")) 
    !istransmat(A) && throw(ArgumentError("Transition matrix A is not valid!")) 

    !all([length(d.p) for d in B] .== length(B[1].p)) && throw(ArgumentError("All distributions must have the same dimensions"))
    !(length(a) == size(A,1) == length(B)) && throw(ArgumentError("Length of initial state vector a, dimension of transition matrix A, and number of distributions B are not the same!"))
    return true
end

function Base.show(io::IO, hmm::BHMM)
    println(io, "Background BHMM")
    println(io, "State Initial and Transition Probabilities")
    print(io, "a: ")
    show(io, hmm.a)
    println(io)
    print(io, "A: ")
    display(hmm.A)
    println(io)
    println(io, "INFORMATIVE SYMBOLS BY STATE")
    for (n,d) in enumerate(hmm.B)
        print(io, "K$n ")
        print_emitters(d)
    end
end