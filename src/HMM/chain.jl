struct Chain_ID
    obs_id::String
    K::Int64
    order::Int64
    replicate::Int64
end

struct EM_step{F<:AbstractFloat,I<:Integer}
    iterate::I
    hmm::HMM
    log_p::F
    delta::F
    converged::Bool
end

