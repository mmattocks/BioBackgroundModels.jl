struct Chain_ID{I<:Integer}
    obs_id::String
    K::I
    order::I
    replicate::I
end

struct EM_step{F<:AbstractFloat,I<:Integer}
    iterate::I
    hmm::HMM
    log_p::F
    delta::F
    converged::Bool
end

