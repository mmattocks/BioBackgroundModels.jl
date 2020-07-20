struct Chain_ID
    obs_id::AbstractString
    K::Integer
    order::Integer
    replicate::Integer
    Chain_ID(obs_id,K,order,replicate) = assert_chain_id(K,order,replicate) && new(obs_id,K,order,replicate)
end

function assert_chain_id(K, order, replicate)
    K < 1 && throw(ArgumentError("Chain_ID K (# of hmm states) must be a positive integer!"))
    order < 0 && throw(ArgumentError("Chain_ID symbol order must be zero or a positive integer!"))
    replicate < 1 && throw(ArgumentError("Chain_ID replicate must be a positive integer!"))
    return true
end

struct EM_step
    iterate::Integer
    hmm::BHMM
    log_p::AbstractFloat
    delta::AbstractFloat
    converged::Bool
    EM_step(iterate,hmm,log_p,delta,converged)=assert_step(iterate,hmm,log_p) && new(iterate, hmm, log_p, delta, converged)
end

function assert_step(iterate, hmm, log_p)
    iterate < 1 && throw(ArgumentError("EM_step iterate number must be a positive integer!"))
    assert_hmm(hmm.a, hmm.A, hmm.B)
    log_p > 0.0 && throw(ArgumentError("EM_step log probability value must be 0 or negative!"))
    return true
end