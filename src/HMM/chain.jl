struct Chain_ID{I<:Integer}
    obs_id::String
    K::I
    order::I
    replicate::I
    Chain_ID{I}(ob,K,or,r) where {I} = assert_chain_id(K,or,r) && new(ob,K,or,r)
end

function assert_chain_id(K, order, replicate)
    K < 1 && throw(ArgumentError("Chain_ID K (# of hmm states) must be a positive integer!"))
    order < 0 && throw(ArgumentError("Chain_ID symbol order must be zero or a positive integer!"))
    replicate < 1 && throw(ArgumentError("Chain_ID replicate must be a positive integer!"))
end

struct EM_step{F<:AbstractFloat,I<:Integer}
    iterate::I
    hmm::HMM
    log_p::F
    delta::F
    converged::Bool
    EM_step{F,I}(i,h,lp,d,c) where {F,I} = assert_step(i,h,lp) && new(i, h, lp, d, c)
end

function assert_step(iterate, hmm, log_p)
    iterate < 1 && throw(ArgumentError("EM_step iterate number must be a positive integer!"))
    assert_hmm(hmm)
    log_p > 0.0 && throw(ArgumentError("EM_step log probability value must be 0 or negative!"))
    return true
end