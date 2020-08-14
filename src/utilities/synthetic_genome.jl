function synthetic_sample(no_obs::Integer, obsl, bhmm_vec::AbstractVector{<:BHMM}, bhmm_dist::Categorical, spikes::AbstractVector{<:AbstractMatrix{<:AbstractFloat}}, spike_instruct::AbstractVector{<:Tuple{<:Bool, <:Tuple}})
    !(typeof(obsl) <:UnitRange || typeof(obsl) <:Integer) && throw(ArgumentError("obsl must be Integer or UnitRange"))
    length(spike_instruct)!=length(spikes) && throw(ArgumentError("spike_instruct must be as long as spikes"))

    obs, hmm_truth=obs_array_from_bhmms(no_obs,obsl,bhmm_vec,bhmm_dist)

    spike_truth = spike_obs!(obs, spikes, spike_instruct)

    bg_scores = score_synthetic(obs, bhmm_vec, hmm_truth)

    return obs, bg_scores, hmm_truth, spike_truth
end

function obs_array_from_bhmms(no_obs, obsl, bhmm_vec, bhmm_dist)
    obs=zeros(UInt8,max(obsl)+1, no_obs)
    truthvec=zeros(UInt8,no_obs)
    for o in 1:no_obs
        hmm_idx=rand(bhmm_dist)
        truthvec[o]=hmm_idx
        bhmm=bhmm_vec[hmm_idx]
        typeof(obsl)<:UnitRange ? (l=rand(obsl)) : (l=obsl)
        obs[1:l,o]=rand(bhmm, l)
    end
    return obs, truthvec
end
        
function spike_obs!(obs, spikes, spike_instruct)
    truthmat=falses(no_obs,length(spikes))
    for (s,spike) in enumerate(spikes)
        structural,ins=spike_instruct[s]
        structural ? (truthmat[:,s] = spike_struc!(obs, spike, ins...)) : (truthmat[:,s] = spike_irreg!(obs, o, spike, ins...))
    end
end

function spike_struc!(obs, spike, frac_obs, periodicity)
    truth=falses(size(obs,2))
    for o in 1:size(obs,2)
        if rand() < frac_obs
            truth[o]=true
            rand()<.5 ? (source=revcomp_pwm(spike)) : (source=spike)
            pos=rand(1:periodicity)
            oidx=findfirst(iszero,obs[:,o])
            while pos<oidx
                pos_ctr=pos
                pwm_ctr=1
                while pos_ctr<oidx&&pwm_ctr<=size(source,1)
                    obs[pos_ctr,o]=rand(Categorical(source[pwm_ctr,:]))
                    pos_ctr+=1
                    pwm_ctr+=1
                end
                pos+=periodicity+pwm_ctr
            end
        end
    end
    return truth
end

function spike_irreg!(obs, spike, frac_obs, recur)
    truth=falses(size(obs,2))
    for o in 1:size(obs,2)
        oidx=findfirst(iszero,obs[:,o])
        if rand()<frac_obs
            truth[o]=true
            for r in 1:rand(recur)
                rand()<.5 ? (source=nnlearn.revcomp_pwm(spike)) : source=spike
                pos=rand(1:oidx-1)
                pwm_ctr=1
                while pos<oidx && pwm_ctr<=size(source,1)
                    obs[pos,o]=rand(Categorical(source[pwm_ctr,:]))
                    pos+=1
                    pwm_ctr+=1
                end
            end
        end
    end
    return truth
end

function score_synthetic(obs, bhmm_vec, hmm_truth)
    lh_mat=zeros(size(obs,1)-1, size(obs,2))
    for o in 1:size(obs,2)
        oidx=findfirst(iszero,obs[:,o])
        lh_mat[1:oidx-1,o]=BGHMM_symbol_lh(obs[1:oidx,o:o],bhmm_vec[hmm_truth[o]])
    end
    return lh_mat
end