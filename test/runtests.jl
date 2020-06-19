using BioBackgroundModels,BioSequences,DataFrames,Distributed,Distributions,FASTX,ProgressMeter,Random,StatsFuns,Test
import BioBackgroundModels: autotransition_init, load_balancer, CompoundAlphabet, get_order_n_seqs, code_seqs, code_job_obs, make_padded_df, add_partition_masks!, partition_genome_coordinates, setup_sample_jobs, meta_to_feature_coord, get_strand_dict, build_scaffold_seq_dict, get_feature_params_from_metacoord, determine_sample_window, fetch_sequence, get_sample_set, assert_hmm, linear_step, get_BGHMM_symbol_lh,make_padded_df,add_partition_masks!,BGHMM_likelihood_calc,linear_likelihood

include("synthetic_sequence_gen.jl")
include("ref_fns.jl")

#JOB FILEPATHS
#GFF3 feature database, FASTA genome and index paths
Sys.islinux() ? genome =  (@__DIR__) * "/synthetic.fna" : genome = (@__DIR__) * "\\synthetic.fna"
Sys.islinux() ? index =  (@__DIR__) * "/synthetic.fna.fai" : index = (@__DIR__) * "\\synthetic.fna.fai"
Sys.islinux() ? gff =  (@__DIR__) * "/synthetic.gff3" : gff = (@__DIR__) * "\\synthetic.gff3"
Sys.islinux() ? posfasta =  (@__DIR__) * "/syntheticpos.fa" : posfasta = (@__DIR__) * "\\syntheticpos.fa"

!isfile(genome) && print_synthetic_fasta(genome)
!isfile(index) && print_synthetic_index(index)
!isfile(gff) && print_synthetic_gff(gff)
!isfile(posfasta) && print_synthetic_position(posfasta)

sample_record_dfs = Dict{String,DataFrame}()

Random.seed!(1)

@testset "HMM/HMM.jl constructors" begin
    good_π0 = [.20,.30,.20,.30]
    bad_probvec_π0 = [.25,.25,.25,.27]
    good_π = fill(.25,4,4)
    bad_probvec_π = deepcopy(good_π)
    bad_probvec_π[:,3].=.27
    bad_size_π = fill(.25,5,4)
    good_D = Categorical([.25,.25,.25,.25])
    bad_size_D = Categorical([.2,.2,.2,.2,.2])
    good_Dvec = [good_D for i in 1:4]
    bad_size_Dvec =[good_D, good_D, good_D, bad_size_D]
    bad_length_Dvec=[good_D for i in 1:5]

    @test_throws ArgumentError HMM(bad_probvec_π0, good_π, good_Dvec)
    @test_throws ArgumentError HMM(good_π0, bad_probvec_π, good_Dvec)
    @test_throws ArgumentError HMM(good_π0, bad_size_π, good_Dvec)
    @test_throws ArgumentError HMM(good_π0, good_π, bad_size_Dvec)
    @test_throws ArgumentError HMM(good_π0, good_π, bad_length_Dvec)

    hmm=HMM(good_π0, good_π, good_Dvec)
    @test typeof(hmm)==HMM{Univariate, Float64}
    @test hmm.π0==good_π0
    @test hmm.π==good_π
    @test hmm.D==good_Dvec 
    @test size(hmm) == (4,4)

    hmm_copy=copy(hmm)
    @test hmm.π0==hmm_copy.π0
    @test hmm.π==hmm_copy.π
    @test hmm.D==hmm_copy.D

    hmm2=HMM(good_π, good_Dvec)
    @test hmm2.π0==[.25,.25,.25,.25]
    @test hmm2.π==good_π
    @test hmm2.D==good_Dvec
end

@testset "HMM/chain.jl Chain_ID and EM_step constructors" begin
    obs_id="test"
    K,zero_k,neg_K=4,0,-1
    order,badorder=0,-1
    replicate,zero_rep,neg_rep=1,0,-1
    
    @test_throws ArgumentError Chain_ID(obs_id, zero_k, order, replicate)
    @test_throws ArgumentError Chain_ID(obs_id, neg_K, order, replicate)
    @test_throws ArgumentError Chain_ID(obs_id, K, badorder, replicate)
    @test_throws ArgumentError Chain_ID(obs_id, K, order, zero_rep)
    @test_throws ArgumentError Chain_ID(obs_id, K, order, neg_rep)

    id=Chain_ID(obs_id, K, order, replicate)
    @test typeof(id)==Chain_ID
    @test id.obs_id==obs_id
    @test id.K==K
    @test id.order==order
    @test id.replicate==replicate

    good_π = fill(.25,4,4)
    good_D = Categorical([.25,.25,.25,.25])
    good_Dvec = [good_D for i in 1:4]
    hmm=HMM(good_π,good_Dvec)
    log_p=-.001

    @test_throws ArgumentError EM_step(0,hmm,0.0,0.0,true)
    @test_throws ArgumentError EM_step(-1,hmm,0.0,0.0,true)
    @test_throws ArgumentError EM_step(1,hmm,1.0,0.0,true)

    steppe=EM_step(1,hmm,log_p,0.0,false)
    @test typeof(steppe)==EM_step
    @test steppe.iterate==1
    @test steppe.hmm==hmm
    @test steppe.log_p==log_p
    @test steppe.delta==0.0
    @test steppe.converged==false
end

@testset "utilities/log_prob_sum.jl Log probability summation functions" begin
    @test isapprox(lps([.1, .1, .1]), .3)
    @test lps([-Inf, .1, .1]) == -Inf
    @test isapprox(lps(.1, .1, .1), .3)
    @test lps(-Inf, .1, .1) == -Inf
    @test lps(.1, -Inf, .1) == -Inf
end

@testset "utilities/HMM_init.jl initialization functions" begin
    K=4
    order=2
    hmm=autotransition_init(K,order)
    @test typeof(hmm)==HMM{Univariate,Float64}
    @test size(hmm)==(4,64)
end

@testset "utilities/load_balancer.jl EM_converge load balancing structs and functions" begin
    @test_throws ArgumentError LoadConfig(0:5,1:5,"")
    @test_throws ArgumentError LoadConfig(1:5,-1:5,"")

    job_ids=[Chain_ID("test",6,2,1),Chain_ID("test",4,2,1),Chain_ID("test",2,0,1)]
    obs_sets=Dict("test"=>[LongSequence{DNAAlphabet{2}}("AAAAAAAAAAAAAA")])

    K_exclusive_config=LoadConfig(1:1,0:0,[""])
    O_exclusive_config=LoadConfig(1:6,3:4,[""])
    blacklist_config=LoadConfig(1:6,0:2,["test"])
    high_config=LoadConfig(4:6,0:2,[""])

    no_input_hmms, chains, input_channel, output_channel = setup_EM_jobs!(job_ids,obs_sets)
    @test load_balancer(no_input_hmms, input_channel, K_exclusive_config) == (0,0,0,0,0)

    no_input_hmms, chains, input_channel, output_channel = setup_EM_jobs!(job_ids,obs_sets)
    @test load_balancer(no_input_hmms, input_channel, O_exclusive_config) == (0,0,0,0,0)

    no_input_hmms, chains, input_channel, output_channel = setup_EM_jobs!(job_ids,obs_sets)
    @test load_balancer(no_input_hmms, input_channel, blacklist_config) == (0,0,0,0,0)

    no_input_hmms, chains, input_channel, output_channel = setup_EM_jobs!(job_ids,obs_sets)
    high_set=load_balancer(no_input_hmms, input_channel, high_config) 
    
    @test high_set[1]==Chain_ID("test",6,2,1)
    @test high_set[2]==1
    @test typeof(high_set[3])==HMM{Univariate,Float64}
    @test high_set[4]==0.0
    @test typeof(high_set[5])==Matrix{UInt8}
end

@testset "utilities/observation_coding.jl Order coding structs and functions" begin
    compound_alphabet=CompoundAlphabet(ACGT, 2)
    @test typeof(compound_alphabet.symbols)==Dict{Mer{DNAAlphabet{2}},Integer}
    @test length(compound_alphabet.symbols)==64
    @test compound_alphabet.symbols[Mer{DNAAlphabet{2}}("AAA")]==1
    @test compound_alphabet.symbols[Mer{DNAAlphabet{2}}("TTT")]==64
    
    test_seqs = [BioSequences.LongSequence{DNAAlphabet{2}}("ACGTACGTACGTACGT"),BioSequences.LongSequence{DNAAlphabet{2}}("TTTTTTT")]
    target0= [1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 0; 4 4 4 4 4 4 4 0 0 0 0 0 0 0 0 0 0]
    target2= [37 58 15 20 37 58 15 20 37 58 15 20 37 58 0; 64 64 64 64 64 0 0 0 0 0 0 0 0 0 0]

    order0_seqs = get_order_n_seqs(test_seqs,0)
    @test typeof(order0_seqs)==BioBackgroundModels.N_Order_ntSequence
    for (mer, code) in CompoundAlphabet(ACGT,0).symbols
        @test order0_seqs.alphabet.symbols[mer]==code
    end
    @test order0_seqs.seq_lengths==[16,7]
    @test order0_seqs.order_kmers[1][1]==Mer{DNAAlphabet{2}}("A")
    code0_seqs = code_seqs(order0_seqs)
    @test target0 == code0_seqs

    order2_seqs = get_order_n_seqs(test_seqs,2)
    @test typeof(order2_seqs)==BioBackgroundModels.N_Order_ntSequence
    for (mer, code) in CompoundAlphabet(ACGT,2).symbols
        @test order2_seqs.alphabet.symbols[mer]==code
    end
    @test order2_seqs.seq_lengths==[16,7]
    @test order2_seqs.order_kmers[1][1]==Mer{DNAAlphabet{2}}("ACG")
    code2_seqs = code_seqs(order2_seqs)
    @test target2 == code2_seqs

    bw_o0_seqs = get_order_n_seqs(test_seqs[end:-1:1], 0)
    sorted_code_seqs = code_seqs(bw_o0_seqs, sorted=true)
    @test target0 == sorted_code_seqs

    obs_sets=Dict("test1"=>test_seqs,
                    "test2"=>test_seqs)
    job_ids=[Chain_ID("test1",4,0,1),Chain_ID("test2",4,0,1),Chain_ID("test1",2,2,1)]

    code_dict=code_job_obs(job_ids,obs_sets)
    @test ("test1",0) in keys(code_dict)
    @test ("test1",2) in keys(code_dict)
    @test ("test2",0) in keys(code_dict)
    @test !(("test2",2) in keys(code_dict))
    @test code_dict[("test1",0)]==target0
    @test code_dict[("test1",2)]==target2
    @test code_dict[("test2",0)]==target0
end

@testset "genome_sampling/partition_masker.jl functions" begin
    synthetic_seq = BioSequences.LongSequence{DNAAlphabet{2}}(generate_synthetic_seq())
    position_start = 501
    position_length=141
    position_pad=350
    perigenic_pad = 250

    position_df = make_padded_df(posfasta, gff, genome, index, position_pad)
    add_partition_masks!(position_df, gff, perigenic_pad)

    @test position_df.SeqID[1]=="1"
    @test length(position_df.PadSeq[1]) == position_pad+position_length == length(position_df.PadStart[1]:position_df.PadStart[1]+position_length+position_pad-1) == position_df.End[1]-position_df.Start[1]+1+position_pad == size(position_df.MaskMatrix[1])[1]
    @test synthetic_seq[position_df.PadStart[1]:position_df.End[1]]==position_df.PadSeq[1]

    mm = position_df.MaskMatrix[1]
    idx=1#pos 151
    @test sum(mm[idx:idx+99,1] .== 1)==length(mm[idx:idx+99,1])
    idx+=100#pos 251
    @test sum(mm[idx:idx+259,1] .== 2)==length(mm[idx:idx+259,1])
    idx+=260 #pos 511
    @test sum(mm[idx:idx+59,1] .== 3)==length(mm[idx:idx+59,1])
    idx+=60#pos 571
    @test sum(mm[idx:idx+29,1] .== 2)==length(mm[idx:idx+29,1])
    idx+=30#pos601
    @test sum(mm[idx:idx+40,1] .== 3)==length(mm[idx:idx+40,1])
end

@testset "genome_sampling/sequence_sampler.jl functions" begin
    partitions = 3 #exonic, periexonic, intragenic
    sample_set_length=100
    min_sample_window=5
    max_sample_window=25
    perigenic_pad=250
    syn_intergenic_starts = [1]
    syn_intergenic_ends = [250]
    syn_periexonic_starts = [251,571,661,761,861,961]
    syn_periexonic_ends = [510,600,700,800,900,1000]
    syn_periexonic_strands = ['-','-','-','-','-','-']
    syn_exonic_starts = [511,601,701,801,901]
    syn_exonic_ends = [570,660,760,860,960]
    syn_exonic_strands = ['-','-','-','-','-']
    partition_lengths = Dict("exon"=>5*60,"intergenic"=>250,"periexonic"=>450)

    # "Testing sequence sampler fns..."
    synthetic_seq = BioSequences.LongSequence{DNAAlphabet{2}}(generate_synthetic_seq())

    # "Partitioning synthetic genome coordinates..."
    coordinate_partitions = partition_genome_coordinates(gff, perigenic_pad)

    @test coordinate_partitions["intergenic"].Start == syn_intergenic_starts
    @test coordinate_partitions["intergenic"].End == syn_intergenic_ends

    @test coordinate_partitions["periexonic"].Start == syn_periexonic_starts
    @test coordinate_partitions["periexonic"].End == syn_periexonic_ends
    @test coordinate_partitions["periexonic"].Strand == syn_periexonic_strands

    @test coordinate_partitions["exon"].Start == syn_exonic_starts
    @test coordinate_partitions["exon"].End == syn_exonic_ends
    @test coordinate_partitions["exon"].Strand == syn_exonic_strands

    # "Checking sampling functions at all synthetic indices..."

    input_test_channel, completed_test_channel = setup_sample_jobs(genome, index, gff, sample_set_length, min_sample_window, max_sample_window, perigenic_pad; deterministic=true)
    while isready(input_test_channel)
        genome_path, genome_index_path, partition_df, partitionid, sample_set_length, sample_window_min, sample_window_max, deterministic = take!(input_test_channel)

        for feature in eachrow(partition_df)
            seqid, scaffold_start, scaffold_end = meta_to_feature_coord(feature.MetaStart, feature.MetaEnd, partition_df)
            @test scaffold_start == feature.Start && scaffold_end == feature.End
        end

        stranded = get_strand_dict()[partitionid]
        @test typeof(stranded) == Bool

        scaffold_sequence_record_dict = build_scaffold_seq_dict(genome, index)
        @test scaffold_sequence_record_dict["1"] == synthetic_seq

        partition_length = partition_df.MetaEnd[end]
        @test partition_length == partition_lengths[partitionid]

        metacoordinate_bitarray = trues(partition_df.MetaEnd[end])

        for bitindex in findall(metacoordinate_bitarray)
            feature_metaStart, feature_metaEnd, strand = get_feature_params_from_metacoord(bitindex, partition_df, stranded)
            @test 1 <= feature_metaStart < feature_metaEnd <= length(metacoordinate_bitarray)
            @test feature_metaStart in partition_df.MetaStart
            @test feature_metaEnd in partition_df.MetaEnd
            if partitionid == "exon" || partitionid == "periexonic"
                @test strand == '-'
            end

            feature_length = length(feature_metaStart:feature_metaEnd)
            window = determine_sample_window(feature_metaStart, feature_metaEnd, bitindex, metacoordinate_bitarray, sample_window_min, sample_window_max) #get an appropriate sampling window around the selected index, given the feature boundaries and params
            @test 1 <= feature_metaStart <= window[1] < window[1]+sample_window_min-1 < window[1]+sample_window_max-1<= window[2] <= feature_metaEnd <= length(metacoordinate_bitarray)
            sample_scaffoldid, sample_scaffold_start, sample_scaffold_end = meta_to_feature_coord(window[1],window[2],partition_df)
            @test sample_scaffoldid == "1"
            @test 1 <= sample_scaffold_start <= sample_scaffold_start+sample_window_min <= sample_scaffold_end <= min(sample_scaffold_start+sample_window_max,1000) <= 1000

            strand == '-' ? target_seq=reverse_complement(synthetic_seq[sample_scaffold_start:sample_scaffold_end]) :
                target_seq = synthetic_seq[sample_scaffold_start:sample_scaffold_end]

            proposal_sequence = fetch_sequence(sample_scaffoldid, scaffold_sequence_record_dict, sample_scaffold_start, sample_scaffold_end, strand; deterministic=deterministic) #get the sequence associated with the sample window

            @test proposal_sequence == target_seq
        end
    end

    # "Verifying sampling channels..."

    input_sample_channel, completed_sample_channel = setup_sample_jobs(genome, index, gff, sample_set_length, min_sample_window, max_sample_window, perigenic_pad; deterministic=true)
    progress_channel = RemoteChannel(()->Channel{Tuple}(20))
    get_sample_set(input_sample_channel, completed_sample_channel, progress_channel)

    #collect sample dfs by partition id when ready
    collected_counter = 0
    while collected_counter < partitions
        wait(completed_sample_channel)
        partition_id, sample_df = take!(completed_sample_channel)
        sample_record_dfs[partition_id] = sample_df
        collected_counter += 1
    end

    
    synthetic_seq = BioSequences.LongSequence{DNAAlphabet{2}}(generate_synthetic_seq())
    for (partid, df) in sample_record_dfs
        for sample in eachrow(df)
            target_seq = synthetic_seq[sample.SampleStart:sample.SampleEnd]
            strand = sample.Strand
            if sample.Strand == '-'
                target_seq = reverse_complement(target_seq)
            end
            @test sample.SampleSequence == target_seq
        end
    end
end

@testset "likelihood_funcs/hmm.jl & bg_lh_matrix.jl functions" begin
    pvec = [.4,.3,.2,.1]
    trans = ones(1,1)
    D = [Categorical(pvec)]
    hmm = HMM(trans, D)

    testseq=zeros(Int64,1,5)
    testseq[1:4] = [1,2,3,4]
    @test isapprox(get_BGHMM_symbol_lh(testseq, hmm)[1], log.(pvec[1]))
    @test isapprox(obs_lh_given_hmm(testseq,hmm),obs_lh_given_hmm(testseq,hmm,linear=false))
    @test isapprox(sum(get_BGHMM_symbol_lh(testseq,hmm)), obs_lh_given_hmm(testseq,hmm))
    @test isapprox(sum(get_BGHMM_symbol_lh(testseq,hmm)), obs_lh_given_hmm(testseq,hmm,linear=false))

    pvec=[.25,.25,.25,.25]
    trans = [.9 .1
         .1 .9]
    D = [Categorical(pvec), Categorical(pvec)]
    hmm = HMM(trans, D)

    @test isapprox(obs_lh_given_hmm(testseq,hmm),obs_lh_given_hmm(testseq,hmm,linear=false))
    @test isapprox(sum(get_BGHMM_symbol_lh(testseq,hmm)), obs_lh_given_hmm(testseq,hmm))
    @test isapprox(sum(get_BGHMM_symbol_lh(testseq,hmm)), obs_lh_given_hmm(testseq,hmm,linear=false))

    testseq=zeros(Int64,1,1001)
    testseq[1,1:1000]=rand(1:4,1000)
    
    @test isapprox(sum(get_BGHMM_symbol_lh(testseq, hmm)),obs_lh_given_hmm(testseq,hmm))

    Dex = [Categorical([.3, .1, .3, .3]),Categorical([.15, .35, .35, .15])]
    Dper = [Categorical([.15, .35, .35, .15]),Categorical([.4, .1, .1, .4])]
    Dint = [Categorical([.4, .1, .1, .4]),Categorical([.45, .05, .05, .45])]
    BGHMM_dict = Dict{String,Tuple{HMM, Int64, Float64}}()
    BGHMM_dict["exon"] = (HMM(trans, Dex), 0, 0)
    BGHMM_dict["periexonic"] = (HMM(trans, Dper), 0, 0)
    BGHMM_dict["intergenic"] = (HMM(trans, Dint), 0, 0)

    position_length=141;perigenic_pad=250;
    position_df = make_padded_df(posfasta, gff, genome, index, position_length)

    reader=open(FASTA.Reader, genome, index=index)
    seq=FASTA.sequence(reader["CM002885.2.1"])
    @test position_df.Start[1]==501
    @test position_df.End[1]==641
    @test position_df.PadSeq[1]==seq[360:641]
    @test position_df.PadStart[1]==360
    @test position_df.RelStart[1]==141
    @test position_df.SeqOffset[1]==0

    offset_position_df = make_padded_df(posfasta, gff, genome, index,600)
    @test offset_position_df.Start[1]==501
    @test offset_position_df.End[1]==641
    @test offset_position_df.PadSeq[1]==seq[1:641]
    @test offset_position_df.PadStart[1]==1
    @test offset_position_df.RelStart[1]==500
    @test offset_position_df.SeqOffset[1]==100

    add_partition_masks!(position_df, gff, perigenic_pad)
    @test all(position_df.MaskMatrix[1][:,2].==-1)
    @test all(position_df.MaskMatrix[1][1:151,1].==2)
    @test all(position_df.MaskMatrix[1][152:211,1].==3)
    @test all(position_df.MaskMatrix[1][212:241,1].==2)
    @test all(position_df.MaskMatrix[1][242:end,1].==3)

    println("PadSeq: $(position_df.PadSeq)")
    println("Mask: $(position_df.MaskMatrix)")

    lh_matrix=BGHMM_likelihood_calc(position_df,BGHMM_dict)
    @test size(lh_matrix)==(position_length*2,1)

    periexonic_frag=LongSequence{DNAAlphabet{2}}(reverse_complement(seq[360:510]))
    pno=get_order_n_seqs([periexonic_frag],0)
    pcode=code_seqs(pno)
    plh=get_BGHMM_symbol_lh(pcode, BGHMM_dict["periexonic"][1])
    println("plh: $plh")
    println("lh_matrix[1:151,1]: $(lh_matrix[1:151,1])")

    @test isapprox(lh_matrix[1:151,1], plh)

    exonic_frag=LongSequence{DNAAlphabet{2}}(reverse_complement(seq[511:570]))
    eno=get_order_n_seqs([exonic_frag],0)
    ecode=code_seqs(eno)
    elh=get_BGHMM_symbol_lh(ecode, BGHMM_dict["exon"][1])
    @test isapprox(lh_matrix[152:211,1],elh)
end

@testset "mle_step functions" begin
    # "Setting up for MLE function tests.."
    π = fill((1/6),6,6)
    D = [Categorical(ones(4)/4), Categorical([.7,.05,.15,.1]),Categorical([.15,.35,.4,.1]), Categorical([.6,.15,.15,.1]),Categorical([.1,.4,.4,.1]), Categorical([.2,.2,.3,.3])]
    hmm = HMM(π, D)
    log_π = log.(hmm.π)

    obs = zeros(Int64,1,250)
    obs[1:249] = rand(1:4,249)
    obs_lengths=[249]
    # "Testing mle_step..."

    #verify that above methods independently produce equivalent output, and that this is true of multiple identical obs, but not true of different obs sets
    mouchet_hmm = mouchet_mle_step(hmm, obs[1:249])

    new_hmm = linear_step(hmm, obs, obs_lengths)

    ms_sng = BioBackgroundModels.bw_step(hmm, Array(transpose(obs)), obs_lengths)

    dblobs = zeros(Int64, 2,250)
    dblobs[1,1:249] = obs[1:249]
    dblobs[2,1:249] = obs[1:249]
    dblobs_lengths=[249,249]
    dbl_hmm = linear_step(hmm, dblobs, dblobs_lengths)


    ms_dbl =BioBackgroundModels.bw_step(hmm, Array(transpose(dblobs)),dblobs_lengths)

    otherobs = deepcopy(dblobs)
    otherobs[2,1:249] = rand(1:4,249)
    if otherobs[1,1]==otherobs[2,1]
        otherobs[1,1]==1&&(otherobs[2,1]=2)
        otherobs[1,1]==2&&(otherobs[2,1]=3)
        otherobs[1,1]==3&&(otherobs[2,1]=4)
        otherobs[1,1]==4&&(otherobs[2,1]=1)
    end
    
    other_hmm = linear_step(hmm, otherobs, dblobs_lengths)

    ms_hmm = BioBackgroundModels.bw_step(hmm, Array(transpose(otherobs)),dblobs_lengths)

    for n in fieldnames(typeof(new_hmm[1]))
        if n == :D
            for (d, dist) in enumerate(D)
            @test isapprox(new_hmm[1].D[d].support,mouchet_hmm[1].D[d].support)
            @test isapprox(new_hmm[1].D[d].p,mouchet_hmm[1].D[d].p)
            @test new_hmm[1].D[d].support==ms_sng[1].D[d].support
            @test isapprox(new_hmm[1].D[d].p,ms_sng[1].D[d].p)
            @test new_hmm[1].D[d].support==dbl_hmm[1].D[d].support
            @test isapprox(new_hmm[1].D[d].p,dbl_hmm[1].D[d].p)
            @test new_hmm[1].D[d].support==ms_dbl[1].D[d].support
            @test isapprox(new_hmm[1].D[d].p,ms_dbl[1].D[d].p)
            @test new_hmm[1].D[d].support==other_hmm[1].D[d].support
            @test ms_hmm[1].D[d].support==other_hmm[1].D[d].support
            @test !isapprox(new_hmm[1].D[d].p, other_hmm[1].D[d].p)
            @test isapprox(ms_hmm[1].D[d].p, other_hmm[1].D[d].p, atol=.005)

            end
        else
            @test isapprox(getfield(new_hmm[1],n), getfield(mouchet_hmm[1],n))
            @test isapprox(getfield(new_hmm[1],n), getfield(ms_sng[1],n))

            @test isapprox(getfield(new_hmm[1],n), getfield(dbl_hmm[1],n))
            @test isapprox(getfield(new_hmm[1],n), getfield(ms_dbl[1],n))

            @test !isapprox(getfield(new_hmm[1],n), getfield(other_hmm[1],n))
            @test isapprox(getfield(ms_hmm[1],n), getfield(other_hmm[1],n),atol=.01)

        end
    end

    @test ms_sng[2] == new_hmm[2] == mouchet_hmm[2] != dbl_hmm[2] == ms_dbl[2] != other_hmm[2] == ms_hmm[2]

    # "Testing fit_mle!..."

    #test fit_mle! function
    input_hmms= RemoteChannel(()->Channel{Tuple}(1))
    output_hmms = RemoteChannel(()->Channel{Tuple}(30))
    chainid=Chain_ID("Test",6,0,1)
    put!(input_hmms, (chainid, 2, hmm, 0.0, obs))
    BioBackgroundModels.EM_converge!(input_hmms, output_hmms, 1; max_iterates=4, verbose=true)
    wait(output_hmms)
    workerid, jobid, iterate, hmm3, log_p, delta, converged = take!(output_hmms)
    @test jobid == chainid
    @test iterate == 3
    @test BioBackgroundModels.assert_hmm(hmm3.π0, hmm3.π, hmm3.D)
    @test size(hmm3) == size(hmm) == (6,4)
    @test log_p < 1
    @test log_p == linear_likelihood(obs, hmm)
    @test converged == false
    wait(output_hmms)
    workerid, jobid, iterate, hmm4, log_p, delta, converged = take!(output_hmms)
    @test jobid == chainid
    @test iterate == 4
    @test BioBackgroundModels.assert_hmm(hmm4.π0, hmm4.π, hmm4.D)
    @test size(hmm4) == size(hmm) == (6,4)
    @test log_p < 1
    @test log_p == linear_likelihood(obs, hmm3)
    @test converged == false

    # "Test convergence.."
    obs=zeros(Int64, 4, 1001)    
    for o in 1:size(obs)[1]
        obsl=rand(100:1000)
        obs[o,1:obsl]=rand(1:4,obsl)
    end
    input_hmms= RemoteChannel(()->Channel{Tuple}(1))
    output_hmms = RemoteChannel(()->Channel{Tuple}(30))
    put!(input_hmms, (chainid, 2, hmm, 0.0, obs))
    BioBackgroundModels.EM_converge!(input_hmms, output_hmms, 1; delta_thresh=.05, max_iterates=100, verbose=true)
    wait(output_hmms)
    workerid, jobid, iterate, hmm4, log_p, delta, converged = take!(output_hmms)
    while isready(output_hmms)
        workerid, jobid, iterate, hmm4, log_p, delta, converged = take!(output_hmms)
    end
    @test converged==1
end

@testset "Baum-Welch subfunctions" begin
    π = [.5 .5
         .5 .5]
    D = [Categorical(ones(4)/4), Categorical([.7,.1,.1,.1])]
    hmm = HMM(π, D)
    log_π = log.(hmm.π)

    obs = zeros(UInt8, 22,1)
    obs[1:21] = [4,3,3,2,3,2,1,1,2,3,3,3,4,4,2,3,2,3,4,3,2]
    obs_lengths=[21]

    lls = BioBackgroundModels.bw_llhs(hmm,obs)
    m_lls = mouchet_log_likelihoods(hmm,obs[1:21])

    K,Tmaxplus1,O = size(lls)
    T=Tmaxplus1-1

    #TEST LOG LIKELIHOOD FN
    @test transpose(lls[:,1:end-1,1]) == m_lls
    
    log_α = BioBackgroundModels.messages_forwards_log(hmm.π0, hmm.π, lls,  obs_lengths)
    m_log_α = mouchet_messages_forwards_log(hmm.π0, hmm.π, m_lls)

    #TEST FORWARD MESSAGES FN
    @test transpose(log_α[:,1:end-1,1]) == m_log_α

    log_β = BioBackgroundModels.messages_backwards_log(hmm.π, lls,  obs_lengths)
    m_log_β = mouchet_messages_backwards_log(hmm.π, m_lls)

    #TEST BACKWARDS  MESSAGES FN
    @test isapprox(transpose(log_β[:,1:end-1,1]),  m_log_β)

    lls = permutedims(lls, [2,3,1]) # from (K,T,O) to (T,O,K)
    log_α = permutedims(log_α, [2,3,1])
    log_β = permutedims(log_β, [2,3,1])

    # "Testing obs probability calcs..."
    #TEST OBSERVATION PROBABILITY CALC
    normalizer = logsumexp(m_log_α[1,:] + m_log_β[1,:])
    log_pobs = logsumexp(lps.(log_α[1,1,:], log_β[1,1,:]))

    @test isapprox(normalizer, log_pobs)

    # "Testing ξ, γ calcs..."

    log_ξ = fill(-Inf, Tmaxplus1,O,K,K)
    log_γ = fill(-Inf, Tmaxplus1,O,K)

    m_log_ξ = zeros(T, K, K)

    for t = 1:T-1, i = 1:K, j = 1:K
        m_log_ξ[t,i,j] = m_log_α[t,i] + log_π[i,j] + m_log_β[t+1,j] + m_lls[t+1,j] - normalizer
    end

    for i = 1:K, j = 1:K, o = 1:O
        obsl = obs_lengths[o]
        for t = 1:obsl-1 #log_ξ & log_γ calculated to T-1 for each o
           log_ξ[t,o,i,j] = lps(log_α[t,o,i], log_π[i,j], log_β[t+1,o,j], lls[t+1,o,j], -log_pobs[o])
           log_γ[t,o,i] = lps(log_α[t,o,i], log_β[t,o,i], -log_pobs[o])
        end
        t=obsl #log_ξ @ T = 0
        log_ξ[t,o,i,j] = 0
        log_γ[t,o,i] = lps(log_α[t,o,i], log_β[t,o,i], -log_pobs)
    end

    ξ = exp.(log_ξ)
    k_ξ = sum(ξ, dims=[3,4])
    nan_mask = k_ξ .== 0; k_ξ[nan_mask] .= Inf #prevent NaNs in dummy renorm arising from multiple sequences indexing
    ξ  ./=  k_ξ #dummy renorm across K to keep numerical creep from causing isprobvec to fail on new new_π during hmm creation

    m_ξ = exp.(m_log_ξ)
    m_ξ ./= sum(m_ξ, dims=[2,3])

    #check equivalency of xi calc methods
    @test isapprox(reshape(ξ,Tmaxplus1,K,K)[1:21,:,:],m_ξ)

    γ = exp.(log_γ)
    k_γ = sum(γ, dims=3); k_γ[nan_mask[:,:,:]] .= Inf #prevent NaNs in dummy renorm
    γ ./= k_γ

    m_γ = exp.((m_log_α .+ m_log_β) .- normalizer)

    #check equivalency of gamma calculations
    @test isapprox(reshape(γ,Tmaxplus1,K)[1:21,:,:],m_γ)

    # "Testing initial and transition matrix calcs..."

    m_new_π = sum(m_ξ[1:end-1,:,:], dims=1)[1,:,:]
    m_new_π ./= sum(m_new_π, dims=2)

    new_π = zeros(K,K)
    for i=1:K, j=1:K
        ∑otξ_vec = zeros(O)
        ∑otγ_vec = zeros(O)
        for o in 1:O
            ∑otξ_vec[o] = sum(ξ[1:obs_lengths[o]-1,o,i,j])
            ∑otγ_vec[o] = sum(γ[1:obs_lengths[o]-1,o,i])
        end
        new_π[i,j] = sum(∑otξ_vec) / sum(∑otγ_vec)
    end
    new_π ./= sum(new_π, dims=[2]) #dummy renorm

    #check equivalency of transition matrix calculations
    @test isapprox(m_new_π,new_π)

    m_new_π0 = exp.((m_log_α[1,:] + m_log_β[1,:]) .- normalizer)
    m_new_π0 ./= sum(m_new_π0)

    new_π0 = (sum(γ[1,:,:], dims=1)./sum(sum(γ[1,:,:], dims=1)))[1,:]
    new_π0 ./= sum(new_π0) #dummy renorm

    @test isapprox(m_new_π0,new_π0)

    # "Testing distribution calcs..."

    F=Univariate
    m_D = Distribution{F}[]
    for (i, d) in enumerate(hmm.D)
        # Super hacky...
        # https://github.com/JuliaStats/Distributions.jl/issues/809
        push!(m_D, fit_mle(eval(typeof(d).name.name), permutedims(obs[1:21]), m_γ[:,i]))
    end

    obs_mask = .!nan_mask
    obs_collection = obs[obs_mask[:,:]]

    D = Distribution{F}[]
    @inbounds for (i, d) in enumerate(hmm.D)
        # Super hacky...
        # https://github.com/JuliaStats/Distributions.jl/issues/809
        γ_d = γ[:,:,i]
        push!(D, fit_mle(eval(typeof(d).name.name), obs_collection, γ_d[obs_mask[:,:]]))
        #slowest call by orders of magnitude
    end

    #test maximization equivalency
    for (d, dist) in enumerate(D)
        @test isapprox(dist.support, m_D[d].support)
        @test isapprox(dist.p, m_D[d].p)
    end

    # "Testing mle_step..."

    #verify that above methods independently produce equivalent output, and that this is true of multiple identical obs, but not true of different obs sets
    mouchet_hmm = mouchet_mle_step(hmm, obs[1:21])

    new_hmm = BioBackgroundModels.bw_step(hmm, obs, obs_lengths)

    dblobs = zeros(UInt8, 22,2)
    dblobs[1:21,1] = [4,3,3,2,3,2,1,1,2,3,3,3,4,4,2,3,2,3,4,3,2]
    dblobs[1:21,2] = [4,3,3,2,3,2,1,1,2,3,3,3,4,4,2,3,2,3,4,3,2]
    dblobs_lengths=[21,21]
    dbl_hmm = BioBackgroundModels.bw_step(hmm, dblobs, dblobs_lengths)

    otherobs = dblobs
    otherobs[1:21,2] = [1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1]
    other_hmm = BioBackgroundModels.bw_step(hmm, otherobs, dblobs_lengths)

    for n in fieldnames(typeof(new_hmm[1]))
        if n == :D
            for (d, dist) in enumerate(D)
                @test isapprox(dist.support,mouchet_hmm[1].D[d].support)
                @test isapprox(dist.p,mouchet_hmm[1].D[d].p)
                @test isapprox(dist.support,dbl_hmm[1].D[d].support)
                @test isapprox(dist.p,dbl_hmm[1].D[d].p)
                @test isapprox(dist.support, other_hmm[1].D[d].support)
                @test !isapprox(dist.p, other_hmm[1].D[d].p)
            end
        else
            @test isapprox(getfield(new_hmm[1],n), getfield(mouchet_hmm[1],n))
            @test isapprox(getfield(new_hmm[1],n), getfield(dbl_hmm[1],n))
            @test !isapprox(getfield(new_hmm[1],n), getfield(other_hmm[1],n))
        end
    end

    @test new_hmm[2] == mouchet_hmm[2] != dbl_hmm[2] != other_hmm[2]

    # "Testing fit_mle!..."

    #test fit_mle! function
    input_hmms= RemoteChannel(()->Channel{Tuple}(1))
    output_hmms = RemoteChannel(()->Channel{Tuple}(30))
    chainid=Chain_ID("Test",2,0,1)
    put!(input_hmms, (chainid, 2, hmm, 0.0, transpose(obs)))
    BioBackgroundModels.EM_converge!(input_hmms, output_hmms, 1; EM_func=BioBackgroundModels.bw_step, max_iterates=4, verbose=true)
    wait(output_hmms)
    workerid, jobid, iterate, hmm3, log_p, epsilon, converged = take!(output_hmms)
    @test jobid == chainid
    @test iterate == 3
    @test assert_hmm(hmm3.π0, hmm3.π, hmm3.D)
    @test size(hmm3) == size(hmm) == (2,4)
    @test log_p < 1
    @test log_p == obs_lh_given_hmm(transpose(obs),hmm, linear=false)
    @test converged == false
    wait(output_hmms)
    workerid, jobid, iterate, hmm4, log_p, epsilon, converged = take!(output_hmms)
    @test jobid == chainid
    @test iterate == 4
    @test assert_hmm(hmm4.π0, hmm4.π, hmm4.D)
    @test size(hmm4) == size(hmm) == (2,4)
    @test log_p < 1
    @test log_p == obs_lh_given_hmm(transpose(obs),hmm3,linear=false)
    @test converged == false

    # "Test convergence.."
    obs=zeros(UInt8, 101, 2)
    for i in 1:size(obs)[2]
        obs[1:100,i]=rand(1:4,100)
    end
    input_hmms= RemoteChannel(()->Channel{Tuple}(1))
    output_hmms = RemoteChannel(()->Channel{Tuple}(30))
    put!(input_hmms, (chainid, 2, hmm, 0.0, transpose(obs)))
    BioBackgroundModels.EM_converge!(input_hmms, output_hmms, 1; EM_func=BioBackgroundModels.bw_step, delta_thresh=.05, max_iterates=100, verbose=true)
    wait(output_hmms)
    workerid, jobid, iterate, hmm4, log_p, epsilon, converged = take!(output_hmms)
    while isready(output_hmms)
        workerid, jobid, iterate, hmm4, log_p, epsilon, converged = take!(output_hmms)
    end
    @test converged==1
end

@testset "setup_EM_jobs API" begin
    #CONSTANTS FOR HMM LEARNING
    replicates = 4 #repeat optimisation from this many seperately initialised samples from the prior
    Ks = [1,2,4,6] #mosaic class #s to test
    order_nos = [0,1,2] #DNA kmer order #s to test
    training_sets, test_sets = split_obs_sets(sample_record_dfs)


    job_ids=Vector{Chain_ID}()
    for (obs_id, obs) in training_sets, K in Ks, order in order_nos, rep in 1:replicates
        push!(job_ids, Chain_ID(obs_id, K, order, rep))
    end

    no_input_hmms, chains, input_hmms, output_hmms = setup_EM_jobs!(job_ids, training_sets)
    @test no_input_hmms == replicates*length(Ks)*length(order_nos)*length(training_sets)

    while isready(input_hmms)
        jobid, start_iterate, hmm, last_norm, observations = take!(input_hmms)
        @test last_norm == 0
        obs_lengths = [findfirst(iszero,observations[o,:])-1 for o in 1:size(observations)[1]]
        #make sure input HMMs are valid and try to mle_step them and ensure their 1-step children are valid
        @test assert_hmm(hmm.π0, hmm.π, hmm.D)
        new_hmm, prob = linear_step(hmm,observations,obs_lengths)
        @test assert_hmm(hmm.π0, hmm.π, hmm.D)
        @test prob < 0
    end
end

@testset "BBM API tests..." begin
    genome_reader = open(FASTA.Reader, genome, index=index)
    seq=LongSequence{DNAAlphabet{2}}(FASTA.sequence(genome_reader["CM002885.2.1"]))
    seqdict = Dict("test"=>[seq for i in 1:3])


    job_ids=Vector{Chain_ID}()
    Ks=[1,2]; order_nos=[0,1]
    for K in Ks, order in order_nos
        push!(job_ids, Chain_ID("test", K, order, 1))
    end

    wkpool=addprocs(2, topology=:master_worker)
    @everywhere using BioBackgroundModels

    em_jobset=setup_EM_jobs!(job_ids, seqdict)
    execute_EM_jobs!(wkpool, em_jobset..., "testchains", delta_thresh=.1, verbose=true)

    rm("testchains")
end