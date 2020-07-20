"""
Utility functions for learning and using background genomic hidden markov models
"""
#function to split random sample dataframe into training and test sets (divide total sequence length by half)
function split_obs_sets(sample_dfs::Dict{String,DataFrame})
    training_sets = Dict{String,Vector{LongSequence{DNAAlphabet{2}}}}()
    test_sets = Dict{String,Vector{LongSequence{DNAAlphabet{2}}}}()

    for (partition_id, partition) in sample_dfs
        partition.sampleLength = (partition.SampleEnd - partition.SampleStart) .+ 1
        midway = sum(partition.sampleLength)÷2
        split_index = 0
        counter = 0
        while split_index == 0
            counter += 1
            length_sum = sum(partition.sampleLength[1:counter])
            if length_sum > midway
                split_index = counter
            end
        end

        training_sets[partition_id]  = partition.SampleSequence[1:split_index-1]
        test_sets[partition_id] = partition.SampleSequence[split_index:end]
    end
    return training_sets, test_sets
end