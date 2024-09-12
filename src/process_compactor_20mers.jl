using Combinatorics
using DataFrames
using CSV

# trim string 
function trim_compactor(input::AbstractString)
    return input[1:length(input)-1]
end

# parse string
function parse_string(input::AbstractString, segment_length::Int)
    seg_num = Int(80/segment_length)
    segments = Vector{String}(undef, seg_num)
    for i in 1:seg_num
        start_idx = (i - 1) * segment_length + 1
        end_idx = i * segment_length
        segments[i] = input[start_idx:end_idx]
    end
    return segments
end

# combine segments based on index_comb
function string_combinations(input::AbstractString, index_comb::Vector{Vector{Int}}, segment_length::Int)
    strings = parse_string(trim_compactor(input), segment_length)
    scrambled_list = []
    for comb in index_comb
        push!(scrambled_list, (strings[comb[1]] * strings[comb[2]], strings[comb[3]] * strings[comb[4]], strings[comb[5]] * strings[comb[6]], strings[comb[7]] * strings[comb[8]]))
    end
    return scrambled_list
end

# get the rank of each segment
function sortperm_rev(vector)
    return sortperm(vector, rev=true)
end

# get base target
function get_base_target(seg_vec, rank_vec)
    return seg_vec[findmin(rank_vec)[2]]
end

# get the weight of target
function get_wgt(support_vec)
    return support_vec ./ sum(support_vec)
end

function main()
    # read in the compactor file
    df =  CSV.read(ARGS[1], DataFrame, delim='\t', header=true)

    # select compactor length of 81
    select!(df, :, [:compactor] => ByRow(length) => [:length])
    df = df[df.length .== 81, :]
    select!(df,Not(:length))

    # Add compactor rank to a new column
    transform!(groupby(df, [:anchor,]), [:support => sortperm_rev => :support_rank])

    # generate scrambled compactor segments, indexed as <Destruction><No.><Combination><No.> and stack them into a new column 'segment
    index_comb = [[1,2,3,4,5,6,7,8], 
                  [1,5,2,6,3,4,7,8], 
                  [1,8,2,7,3,6,4,5], 
                  [1,4,2,3,5,8,6,7], 
                  [1,6,2,5,3,4,7,8]]
    transform!(df, :, :compactor => ByRow(x -> string_combinations(x, index_comb, 10)) => [:D1, :D2, :D3, :D4, :D5])
    pivoted_df = stack(df, r"D")
    rename!(pivoted_df, Dict(:variable => :segment_index, :value => :segment))

    # get target weight 
    transform!(groupby(pivoted_df, [:anchor, :segment_index]), [:support => get_wgt => :target_weight])
    
    # filter compactor abundance >.05
    filter!(pivoted_df -> pivoted_df.target_weight .>= 0.05, pivoted_df)

    # add 'base_target' column and filter rows to remove those where 'base_target' is equal to 'segment'.
    transform!(groupby(pivoted_df, [:anchor, :segment_index]), [[:segment, :support_rank] => get_base_target => :base_target])
    filter!(pivoted_df -> pivoted_df.base_target != pivoted_df.segment, pivoted_df)
    
    # Recalculate 'target_weight' after filtering, grouping by 'anchor' and 'segment_index'.
    transform!(groupby(pivoted_df, [:anchor, :segment_index]), [:support => get_wgt => :target_weight])
        
    # place base_target tuple and segment tuple to seperate columns
    transform!(pivoted_df, :segment => ByRow(x -> (x[1], x[2], x[3], x[4])) => [:S1, :S2, :S3, :S4])
    transform!(pivoted_df, :base_target => ByRow(x -> (x[1], x[2], x[3], x[4])) => [:base_S1, :base_S2, :base_S3, :base_S4])
    select!(pivoted_df, Not([:segment, :base_target]))
    
    # write tsv 
    out_file = string(ARGS[2], "/", "processed_compactors_20mers.tsv")
    CSV.write(out_file, pivoted_df, delim='\t')
end

main()


