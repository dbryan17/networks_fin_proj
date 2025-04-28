include("util.jl")
include("plots.jl")
using Graphs
using GraphPlot
using Plots
using Distributions
#=
This file will be for random graph generation and getting the empricial graphs. 
  

Notes: 
- The paper used: 
  + Erdos Reney x 
  + Watts-Strogatz x 
  + Barbasbasi Albert x 
  ----------------------
  + Twitter Networks NO 
  + Karate Club TODO
- I want to add
  + FB TODO 
  + DC-SBM NO 
  + SBM x 
  + Two scale SBM (friend groups with core perhihary) x 
- My proposal said dc-sbm with friend groups with core-perhiery 
=#


function er_graph(n :: Int, p :: Float64) :: Graphs.SimpleGraph{Int64}
  return erdos_renyi(n, p)
end

# k is the nearest neighbor connection
function ws_graph(n :: Int, k :: Int, p :: Float64) :: Graphs.SimpleGraph{Int64}
  return watts_strogatz(n, k, p)
end

function ba_graph(n :: Int, k :: Int) :: Graphs.SimpleGraph{Int64}
  return barabasi_albert(n, k)
end

#=
z - partition of nodes by community membership
m - mixing matrix (stochasitc block matrx)
    m_r_s gives the probablity a node in community r is connection to a node in community s
    this is symetric
=#
function sbm_core(z :: Vector{Int}, m :: Matrix{Float64}) :: Graphs.SimpleGraph{Int64}
  g = SimpleGraph(length(z))
  for i in 1:length(z)
    # iterate through every node
    i_comm = z[i]
    for j in i+1:length(z)
      # every node "after"
      j_comm = z[j]
      p = m[i_comm, j_comm]
      if rand() < p
        add_edge!(g, i, j)
      end
    end
  end
  return g
end



#= 
this is for assortive mixing with core perhicary inside each assortive community. 
n - number of nodes 
d - mean degree
c_assort - number of communities for the outer assortative mixing 
c_cp - number of communities inside the assortative communities for the core perherairy
e - how strong of assortive mixing. 0 is ER graph, d(1 - 1/c) is all connections within community. c_assort that is
r - scale factor. How much bigger the further out community is from the next closest to the core
r0 - scale factor for degree. How much smaller the average degree of one out is. For ex, c1 = 5 c2 = 10. c1 avg degree is double that of c2 if r0 = 2
=# 

# now also returns the parition for associating with opinion
function sbm_two_scale(n :: Int, d :: Float64, c_assort :: Int, c_cp :: Int, e :: Float64, r :: Float64, r0 :: Float64) :: Tuple{Vector{Int}, Graphs.SimpleGraph{Int64}}

  # I am going for inside community core perhiary, then every node inside that communitiy has same likelyhood of connecting to a node outside the community 
  # this is easier, and I think it makes the most sense


  # get z and mm for regular assortative
  z_top, mm_top = sbm_assoritive_helper(n, d, c_assort, e)

  z_top_sorted = sort(z_top)

  # collect amount of nodes for each 
  counts :: Vector{Int} = foldl(
    ((counts, prev), curr) -> prev == curr ? (counts[end] +=1; (counts, curr)) : (push!(counts, 1), curr), 
    z_top_sorted, 
    init = (Int[], -1)
  )[1]

  counts_poss = unique(counts)

  # calculate n, d for inside c-p
  intra_p = mm_top[1, 1]

  # maps counts to average degree
  d_inners :: Dict{Int, Float64} = Dict(count => count * intra_p for count in counts_poss)

  # maps 
  z_mm_insides :: Dict{Int, Tuple{Vector{Int}, Matrix{Float64}}} = Dict(n => sbm_pp_cp_helper(n, d, c_cp, r, r0) for (n, d) in d_inners)

  num_to_count = Dict{Int, Int}()

  big_z :: Vector{Tuple{Int, Int}} = [] 

  for (i, count) in enumerate(counts)
    z, mm = z_mm_insides[count]
    num_to_count[i] = count
    big_z_part = [(i, j) for j in z]
    big_z = vcat(big_z, big_z_part)
  end

  g = SimpleGraph(n)

  for (i, (i_top, i_inner)) in enumerate(big_z)
    for j in i+1:length(big_z)
      j_top, j_inner = big_z[j]
      # if the outer communities are the same, defer to the inner communities 
      if i_top == j_top
        mm_inside = z_mm_insides[num_to_count[i_top]][2]
        p = mm_inside[i_inner, j_inner]
      # if the outer communites are not the same, the probablity of connection is the outer mm
      else
        p = mm_top[i_top, j_top]
      end
      if rand() < p
        add_edge!(g, i, j)
      end
    end
  end

  return z_top, g

end

#=
planted paritions sbm model as described in the notes for assortiative mixing
n - number of nodes 
d - mean degree
c - number of communities
e - how strong of assortive mixing. 0 is ER graph, d(1 - 1/c) is all connections within community

this returns parition and mixing matrix instead of graph
=#
function sbm_assoritive_helper(n :: Int, d :: Float64, c :: Int, e :: Float64) :: Tuple{Vector{Int}, Matrix{Float64}}

  # solve for d_in 
  d_in = d/c + e 
  d_out = d - (d/c) - e

  # solve for p and q
  p = d_in / (n/c)
  q = d_out/ (n - (n/c))

  # # p - along the diagonoal probablity
  # p = (d/c + e/2)/n
  # # q - off diagonol probablity
  # q = (d/c - e/2)/n

  mm = [i == j ? p : q for i in 1:c, j in 1:c]
  n_per_c :: Int = floor(n / c)
  leftover :: Int = n % c
  # partition
  z :: Vector{Int} = [j for i in 1:n_per_c for j in 1:c]

  # add leftovers
  for i in 1:leftover
    push!(z,i)
  end
  return z,mm
end


#=
planted paritions sbm model as described in the notes for assortiative mixing

returns graph
=#
function sbm_pp_asoritive(n :: Int, d :: Float64, c :: Int, e :: Float64) :: Graphs.SimpleGraph{Int64} 
  z, mm = sbm_assoritive_helper(n, d, c, e)
  return sbm_core(z, mm)
end


#= 
planted partition sbm model for core perhary structure
given number of nodes, number of communities, and scale factor, generate number of nodes per group
n - number of nodes 
d - mean degree
c - number of communities
r - scale factor. How much bigger the further out community is from the next closest to the core
r0 - scale factor for degree. How much smaller the average degree of one out is. For ex, c1 = 5 c2 = 10. c1 avg degree is double that of c2 if r0 = 2
e - how strong of assortive mixing. 0 is ER graph, 2*d is all within communities
=#

# could potentially edit this to be given a list of sizes and a density list 
# this is working, but need to check it because at large c values and small n values and large r and r0 values the groups are too small for it to work
# the groups need to be big enough, then the average degree works, but at small values, it doesn't really work 
function sbm_pp_cp_helper(n :: Int, d :: Float64, c :: Int, r :: Float64, r0 :: Float64) :: Tuple{Vector{Int}, Matrix{Float64}}

  d :: Float64 = d * 2

  # calc sizes of each group
  sizes_f :: Vector{Float64} = foldl(
    (acc, _) -> push!(acc, (acc[end] * r)),
    2:c, 
    init = [n * (r - 1.) / (r^c - 1.)]
  )

  sizes :: Vector{Int} = [floor(size) for size in sizes_f]

  leftover = n - reduce(+, sizes)

  for i in 0:leftover
    if i != leftover 
      idx = c - (i % c)
      sizes[idx] += 1
    end
  end

  z :: Vector{Int} = []

  # make groups
  for (id, size) in enumerate(sizes)
    for i in 1:size 
      push!(z, id)
    end
  end

  # calculate average degree for each block based on total average degree
  mult :: Float64 = 1.
  sum :: Float64 = 0.
  for size in sizes
    sum += size * (1/mult)
    mult = mult * r0
  end

  start_mean_degree :: Float64 = (d * n) / sum

  # calcuate expected num stubs for each group
  mult = 1.
  comm_stubs :: Vector{Float64} = []
  right_hand_sides :: Vector{Float64} = []
  for size in sizes
    stubs = size * (start_mean_degree/mult)
    push!(comm_stubs, stubs)
    push!(right_hand_sides, stubs / (size * 2))
    mult = mult * r0
  end

  # create and solve system of linear equations for the porabiblites 
  # going for mixing matrix like this: 
  #= 
  p0 p1 p2
  p1 p1 p2 
  p2 p2 p2
  =#

  # stubs per group conditioned on the probablity is 
  #= 
  d0 * s0 = 2 * s0 * (p0*s0 + p1*s1 + p2*s2)
  d1 * s1 = 2 * s0 * (p1*s0 + p1*s1 + p2*s2)
  ...
  =# 
  system_matrix = Matrix{Int}(undef, c, c)

  rewrite_num = 0
  for i in 1:c
    if rewrite_num > 0
      row :: Vector{Int} = []
      sum = 0 
      for j in 1:rewrite_num
        sum += sizes[j]
        push!(row, 0)
      end
      push!(row, sum + sizes[rewrite_num + 1])
      for j in rewrite_num+2:length(sizes)
        push!(row, sizes[j])
      end
      vcat 
      system_matrix[i, :] .= row
    else 
      system_matrix[i, :] .= sizes
    end
    
    rewrite_num += 1

  end

  probs_arr :: Vector{Float64} = system_matrix \ right_hand_sides


  mm = Matrix{Float64}(undef, c, c)

  for (i, p) in enumerate(probs_arr)
    for a in 1:i 
      mm[i, a] = p 
      mm[a, i] = p
    end
  end
  return z, mm
end

# returns graph
function sbm_pp_cp(n :: Int, d :: Float64, c :: Int, r :: Float64, r0 :: Float64) :: Graphs.SimpleGraph{Int64} 
  z, mm = sbm_pp_cp_helper(n, d, c, r, r0)
  sbm_core(z, mm)
end 


