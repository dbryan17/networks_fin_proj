include("util.jl")
using Graphs

#=
This file will be for random graph generation and getting the empricial graphs. 
  

Notes: 
- The paper used: 
  + Erdos Reney x 
  + Watts-Strogatz x 
  + Barbasbasi Albert x 
  ----------------------
  + Twitter Networks 
  + Karate Club 
- I want to add
  + DC-SBM 
  + SBM 
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
planted paritions sbm model as described in the notes for assortiative mixing
n - number of nodes 
d - mean degree
c - number of communities
e - how strong of assortive mixing. 0 is ER graph, 2*d is all within communities
=#
function sbm_pp_asoritive(n :: Int, d :: Int, c :: Int, e :: Float64) :: Graphs.SimpleGraph{Int64} 
  # p - along the diagonoal probablity
  p = (d + e/c)/n 
  # q - off diagonol probablity
  q = (d - e/c)/n
  mm = [i == j ? p : q for i in 1:c, j in 1:c]
  n_per_c :: Int = floor(n / c)
  leftover :: Int = n % c
  # partition
  z :: Vector{Int} = [j for i in 1:n_per_c for j in 1:c]

  # add leftovers
  for i in 1:leftover
    push!(z,i)
  end

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
function sbm_pp_cp(n :: Int, d :: Int, c :: Int, r :: Float64, r0 :: Float64) :: Graphs.SimpleGraph{Int64} 

  d :: Float64 = d * 2

  # calc sizes of each group
  sizes = foldl(
    (acc, _) -> push!(acc, floor(Int, acc[end] * r)),
    2:c, 
    init = [floor(Int, n * (r - 1.) / (r^c - 1.))]
  )

  println(reduce(+, sizes))
  leftover = n - reduce(+, sizes)

  println(sizes)

  println(leftover)


  for i in 0:round(leftover)
    idx = c - (i % c)
    sizes[idx] += 1
  end

  println(leftover)

  println(sizes)

  z :: Vector{Int} = []

  # make groups
  start = 1 
  for (id, size) in enumerate(sizes)
    for i in start:start+size 
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



  return sbm_core(z, mm)
end

g = sbm_pp_asoritive(10000, 8, 4, 0.)
println((2 * ne(g)) / nv(g) )

g1 = sbm_pp_cp(10000, 32, 10, 2.1, 2.3)
println((2 * ne(g1)) / nv(g1) )
