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

  println("n, d, c: " * string(n) * " " * string(d) * " " * string(c))
  println("len of k " * string(length(z)))
  for i in z
    println(i)
  end
  print(mm)

  return sbm_core(z, mm)


end

