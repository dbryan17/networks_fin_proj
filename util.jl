```
Utility functions that I use all the time
```
using Graphs

# from edge set
function get_adj_list(edges :: Set{Tuple{Int, Int}}) :: Dict{Int, Set{Int}}
  adj_list = Dict{Int, Set{Int}}()
  for (n1, n2) in edges
    if (haskey(adj_list, n1))
      push!(adj_list[n1], n2)
    else
      adj_list[n1] = Set(n2)
    end

    if (haskey(adj_list, n2))
      push!(adj_list[n2], n1)
    else
      adj_list[n2] = Set(n1)
    end
  end
  return adj_list
end

# from adj dict
function get_vertices(adj_list :: Dict{Int, Set{Int}}) :: Vector{Int}
  keyss = keys(adj_list)
  vertices :: Vector{Int} = []
  for key in keyss
    push!(vertices, key)
  end
  return vertices
end

# create julia graph from my graph
function create_jl_graph(edgeSet :: Set{Tuple{Int, Int}}, length :: Int) :: Graphs.SimpleGraph{Int64}
  g = SimpleGraph(length)

  for (u, v) in edgeSet
    add_edge!(g, u, v)
  end
  return g
end 

function create_my_graph(g :: Graphs.SimpleGraph{Int64}) :: Set{Int, Int} 
  return edges(g)
end

