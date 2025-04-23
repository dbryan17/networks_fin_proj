#= 
DeGroot model 
takes a node's opinion, neigbors opions, and self weight, returns new opinion
=#
function degroot_node_comp(n_x :: Float64, adjs_x :: Vector{Float64}, w_nn :: Float64) :: Float64
  return ((n_x * w_nn) + reduce(+, adjs_x)) / (w_nn + length(adjs_x))
end


#= 
run one time step of degroot model. 
The opinions (ops) is the attributes from 0 to 1 for each node. Represented as a vector from nodes (indices) to node attributes
Selfs is the self value as defined below

self_value = 0 - means people value their own opinion as the same as everyone else's (not sum)
self_value = 1 - means that people value there own opinion as equal to everyone else's (sum)
self_vlaue = 2 - means that people value their own opinion as double of everyone else (sum)

returns new opinions
=#
function degroot_step(ops :: Vector{Float64}, selfs :: Vector{Float64}, g :: Graphs.SimpleGraph{Int64}) :: Vector{Float64}
  ops_t1 = Vector{Float64}(undef, length(ops))
  for n in vertices(g)
    adjs = neighbors(g, n)
    adj_vals :: Vector{Float64} = [ops[i] for i in adjs]
    n_val = ops[n]
    self_val = selfs[n]
    wii :: Float64 = self_val == 0. ? 1. : self_val * length(adjs)
    new_n_val = degroot_node_comp(n_val, adj_vals, wii)
    ops_t1[n] = new_n_val
  end
  return ops_t1
end



#=
n - number of steps
g - graph
ops - initial opinions 
selfs - self values

returns array of opinions for each time step
=# 
function degroot_sim(n :: Int, g :: Graphs.SimpleGraph{Int}, ops :: Vector{Float64}, selfs :: Vector{Float64})  :: Vector{Vector{Float64}}

  if !(nv(g) == length(ops) == length(selfs)) 
    error("In Degroot Sim - lengths do not match") 
  end

  all_ops :: Vector{Vector{Float64}} = [ops]

  for _ in n
    ops = degroot_step(ops, selfs, g)
    push!(all_ops, ops)
  end

  return all_ops
end


# takes a graph, and runs degrot with full random variables from 0 to 1, with self wieght = s for all. and a number of iterations
# returns array of how the opinions progressed
function degroot_full_rand(g :: Graphs.SimpleGraph{Int}, n :: Int, s :: Float64) :: Vector{Vector{Float64}}
  ops_0 = Vector{Float64}(undef, nv(g))
  selfs = Vector{Float64}(undef, nv(g))
  for i in vertices(g)
    ops_0[i] = rand()
    selfs[i] = s

  end
  return degroot_sim(n, g, ops_0, selfs)
end


