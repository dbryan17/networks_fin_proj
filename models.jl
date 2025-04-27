using Random
using Graphs



######### node comps ##########


#= 
DeGroot model 
takes a node's opinion, neigbors opions, and self weight, returns new opinion
=#
function degroot_node_comp(n_x :: Float64, adjs_x :: Vector{Float64}, w_nn :: Float64) :: Float64
  return ((n_x * w_nn) + reduce(+, adjs_x)) / (w_nn + length(adjs_x))
end


#=
BEBA model as discribed in the paper
now goes between -1 and 1
=#

function beba_node_comp(y_i :: Float64, adjs_i :: Vector{Float64}, w_ii :: Float64, b :: Float64) :: Float64
  weights :: Vector{Float64} = [b * y_i * y_j + 1 for y_j in adjs_i]

  # equation 5 in paper, to normalize the thing
  if w_ii + reduce(+, weights) <= 0
    # do the normalization thing to avoid weirdness
    return sign(y_i)
  else 
    # normal version
    weights_and_ops = [wij * yj for (wij, yj) in zip(weights, adjs_i)]
    # we clamp the nodes new value at -1 to 1, because that is the range it exists on. Because of the bias, the computation 
    # can go above 1 and below -1
    return clamp(( (w_ii * y_i) + (reduce(+, weights_and_ops)) ) / ( w_ii + reduce(+, weights) ), -1., 1.)
  end

end


#= 
my model
goes between -1 and 1
this is the full version that genralizes all the other versions
all we are chaning is the weight calcuation of beba 
=#

# b is β and g is γ 
function mine_weight_comp(b_i :: Float64, g_i_tup :: Tuple{Float64, Float64}, y_i :: Float64, y_j :: Float64) :: Float64
  # very unlikely with 64 pt prcision
  if y_i == 0.0 || y_j == 0.0
    error("in mine weight comp... zero op")
  end
  g_i, p_i = g_i_tup
  r = rand()
  prod :: Float64 = y_i * y_j
  diff_abs :: Float64 = abs(y_i - y_j)
  # backfire effect activated
  if y_i != 0. && 0. > prod && r < p_i
    return (diff_abs / 2) * (g_i * -1)
  # confirmation bias comp on the same side
  elseif 0. < prod  
    return (b_i * (1 - diff_abs)) + 1
  # confirmation bias comp on different sides
  else 
    return max( (b_i*((-1*diff_abs) / 2)) + 1 , 0)
  end
end

function mine_comp(y_i :: Float64, adjs_i :: Vector{Float64}, w_ii :: Float64, b_i :: Float64, g_i_tup :: Tuple{Float64, Float64}) :: Float64
  weights :: Vector{Float64} = [mine_weight_comp(b_i, g_i_tup, y_i, y_j) for y_j in adjs_i]
  numer :: Float64 = w_ii * y_i + reduce(+, [ w_ij * y_j for (w_ij, y_j) in zip(weights, adjs_i) ])
  denomer :: Float64 = w_ii + reduce(+, [ abs(w_ij) for w_ij in weights ])
  return numer / denomer

end


###################################################


############# step functions ######################



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


function beba_step(ops :: Vector{Float64}, selfs :: Vector{Float64}, bs :: Vector{Float64}, g :: Graphs.SimpleGraph{Int64}) :: Vector{Float64}
  ops_t1 = Vector{Float64}(undef, length(ops))
  for n in vertices(g)
    adjs :: Vector{Int} = neighbors(g, n)
    adj_vals :: Vector{Float64} = [ops[i] for i in adjs]
    n_val :: Float64 = ops[n]
    self_val :: Float64 = selfs[n]
    b :: Float64 = bs[n]
    wii :: Float64 = self_val == 0. ? 1. : self_val * length(adjs)
    new_n_val :: Float64 = beba_node_comp(n_val, adj_vals, wii, b)
    ops_t1[n] = new_n_val
  end
  return ops_t1
end


function mine_step(ops :: Vector{Float64}, selfs :: Vector{Float64}, bs :: Vector{Float64}, g_ps :: Vector{Tuple{Float64, Float64}}, g) :: Vector{Float64}
  ops_t1 = Vector{Float64}(undef, length(ops))
  for n in vertices(g)
    adjs :: Vector{Int} = neighbors(g, n)
    adj_vals :: Vector{Float64} = [ops[i] for i in adjs]
    n_val :: Float64 = ops[n]
    self_val :: Float64 = selfs[n]
    b :: Float64 = bs[n]
    wii :: Float64 = self_val == 0. ? 1. : self_val * length(adjs)
    g_p :: Tuple{Float64, Float64} = g_ps[n]
    new_n_val :: Float64 = mine_comp(n_val, adj_vals, wii, b, g_p)
    if new_n_val > 1. || new_n_val < -1.
      error("bad calc in mine step")
    end
    ops_t1[n] = new_n_val
  end
  return ops_t1
end



#####################################################




############## simulate functions ###################
# runs many steps



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

  for _ in 1:n
    ops = degroot_step(ops, selfs, g)
    push!(all_ops, ops)
  end

  return all_ops
end

function beba_sim(n :: Int, g :: Graphs.SimpleGraph{Int}, ops :: Vector{Float64}, selfs :: Vector{Float64}, bs :: Vector{Float64})  :: Vector{Vector{Float64}}
  if !(nv(g) == length(ops) == length(selfs) == length(bs)) 
    error("In Degroot Sim - lengths do not match") 
  end

  all_ops :: Vector{Vector{Float64}} = [ops]

  for _ in 1:n 
    ops = beba_step(ops, selfs, bs, g)
    push!(all_ops, ops)
  end
  return all_ops
end

function mine_sim(n :: Int, g :: Graphs.SimpleGraph{Int}, ops :: Vector{Float64}, selfs :: Vector{Float64}, bs :: Vector{Float64}, g_ps :: Vector{Tuple{Float64, Float64}}) :: Vector{Vector{Float64}}
  if !(nv(g) == length(ops) == length(selfs) == lengths(bs) == length(g_ps))
    error("in mine sim - lengths do not match")
  end
  if maximum(ops) > 1. || maximum([p for (g, p) in g_ps]) > 1. || minimum(bs) < -.1 || minimum([g for (g,p) in g_ps]) < -1.
    error("bad vals in mine sim")
  end
  all_ops :: Vector{Vector{Float64}} = [ops]
  for _ in 1:n 
    ops = mine_step(ops, selfs, bs, g_ps, g)
    push!(all_ops, ops)
  end
  return all_ops
end
