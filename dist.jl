using Distributions
using Graphs

#= 
This file is for making distrubtions
=# 


abstract type MyDistTyp end

struct MyNormal <: MyDistTyp 
  mu :: Float64
  std :: Float64
end 

struct MyUniform <: MyDistTyp
  low :: Float64
  high :: Float64
end

struct MyNumber <: MyDistTyp
  num :: Float64
end 


function convert_to_dist(norm :: MyNormal) 
  return Normal(norm.mu, norm.std)
end
function convert_to_dist(uni :: MyUniform)
  return Uniform(uni.low, uni.high)
end
function convert_to_dist(num :: MyNumber)
  return num.num
end






function make_dists(dists :: Vector{<:MyDistTyp}, fractions :: Vector{Float64}, n :: Int)
  if reduce(+, fractions) != 1.0
    error("make dists - fractions do not add up to 1")
  end
  if length(dists) != length(fractions)
    error("lenghs do not match")
  end

  # takes my dist type and creates the real distrubtion
  dists = [convert_to_dist(d) for d in dists]

  sizes :: Vector{Int} = [floor(n * f) for f in fractions]

  leftover :: Int = n - reduce(+, sizes)

  for i in 1:leftover
    sizes[i] += 1
  end

  if reduce(+, sizes) != n
    error("make dists failing")
  end

  return dists, sizes
  
end

# use the distributions and fractions given to give that fraction of COMMUNITIES (not nodes like the other one)
# that disutrbution

# this is make_dists and dists_to_vals

# this relies on the commmunites being the same length which is true +- 1

# rewrite is the percent of all ops that get rewired (added level of randomness)
function make_assoc_dists(dists :: Vector{<:MyDistTyp}, fractions :: Vector{Float64}, z :: Vector{Int}, minVal :: Float64, maxVal :: Float64, rewire :: Float64) :: Vector{Float64}
  if reduce(+, fractions) != 1.0
    error("make assoc dists - fractions do not add up to 1")
  end

  if length(dists) != length(fractions)
    error("lengths dont match")
  end

  cs :: Vector{Int} = unique(z)

  num_c :: Int = length(cs)

  # calculate number of communites that each dist will have 
  nums :: Vector{Float64} = [f * num_c for f in fractions]

  if !all(x -> isinteger(x), nums)
    error("some fractions given do not work for the number of communities. Bad idxs are $(string(join(findall(x -> !isinteger(x), nums), ", ")))")
  end

  # the index can be the community that has that number
  int_nums :: Vector{Int} = Int.(nums)

  vals = Vector{Float64}()
  offset = 1

  for (i, num_coms) in enumerate(nums)
    for com in 1:num_coms
      com = offset + com - 1
      num_in_com = length( filter(x -> x == com, z) )
      
      dist = convert_to_dist(dists[i])

      if dist isa Number 
        rands = fill(dist, num_in_com)
      else
        rands = rand(dist, num_in_com)
      end

      append!(vals, rands)
    end
    offset += num_coms

  end

  vals = clamp.(vals, minVal, maxVal)

  if length(vals) != length(z)
    error("lengths dont match at end of make assoc dists")
  end


  if rewire != 0. && rewire <= 1.
    all_nodes :: Vector{Int} = collect(1:length(z))
    nodes :: Vector{Int} = collect(1:length(z))
    num_to_rewire :: Int = floor(rewire * length(z))

    for _ in 1:num_to_rewire
      node = rand(nodes)
      node_val = vals[node]
      filter!(n -> n != node, nodes)
      swap_node = rand(all_nodes)
      swap_node_val = vals[swap_node]

      # do the swap 
      vals[node] = swap_node_val
      vals[swap_node] = node_val
    end

  end

  return vals

end


# convert distrubtions to a list of values
function dists_to_vals(g, dists, dists_sizes, minVal :: Float64, maxVal :: Float64) :: Vector{Float64}
  rand_vs = shuffle(vertices(g))
  zipped = zip(dists, dists_sizes)
  vals = Vector{Float64}(undef, nv(g))
  offset = 1
  for (dist, size) in zipped 
    if dist isa Number
      rands = fill(dist, size)
    else
      rands = rand(dist, size)
      rands = clamp.(rands, minVal, maxVal)
    end
    for i in offset:offset + size - 1
      vals[rand_vs[i]] = rands[i - offset + 1]
    end
    offset += size
  end
  return vals
end

# convert distrubtions to a list of values but when there are two values (for the gamma and p tuple)
function tup_dists_to_vals(g, dists_tup, dists_sizes, minVal1 :: Float64, maxVal1 :: Float64, minVal2 :: Float64, maxVal2 :: Float64) :: Vector{Tuple{Float64, Float64}}
  rand_vs = shuffle(vertices(g))
  zipped = zip(dists_tup, dists_sizes)
  vals = Vector{Tuple{Float64, Float64}}(undef, nv(g))
  offset = 1 
  for ((dist1, dist2), size) in zipped
    if dist1 isa Number
      rands1 = fill(dist1, size)
    else
      rands1 = rand(dist1, size)
      rands1 = clamp.(rands1, minVal1, maxVal1) 
    end 

    if dist2 isa Number 
      rands2 = fill(dist2, size)
    else 
      rands2 = rand(dist2, size)
      rands2 = clamp.(rands2, minVal2, maxVal2)
    end 
   
    for i in offset:offset + size - 1
      vals[rand_vs[i]] = (rands1[i - offset + 1], rands2[i - offset + 1])
    end
    offset += size
  end
  return vals
end

