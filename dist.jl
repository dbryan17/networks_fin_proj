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

