using Distributions

#= 
This file is for making distrubtions
=# 

# fractions is fraction of all that each distrubtion should be
# n is the total number to make
function n_modal_normal(mu_stds :: Vector{Tuple{Float64, Float64}}, fractions :: Vector{Float64}, n :: Int)
  if reduce(+, fractions) != 1.0
    error("n modal - fractions do not add up to 1")
  end

  dists = [Normal(mu, std) for (mu, std) in mu_stds]

  sizes :: Vector{Int} = [floor(n * f) for f in fractions]

  leftover :: Int = n - reduce(+, sizes)

  for i in 1:leftover
    sizes[i] += 1
  end
  return dists, sizes

end