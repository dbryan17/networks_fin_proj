using Graphs

include("models.jl")
include("dist.jl")


# helper 
function check_lens(l1, l2)
  if length(l1) != length(l2)
    error("length of sizes does not match length of distributions")
  end
  return true

end

function check_adds(n :: Int, lens_l :: Vector{Vector{Int}})
  if !all( n == reduce(+, lens) for lens in lens_l )
    error("total sizes of distrubtions does not match size of graph")
  end
  return true
end

function ge_0(v :: Float64) 
  if v < 0.
    error("value must be greater than or equal to 0")
  end
  return true
end



################# DEGROOT ##################

# self weights along a distrubtion, ops along a distrubtion
function degroot_dists(g :: Graphs.SimpleGraph{Int}, n :: Int, dist_op, dist_lens_op, dist_s, dist_lens_s) :: Vector{Vector{Float64}}
  check_lens(dist_op, dist_lens_op)
  check_lens(dist_s, dist_lens_s)
  check_adds(nv(g), [dist_lens_op, dist_lens_s])

  ops_0 = dists_to_vals(g, dist_op, dist_lens_op, 0., 1.)
  selfs = dists_to_vals(g, dist_s, dist_lens_s, 0., Inf)

  return degroot_sim(n, g, ops_0, selfs)
end

# self weights fixed, opinions along distributions
function degroot_selfs_fixed(g :: Graphs.SimpleGraph{Int}, n :: Int, s :: Float64, ops_dists, ops_dists_sizes) :: Vector{Vector{Float64}}
  check_lens(dists, dists_sizes)
  check_adds(nv(g), [dists_sizes])
  ge_0(s)

  ops_0 = dists_to_vals(g, ops_dists, ops_dists_sizes, 0., 1.)
  selfs :: Vector{Float64} = fill(s, nv(g))

  return degroot_sim(n, g, ops_0, selfs)
end

# to get just a random number between 0 and 1 give it a uniform distrubution

#################################



########### beba ################

# self weights along disturctions, ops along distrubtions, bs along a distributions
# max b will proably either be 1 or Inf. It's just if we want no chance of backfire effect
function beba_dists(g :: Graphs.SimpleGraph{Int64}, n :: Int, dists_ops, dists_lens_ops, dists_bs, dists_lens_bs, maxB :: Float64, dists_s, dists_lens_s) :: Vector{Vector{Float64}}
  check_lens(dists_ops, dists_lens_ops)
  check_lens(dists_bs, dists_lens_bs)
  check_lens(dists_s, dists_lens_s)
  check_adds(nv(g), [dists_lens_ops, dists_lens_bs, dists_lens_s])

  ops_0 :: Vector{Float64} = dists_to_vals(g, dists_ops, dists_lens_ops, -1., 1.)
  selfs = dists_to_vals(g, dists_s, dists_lens_s, 0., Inf)
  bs = dists_to_vals(g, dists_bs, dists_lens_bs, 0., maxB)

  return beba_sim(n, g, ops_0, selfs, bs)
end

# selfs fixed, ops and bs along distirubtions
function beba_selfs_fixed(g, n :: Int, s :: Float64, op_dists, op_dists_sizes, b_dists, b_dists_sizes, maxB :: Float64) :: Vector{Vector{Float64}}
  check_lens(op_dists_sizes, op_dists)
  check_lens(b_dists, b_dists_sizes)
  check_adds(nv(g), [op_dists_sizes, b_dists_sizes])
  ge_0(s)

  selfs :: Vector{Float64} = fill(s, nv(g))
  ops_0 = dists_to_vals(g, op_dists, op_dists_sizes, -1., 1.)
  bs = dists_to_vals(g, b_dists, b_dists_sizes, 0. , maxB)

  return beba_sim(n, g, ops_0, selfs, bs)
end

function beba_selfs_bs_fixed(g, n :: Int, s :: Float64, b :: Float64, op_dists, op_dists_sizes) :: Vector{Vector{Float64}}
  check_lens(op_dists, op_dists_sizes)
  check_adds(nv(g), [op_dists_sizes])
  ge_(b)

  selfs :: Vector{Float64} = fill(s, nv(g))
  bs :: Vector{Float64} = fill(b, nv(g))
  ops_0 :: Vector{Float64} = dists_to_vals(g, op_dists, op_dists_sizes, -1., 1.)

  return beba_sim(n, g, ops_0, selfs, bs)
 end


################################


############# my model ########

# everything along distrubtions. 
# this is going to be the only one because now I am making it so you can just give the distributions a number to fill it all as the same number



###############################
