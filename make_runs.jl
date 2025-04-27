include("graphs.jl")

# this is a general recipe to create the runs

# this is one where the opinions are not initially aossicatied with groups, just random across the whole graph
# takes only a graph and num iters the idea here is you jsut edit stuff instead this function 
function mine_non_assoc(g, n) :: Vector{Vector{Float64}}
  g = g 
  n = n 
  ops_dist, ops_dists_lens = make_dists([MyNormal(-.9, .3), MyNormal(.9, .3)], [.5, .5] , nv(g))
  bs_dist, bs_dists_lens = make_dists([MyNumber(1.)], [1.], nv(g))
  max_B :: Float64 = 1.
  gs_dist, gs_dists_lens = make_dists([MyNumber(0.), MyNumber(1.)], [.05, .95], nv(g))
  max_G :: Float64 = 1. 
  ps_dist, ps_dists_lens = make_dists([MyNumber(1.), MyNumber(1.)], [.05, .95], nv(g))
  g_ps_dist = zip(gs_dist, ps_dist)
  if ps_dists_lens != gs_dists_lens 
    error("G and P not matching")
  end
  ss_dist, ss_dists_lens = make_dists([MyNumber(0.)], [1.], nv(g))
  # this is to change the model to the new version if it is set to true
  true_backfire = false
  all_ops = mine_dists(g, n, ops_dist, ops_dists_lens, bs_dist, bs_dists_lens, max_B, g_ps_dist, ps_dists_lens, max_G, ss_dist, ss_dists_lens, true_backfire)

  return all_ops


end

# so one



# TODO make one that keeps running until some level of polarization or convergance is reached 


# TODO alter the model so that it doesn't go to zero. Instead of backfire just being kind of mis interpreting 
# someones indeas as for your side to some extend, make it actually be just further than you are on the spectrum only
# it would start as weight 1 for your exact opinion and only grow from there
# this would involved an update function that goes beyond just weights 
