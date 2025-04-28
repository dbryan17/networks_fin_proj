include("graphs.jl")
include("models.jl")
include("runs.jl")
include("dist.jl")

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

# for everything being associative or non associativaite
function mine_assoc(g, n, z) :: Vector{Vector{Float64}}
  g = g 
  n = n 

  ####### ops #######

  # non assoc ops
  ops_dist, ops_dists_lens = make_dists([MyNormal(-.1, .3), MyNormal(.1, .3)], [.5, .5] , nv(g))
  ops_0 :: Vector{Float64} = dists_to_vals(g, ops_dist, ops_dists_lens, -1., 1.)

  # assoc ops
  # ops_0 :: Vector{Float64} = make_assoc_dists([MyNormal(-.9, .3), MyNormal(.9, .3)], [.5, .5], z, -1., 1., 1.)
  ##################

  
  max_B :: Float64 = 1.
  # non assoc bs
  bs_dist, bs_dists_lens = make_dists([MyNumber(0.)], [1.], nv(g))
  bs = dists_to_vals(g, bs_dist, bs_dists_lens, 0., max_B)

  max_G :: Float64 = 1. 
  # non assoc gamma, p
  # gs_dist, gs_dists_lens = make_dists([MyNumber(0.), MyNumber(100.)], [.05, .95], nv(g))
  # ps_dist, ps_dists_lens = make_dists([MyNumber(1.), MyNumber(1.)], [.05, .95], nv(g))
  gs_dist, gs_dists_lens = make_dists([MyNumber(100.)], [1.], nv(g))
  ps_dist, ps_dists_lens = make_dists([MyNumber(1.)], [1.], nv(g))
  if ps_dists_lens != gs_dists_lens 
    error("G and P not matching")
  end
  g_ps_dist = zip(gs_dist, ps_dist)
  g_ps = tup_dists_to_vals(g, g_ps_dist, gs_dists_lens, 0., max_G, 0., 1.)

  # non assoc selfs
  ss_dist, ss_dists_lens = make_dists([MyNumber(0.)], [1.], nv(g))
  selfs = dists_to_vals(g, ss_dist, ss_dists_lens, 0., Inf)

  # this is to change the model to the new version if it is set to true
  true_backfire = false

  all_ops = mine_sim(n, g, ops_0, selfs, bs, g_ps, true_backfire)

  return all_ops
end




# for everything being associative or non associativaite.... same as above except 
function mine_assoc_fix(g, z) :: Vector{Vector{Float64}}
  g = g 
  # n = n 
  z = z

  ####### ops #######

  # non assoc ops
  ops_dist, ops_dists_lens = make_dists([MyNormal(-.3, .3), MyNormal(.3, .3)], [.5, .5] , nv(g))
  ops_0 :: Vector{Float64} = dists_to_vals(g, ops_dist, ops_dists_lens, -1., 1.)

  # assoc ops
  # ops_0 :: Vector{Float64} = make_assoc_dists([MyNormal(-.9, .3), MyNormal(.9, .3)], [.5, .5], z, -1., 1., 1.)
  ##################

  
  max_B :: Float64 = 1.
  # non assoc bs
  bs_dist, bs_dists_lens = make_dists([MyNumber(0.)], [1.], nv(g))
  bs = dists_to_vals(g, bs_dist, bs_dists_lens, 0., max_B)

  max_G :: Float64 = 1. 
  # non assoc gamma, p
  # gs_dist, gs_dists_lens = make_dists([MyNumber(0.), MyNumber(100.)], [.05, .95], nv(g))
  # ps_dist, ps_dists_lens = make_dists([MyNumber(1.), MyNumber(1.)], [.05, .95], nv(g))
  gs_dist, gs_dists_lens = make_dists([MyNumber(1000.)], [1.], nv(g))
  ps_dist, ps_dists_lens = make_dists([MyNumber(1.)], [1.], nv(g))
  if ps_dists_lens != gs_dists_lens 
    error("G and P not matching")
  end
  g_ps_dist = zip(gs_dist, ps_dist)
  g_ps = tup_dists_to_vals(g, g_ps_dist, gs_dists_lens, 0., max_G, 0., 1.)

  # non assoc selfs
  ss_dist, ss_dists_lens = make_dists([MyNumber(0.)], [1.], nv(g))
  selfs = dists_to_vals(g, ss_dist, ss_dists_lens, 0., Inf)

  # this is to change the model to the new version if it is set to true
  true_backfire = false

  n, fin_ops, all_ops = mine_sim_fix(g, ops_0, selfs, bs, g_ps, true_backfire)

  println(n)
  println(fin_ops)

  return all_ops
end

#= 
The plots I will have for my presentation
the basic idea is seeing: 

1) how much these models cause polarization
2) how much these models reduce ploarizaion 

So, I will fix the graph (thinking one for sbm two scale and one for ba scale free)

For each graph model, I will do expirments for each model, there will be two classes of experiments for each model (reversing and creating polarization)

Creating polarization: 

starting opinion distrubtions: 
- uniform
- normal around one point
- bimodal normal
- bimodal uniform

Then, tweak parameters to see what each configuration converges to, and how long it takes to converge to there


Reversing polarization: 

opinion distrubution will be associated with parition

for each starting opinion distrubtion, we will measure it on some level of randomness injected. Fully randomized turns into the expirments for creating polarization (can maybe leverge this for reducing the number of plots)

- 4 unfiorms distrubtions evenly distrubted (must check that c % 4 == 0)
- 2 uniform distrubtions evenly distrubted (c % 2 == 0)
- 2 nomral distrubtions (bimodal, but nomral within each community)

Together, these will just make uniform and bimodal distrubutions

and remember, for each of these, we will inject differewnt levels of randomness

For all of this, the "truth" (mean value, will be 0, will then run other experiments when the truth value is skewed to one side)
=# 


#=
for each model/ graph type

one plot for what parameter level is needed for the given polarization or convergance 
one plot for how long it took 

UUUU I am thining I actually only need one plot. It is parameter level on the x axis, one the 
y axis, it is how long it took to reach converangnce. Then, I color code the boxes based on what type of 
convergance was reached 
if there is two types of convergance for a given value, it gets two boxes. 

for the two parameter models, it will be a heat map. the x and y axis is the two parameters, and the color is convergance type

For each of the two parameter models, we will also have the thing above (time data) for key values for one of the parameters fixed 
like 0 for both, and 1 for both and maybe .5 for both

The ocisinally, we will pull out the varying opinion distrubtion plots

These plots work well for the creating polarization and honestly also for reversing it

For reversing it, another cool plot would be the x and y axis are parameters, then dots in the middle for 
each initial opinion confiugration (could even color code those dots based on graph type). The dots represent the 
average "tipping point" where the convergance went from polarized to truth value 

Need to watch out to cnverging on 0 vs convergening on the truth

OK I think I have all the plots I need for this

aother FUTURE WORK is varying the bias per group. So like people on the negative side are more likely to have bias

I would agruge that the related work isn't true backfire, because it was this forgetting thing - so mine is the first for true backfire

FUTURE WORK associate the biases with different opinioin groups using the same process as for the opinions
this makes sense because it has actually been shown that certain conginative biases (especially confirmation bias and backfire)
are more prominent in conservatives as opposed to liberals

Also the future work of trying this with different ways to model opinions like linear threashold 
also the future work of adding random external opinions like news sources that add in to the update function
can associate certain ones with nodes (conservatives watch this new source)

Research suggests that conservatives are more likely than liberals to resist changing their beliefs in the face of contradictory evidence, particularly when such evidence threatens their worldview. This is not to say that liberals are immune to bias — far from it — but that the psychological mechanisms that underlie belief formation and protection may be more pronounced in conservatives."
it is at the very least true that on certain issues certain sides are more prone to cnfirmation bias and backfire effect 

=# 

# To see how long it takes to converge, get all opinions and lose some level of precesion, then get only
# the unquie values, see how many there are, if it is 2 or 1, then boom! 

# could also do plots that given a fixed parameter, see how much randomness we need to add to cause non polarization to the start
# thats kind of the only interesting thing for randomness. otherwise fix it at like .5

# can also do this for the starting values that are not assoc, how far we have to pull it apart to get stuff

# also have a graph thats how many people need to stop their bias to return to normal

# this is the whole power of one person

# interesting.... for only gamma, the bigger the gamma the further from the middle the converange is test with bi modal



