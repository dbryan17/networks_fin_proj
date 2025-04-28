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

function mine_assoc(g, n, )

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


=# 

# To see how long it takes to converge, get all opinions and lose some level of precesion, then get only
# the unquie values, see how many there are, if it is 2 or 1, then boom! 

# TODO make one that keeps running until some level of polarization or convergance is reached 


