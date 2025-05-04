# this is for making the plots for how much we need to crank a parameter to reach polarization for many diffeferent graph configurations
include("graphs.jl")
include("plots.jl")
include("models.jl")
include("dist.jl")
include("runs.jl")
include("make_runs.jl")


using StatsBase


### THE first one was for 



# for everything being associative or non associativaite
function mine_assoc(g, n, ops_0) :: Vector{Vector{Float64}}
  g = g 
  n = n 

  ####### ops #######

  # non assoc ops
  # ops_dist, ops_dists_lens = make_dists([MyNormal(-.1, .3), MyNormal(.1, .3)], [.5, .5] , nv(g))
  # ops_0 :: Vector{Float64} = dists_to_vals(g, ops_dist, ops_dists_lens, -1., 1.)

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
  gs_dist, gs_dists_lens = make_dists([MyNumber(1.)], [1.], nv(g))
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

function mine_assoc_fix(g, ops_0)
  g = g 
  # n = n 

  ####### ops #######



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
  gs_dist, gs_dists_lens = make_dists([MyNumber(1.), MyNumber(0.)], [1., 0.], nv(g))
  ps_dist, ps_dists_lens = make_dists([MyNumber(1.), MyNumber(1.)], [1., 0.], nv(g))
  if ps_dists_lens != gs_dists_lens 
    error("G and P not matching")
  end
  g_ps_dist = zip(gs_dist, ps_dist)
  g_ps = tup_dists_to_vals(g, g_ps_dist, gs_dists_lens, 0., max_G, 0., 1.)

  # non assoc selfs
  ss_dist, ss_dists_lens = make_dists([MyNumber(0.)], [1.], nv(g))
  selfs = dists_to_vals(g, ss_dist, ss_dists_lens, 0., Inf)

  # this is to change the model to the new version if it is set to true
  true_backfire = true

  always_ops_0 = deepcopy(ops_0)

  # b_vals = range(xmin, xmax; length=nx) |> collect
  # yvals = range(ymin, ymax; length=ny) |> collect

  n, fin_val, all_ops = mine_sim_fix(g, always_ops_0, selfs, bs, g_ps, true_backfire)
  # println(n)
  println(fin_val)
  # error("quite")
  # all_ops = mine_sim(1000, g, ops_0, selfs, bs, g_ps, true_backfire)

  # println(n)
  # println(fin_val)

  # println(mean(all_ops[length(all_ops)]))
  # println(mean(all_ops[1]))

  # op_heatmap("test")

  return fin_val, all_ops[length(all_ops)]


end

function diff_index(ops) 
  total = sum(abs(ops[i] - ops[j]) for i in eachindex(ops), j in eachindex(ops))
  return total
end

function make_plot(g)
  # non assoc ops


  all_init_ops = []

  ops_dist, ops_dists_lens = make_dists([MyNormal(-.9, .2), MyNormal(.1, .2)], [.5, .5] , nv(g))
  ops_0  = dists_to_vals(g, ops_dist, ops_dists_lens, -1., 1.)
  push!(all_init_ops, ops_0)
  ops_dist, ops_dists_lens = make_dists([MyNormal(-.7, .2), MyNormal(.2, .2)], [.5, .5] , nv(g))
  ops_0 = dists_to_vals(g, ops_dist, ops_dists_lens, -1., 1.)
  push!(all_init_ops, ops_0)
  ops_dist, ops_dists_lens = make_dists([MyNormal(-.9, .2), MyNormal(.8, .2)], [.5, .5] , nv(g))
  ops_0  = dists_to_vals(g, ops_dist, ops_dists_lens, -1., 1.)
  push!(all_init_ops, ops_0)
  ops_dist, ops_dists_lens = make_dists([MyNormal(-.2, .2), MyNormal(.2, .2)], [.5, .5] , nv(g))
  ops_0 = dists_to_vals(g, ops_dist, ops_dists_lens, -1., 1.)
  push!(all_init_ops, ops_0)
  ops_dist, ops_dists_lens = make_dists([MyNormal(-.3, .2), MyNormal(.6, .2)], [.5, .5] , nv(g))
  ops_0  = dists_to_vals(g, ops_dist, ops_dists_lens, -1., 1.)
  push!(all_init_ops, ops_0)
  ops_dist, ops_dists_lens = make_dists([MyNormal(-.3, .2)], [1.] , nv(g))
  ops_0  = dists_to_vals(g, ops_dist, ops_dists_lens, -1., 1.)
  push!(all_init_ops, ops_0)
  ops_dist, ops_dists_lens = make_dists([MyNormal(0, .2)], [1.] , nv(g))
  ops_0 = dists_to_vals(g, ops_dist, ops_dists_lens, -1., 1.)
  push!(all_init_ops, ops_0)
  ops_dist, ops_dists_lens = make_dists([MyNormal(5, .2)], [1.] , nv(g))
  ops_0 = dists_to_vals(g, ops_dist, ops_dists_lens, -1., 1.)
  push!(all_init_ops, ops_0)
  ops_dist, ops_dists_lens = make_dists([MyUniform(-1., 1.)], [1.] , nv(g))
  ops_0 = dists_to_vals(g, ops_dist, ops_dists_lens, -1., 1.)
  push!(all_init_ops, ops_0)
  ops_dist, ops_dists_lens = make_dists([MyNumber(-1.), MyNumber(1.)], [.5, .5] , nv(g))
  ops_0= dists_to_vals(g, ops_dist, ops_dists_lens, -1., 1.)
  push!(all_init_ops, ops_0)
  ops_dist, ops_dists_lens = make_dists([MyNumber(-.5), MyNumber(.5)], [.5, .5] , nv(g))
  ops_0 = dists_to_vals(g, ops_dist, ops_dists_lens, -1., 1.)
  push!(all_init_ops, ops_0)
  ops_dist, ops_dists_lens = make_dists([MyNumber(-.2), MyNumber(.7)], [.5, .5] , nv(g))
  ops_0 = dists_to_vals(g, ops_dist, ops_dists_lens, -1., 1.)
  push!(all_init_ops, ops_0)


  d_i_start_is = []
  d_end_is = []
  init_means = []
  end_vals = []
  end_fracs = []
  fin_mean = []
  for op in all_init_ops
    d_i = diff_index(op) / nv(g)^2
    push!(d_i_start_is, d_i)
    push!(init_means, sum(op) / length(op))
    inner_fin_vals = []
    inner_d_ends = []
    inner_frac = []
    inner_fin_mean = []
    always_0 = deepcopy(op)
    for _ in 1:3
      always_0 = deepcopy(op)

      print("here")



      # fin_vals, end_ops = mine_assoc_fix(g, always_0)

      all_ops = mine_assoc_fix = mine_assoc(g, 5000, always_0)



      ###### delte this 
      # op_dist_plot_11_avg(all_ops[length(all_ops)], "end", "Opinions after 5000 iterations")
      # op_dist_plot_11_avg(all_ops[1], "start", "Starting Opinion Values")
      epsilon = 0.01
      end_oppp = deepcopy(all_ops[length(all_ops)])
      roundeds = round.(end_oppp ./ epsilon) .* epsilon
      p = abs(mode(roundeds))
      push!(inner_fin_vals, p)
      d_i_end = diff_index(all_ops[length(all_ops)]) / nv(g)^2
      push!(inner_d_ends, d_i_end)
      end_ops_abs = [abs(o) for o in roundeds]
      all_p = filter(o -> o == p, end_ops_abs)
      all_p_frac = length(all_p) / nv(g)
      push!(inner_frac, all_p_frac)
      push!(inner_fin_mean, sum(all_ops[length(all_ops)]) / nv(g))

      

      # error("end")
      ########
      # p = abs(fin_vals[1])
      # push!(inner_fin_vals, p)
      # d_i_end = diff_index(end_ops) / nv(g)
      # push!(inner_d_ends, d_i_end)
    end
    push!(d_end_is, sum(inner_d_ends) / length(inner_d_ends))
    push!(end_vals, sum(inner_fin_vals) / length(inner_fin_vals))
    push!(end_fracs, sum(inner_frac) / length(inner_frac))
    push!(fin_mean, sum(inner_fin_mean) / length(inner_fin_mean))
  end


  return d_i_start_is, d_end_is, init_means, end_vals, end_fracs, fin_mean



end


gs = []
z, g = sbm_two_scale(2000, 40., 20, 2, 20., 4., 4.)
push!(gs, g)


g = er_graph(2000, 0.002)
push!(gs, g)

g = er_graph(2000, .01)
push!(gs, g)
# make_plot(g)


g = ws_graph(2000, 20, .001)
push!(gs, g)

g = ws_graph(2000, 20, .1)
push!(gs, g)

g = ws_graph(2000, 20, .4)
push!(gs, g)


g = ws_graph(2000, 40, .001)
push!(gs, g)

g = ws_graph(2000, 40, .1)
push!(gs, g)


g = ws_graph(2000, 40, .4)
push!(gs, g)


g = ba_graph(2000, 40)
push!(gs, g)

g = ba_graph(2000, 20)
push!(gs, g)

g = ba_graph(200, 10)
push!(gs, g)



g = sbm_pp_asoritive(2000, 40., 8, 15.)
push!(gs, g)



g = sbm_pp_asoritive(2000, 40., 8, 30.)
push!(gs, g)


g = sbm_pp_asoritive(2000, 40., 30, 30.)
push!(gs, g)


g = sbm_pp_asoritive(2000, 40., 30, 15.)
push!(gs, g)


dss = []
des = []
imeans = []
endvals = []
endfracs = []
endmeans = []
i = 1
for g in gs
  global i
  d_i_start_is, d_end_is, init_means, end_vals, fracs, end_mean = make_plot(g)
  println(i)
  i += 1
  push!(dss, d_i_start_is)
  push!(des, d_end_is)
  push!(imeans, init_means)
  push!(endvals, end_vals)
  push!(endfracs, fracs)
  push!(endmeans, end_mean)
end

# start -> end
scatter_plot_l(dss, des, "0-ds-de", "Disagreement Index Relationship", "Starting Disagreement", "End Disagreement")
scatter_plot_l(dss, endvals, "0-ds-ev", "Disagreement Index vs. |Mode|", "Starting Disagreement", "Absolute Value of Mode")
scatter_plot_l(dss, endfracs, "0-ds-fracs", "Disagreement Index vs. Polarization Fraction", "Starting Disagreement", "Fraction of Nodes with |mode|")
scatter_plot_l(dss, endmeans, "0-ds-em", "Disagreement Index vs. Final Mean", "Starting Disagreement", "Final Mean")

scatter_plot_l(imeans, des, "0-im-des", "Initial Mean vs. Final Disagreement", "Initial Mean", "Final Disagreement")
scatter_plot_l(imeans, endvals, "0-im-ev", "Initial Mean vs. |Mode|", "Initial Mean", "Absolute Value of Mode")
scatter_plot_l(imeans, endfracs, "0-im-fracs", "Initial Mean vs. Polarization Fraction", "Initial Mean", "Fraction of Nodes with |mode|")
scatter_plot_l(imeans, endmeans, "0-im-em", "Initial Mean vs. Final Mean", "Initial Mean", "Final Mean")

scatter_plot_l(endvals, endfracs, "0-ev-fracs", "Final |Mode| vs. Polarization Fraction", "|mode|", "Fraction of Nodes with |mode|")
scatter_plot_l(endmeans, endvals, "0-em-ev", "Final mean vs. Final |Mode|", "Final Mean", "Final |mode|")
scatter_plot_l(endmeans, endfracs, "0-em-fracs", "Final mean vs. Polarization Fraction", "Final mean", "Fraction of Nodes with |mode|")






