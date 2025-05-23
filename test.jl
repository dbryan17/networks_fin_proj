using Distributions


include("graphs.jl")
include("plots.jl")
include("models.jl")
include("dist.jl")
include("runs.jl")
include("make_runs.jl")

# err_graph = er_graph(2000, 0.008)
# degree_dist_plot(err_graph, "test1", "er graph")
# wss_graph = ws_graph(2000, 10, 0.5)
# degree_dist_plot(wss_graph, "test2", "ws graph")
# baa_graph = ba_graph(2000, 10)
# degree_dist_plot(baa_graph, "test3", "ba graph")


# TESTING

# e : 0 -> d(1 - 1/c) 
# g = sbm_pp_asoritive(1000, 20., 10, 7.2)
# g = sbm_pp_cp(10000, 100., 2, 4.0, 2.0)

my_heatmap()

error("exit")

# for e... 0 is ER, d(1 - 1/c) is all within com
z, g = sbm_two_scale(2000, 40., 20, 2, 20., 4., 4.)

println((2 * ne(g)) / nv(g) )
degree_dist_plot(g, "test1", "smb assortive two scale dd")

error("exit")
# # all_ops = degroot_full_rand(g, 5, 0.)
# all_ops = degroot_dist(g, 5, 0., Normal(.5, .1))


# all_ops = beba_full_rand(g, 15, 0., 3.1)

all_ops = beba_dists(g, 10, [Uniform(-1., 1.)], [nv(g)], [4.], [nv(g)], 0., [8.], [nv(g)])
ops_dist, ops_dist_len = make_dists([MyUniform(-1. , 1.)], [1.], nv(g))
all_ops = beba_dists(g, 10, ops_dist, ops_dist_len, [4.], [nv(g)], 0., [8.], [nv(g)])

# dists, sizes = n_modal_normal([(.3, .1), (.8, .05)], [.5, .5], nv(g))
# all_ops = degroot_n_dist(g, 15, 2., dists, sizes)

# g = g
# n = 100
# ops_dist, ops_dists_lens = make_dists([MyNormal(0., 1.), MyNormal(0., 1.)], [.5, .5] , nv(g))
# bs_dist, bs_dists_lens = make_dists([MyNumber(1.)], [1.], nv(g))
# max_B :: Float64 = 1.
# gs_dist, gs_dists_lens = make_dists([MyNumber(0.), MyNumber(1.)], [.001, .999], nv(g))
# max_G :: Float64 = 1. 
# ps_dist, ps_dists_lens = make_dists([MyNumber(1.), MyNumber(1.)], [.001, .999], nv(g))
# g_ps_dist = zip(gs_dist, ps_dist)
# if ps_dists_lens != gs_dists_lens 
#   error("G and P not matching")
# end
# ss_dist, ss_dists_lens = make_dists([MyNumber(0.)], [1.], nv(g))
# all_ops = mine_dists(g, n, ops_dist, ops_dists_lens, bs_dist, bs_dists_lens, max_B, g_ps_dist, ps_dists_lens, max_G, ps_dist, ps_dists_lens)

# all_ops = mine_non_assoc(g, 100)
# all_ops = mine_assoc(g, 500, z)
all_ops = mine_assoc_fix(g, z)

op_dist_plot_11(all_ops[1], "test3", "first op")
for op in all_ops[length(all_ops)]
  if op > 1. 
    print("BAD")
  end
end
op_dist_plot_11(all_ops[length(all_ops)], "test2", "last op")


# g1 = sbm_pp_cp(10000, 100., 2, 4.0, 2.0)
# println((2 * ne(g1)) / nv(g1) )

# # e : 0 -> d(1 - 1/c) 
# g2 = sbm_two_scale(10000, 100., 10, 2, 90., 4., 2.)
# println((2 * ne(g2)) / nv(g2) )

# # layout = spring_layout(g2)

# # plot(layout, nodecolor=:lightblue, linewidth=2, marker=:circle, nodelegend=false)
# # savefig("graph_plot.png")

# degree_dist_plot(g2, "testtt", "testttt")