using Distributions


include("graphs.jl")
include("plots.jl")
include("models.jl")

# err_graph = er_graph(2000, 0.008)
# degree_dist_plot(err_graph, "test1", "er graph")
# wss_graph = ws_graph(2000, 10, 0.5)
# degree_dist_plot(wss_graph, "test2", "ws graph")
# baa_graph = ba_graph(2000, 10)
# degree_dist_plot(baa_graph, "test3", "ba graph")


# TESTING

# e : 0 -> d(1 - 1/c) 
g = sbm_pp_asoritive(1000, 20., 10, 7.2)
println((2 * ne(g)) / nv(g) )
degree_dist_plot(g, "test3", "smb assortive dd")

# all_ops = degroot_full_rand(g, 5, 0.)
all_ops = degroot_dist(g, 5, 0., Normal(.5, .1))


op_dist_plot_01(all_ops[1], "test3", "first op")

op_dist_plot_01(all_ops[length(all_ops)], "test2", "last op")




# g1 = sbm_pp_cp(10000, 100., 2, 4.0, 2.0)
# println((2 * ne(g1)) / nv(g1) )

# # e : 0 -> d(1 - 1/c) 
# g2 = sbm_two_scale(10000, 100., 10, 2, 90., 4., 2.)
# println((2 * ne(g2)) / nv(g2) )

# # layout = spring_layout(g2)

# # plot(layout, nodecolor=:lightblue, linewidth=2, marker=:circle, nodelegend=false)
# # savefig("graph_plot.png")

# degree_dist_plot(g2, "testtt", "testttt")