include("graphs.jl")
include("plots.jl")

err_graph = er_graph(2000, 0.008)
degree_dist_plot(err_graph, "test1", "er graph")
wss_graph = ws_graph(2000, 10, 0.5)
degree_dist_plot(wss_graph, "test2", "ws graph")
baa_graph = ba_graph(2000, 10)
degree_dist_plot(baa_graph, "test3", "ba graph")
