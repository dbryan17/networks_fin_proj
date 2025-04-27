using Plots
using Graphs
using StatsBase
using KernelDensity


# plot degree distrubution as a histogram
function degree_dist_plot(g :: Graphs.SimpleGraph{Int64}, filename :: String, title :: String)
  degrees = [degree(g, v) for v in vertices(g)]
  p = histogram(degrees, bins = 20, title = title, xlabel="Degree", ylabel="Frequency")
  savefig(p, filename * ".png")
end


# plot opinion distribution as a historgram
function op_dist_plot_01(ops :: Vector{Float64}, filename :: String, title :: String) 
  p = histogram(ops, bins = 20, title = title, xlabel ="Opinion", ylabel="Frequency", legend=false)
  xlims!(p, -.1, 1.1)
  savefig(p, filename * ".png")
end


function op_dist_plot(ops :: Vector{Float64}, filename :: String, title :: String) 
  p = histogram(ops, bins = 20, title = title, xlabel ="Opinion", ylabel="Frequency", legend=false)
  savefig(p, filename * ".png")
end


# for -1 to 1
function op_dist_plot_11(ops :: Vector{Float64}, filename :: String, title :: String) 
  p = histogram(ops, bins = 20, title = title, xlabel ="Opinion", ylabel="Frequency", legend=false)
  xlims!(p, -1.1, 1.1)
  savefig(p, filename * ".png")
end


