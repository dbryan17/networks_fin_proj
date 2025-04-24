using Plots
using Graphs

# plot degree distrubution as a histogram
function degree_dist_plot(g :: Graphs.SimpleGraph{Int64}, filename :: String, title :: String)
  degrees = [degree(g, v) for v in vertices(g)]
  p = histogram(degrees, bins = 20, title = title, xlabel="Degree", ylabel="Frequency")
  savefig(p, filename * ".pdf")
end


# plot opinion distribution as a historgram
function op_dist_plot(ops :: Vector{Float64}, filename :: String, title :: String) 
  p = histogram(ops, bins = 20, title = title, xlabel ="Opinion", ylabel="Frequency")
  savefig(p, filename * ".pdf")
end