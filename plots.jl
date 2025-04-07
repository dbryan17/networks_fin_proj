using Plots
using Graphs

function degree_dist_plot(g :: Graphs.SimpleGraph{Int64}, filename :: String, title :: String)
  degrees = [degree(g, v) for v in vertices(g)]
  p = histogram(degrees, bins = 20, title = title, xlabel="Degree", ylabel="Frequency")
  savefig(p, filename * ".pdf")
end