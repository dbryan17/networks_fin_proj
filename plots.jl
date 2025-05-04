using Plots
using Graphs
using StatsBase
using KernelDensity
using LaTeXStrings


# plot degree distrubution as a histogram
function degree_dist_plot(g :: Graphs.SimpleGraph{Int64}, filename :: String, title :: String)
  degrees = [degree(g, v) for v in vertices(g)]
  p = histogram(degrees, bins = 20, title = title, xlabel="Degree", ylabel="Frequency")
  savefig(p, filename * ".png")
end

function degree_dist_plot_w_avg(g :: Graphs.SimpleGraph{Int64}, avg :: Float64, filename :: String, title :: String)
  degrees = [degree(g, v) for v in vertices(g)]
  p = histogram(degrees, bins = 20, title = title, xlabel="Degree", ylabel="Frequency", label = "", legend = true)
  vline!([avg], color=:red, label="average degree")
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


function line_graph_mult_y(xs, ys :: Vector{Vector{Any}}, labels :: Vector{String}, filename :: String, title)
  p = plot(xs, ys[1], label=labels[1], lw=2, title = title, xlabel = L"y_{2,3,4}(t)", ylabel = L"y_1(t+1)")
  for i in 2:length(ys) 
    print(labels[i])
    plot!(xs, ys[i], label=labels[i], lw=2)
  end
  savefig(p, filename * ".png")
end

function op_dist_plot_11_avg(ops :: Vector{Float64}, filename :: String, title :: String) 
  avg = mean(ops)
  p = histogram(ops, bins = 20, title = title, xlabel ="Opinion", ylabel="Frequency", legend=true, label ="")
  vline!([avg], color=:red, label="average value")
  xlims!(p, -1.1, 1.1)
  savefig(p, filename * ".png")
end

function get_color(val)
  x, y = val
  diff = abs(x - y)  # Calculate the absolute difference between the two values
  
  # Map the difference (0 -> green, 2 -> red)
  g = max(0.0, 1.0 - diff / 2.0)  # Green decreases as difference increases
  r = min(1.0, diff / 2.0)       # Red increases as difference increases
  
  return RGB(r, g, 0.0)           # No blue component
end


function scatter_plot(xs, ys, filename, title, xlab, ylab)
  p = scatter(xs, ys, title = title, xlabel = xlab, ylabel = ylab)
  savefig(p, filename * ".png")
end

function scatter_plot_l(xss, yss, filename, title, xlab, ylab)
  p = scatter([], [], title = title, xlabel = xlab, ylabel = ylab, legend=false)
  for (xs, ys) in zip(xss, yss)
    scatter!(p, xs, ys, color=:auto)
  end
  savefig(p, filename * ".png")
end



# function op_heatmap(x_vals :: Vector{Float64}, y_vals :: Vector{Float64}, converged_to, time_to, filename :: String)
#   # z = [i * j for i in 1:10, j in 1:10]
#   # p = heatmap(z)
#   rounded_converged_to = map(v -> 0.1 .* round.(v ./ 0.1), converged_to)
#   normalized_converged_to = map(v -> map(x -> x == -0.0 ? 0.0 : x, v), rounded_converged_to)
  
#   poss_vals = unique(normalized_converged_to)

#   pair_colors = Dict()
#   for (i, pair) in enumerate(poss_vals)
#     pair_colors[tuple(pair)] = abs(pair[1] - pair[2])
#   end
  


#   # Normalize times to [0, 1] based on global max
#   max_time = maximum(time_to)
#   normalized_time = time_to ./ max_time


#   heatmap_data = zeros(length(x_vals), length(y_vals))
#   for i in 1:length(x_vals)
#       for j in 1:length(y_vals)
#           idx = (i - 1) * length(y_vals) + j
#           pair_val = normalized_converged_to[idx]
#           time_val = normalized_time[idx]
#           color_value = pair_colors[tuple(pair_val)] * (1 - time_val)
#           heatmap_data[i, j] = color_value
#       end
#   end

#   # Step 6: Plot heatmap with color scaling based on `heatmap_data`
#   p = plot()
#   heatmap!(heatmap_data, color=:auto, xlabel="X", ylabel="Y", title="Convergence Heatmap")
#   savefig(p, filename * ".png")
# end



function op_heatmap(x_vals :: Vector{Float64}, y_vals :: Vector{Float64}, converged_to, time_to, filename :: String)
  # z = [i * j for i in 1:10, j in 1:10]
  # p = heatmap(z)
  rounded_converged_to = map(v -> 0.05 .* round.(v ./ 0.05), converged_to)
  normalized_converged_to = map(v -> map(x -> x == -0.0 ? 0.0 : x, v), rounded_converged_to)
  
  poss_vals = unique(normalized_converged_to)

  println(poss_vals)


  # Normalize times to [0, 1] based on global max
  max_time = maximum(time_to)
  normalized_time = 1 .- (time_to ./ max_time)

  # scaled_time = 0.5 .+ 0.5 .* normalized_time
  scaled_time = 0.5 .+ 0.5 .* (normalized_time .- 0.5)

  heatmap_data = Matrix{RGB{Float64}}(undef, length(x_vals), length(y_vals))

  for i in 1:length(x_vals)
      for j in 1:length(y_vals)
          idx = (i - 1) * length(y_vals) + j
          pair_val = normalized_converged_to[idx]
          time_val = scaled_time[idx]
          base_color = RGB(1,1,1)  # default white
          if length(pair_val) == 1
            pair_val = [pair_val[1], pair_val[1]]
          end
          if length(pair_val) == 2
            one = maximum(pair_val)
            two = minimum(pair_val)
            # 0
            if (isapprox(one, 0.0; atol=1e-3) ||isapprox(one, 0.05; atol=1e-3) ) && (isapprox(two, 0.0; atol=1e-3) || isapprox(two, -0.05; atol=1e-3))
                base_color = RGB(0.0, 0.0, 1.0)  # blue

            # -.2
            elseif (isapprox(one, -0.2; atol=1e-3) || isapprox(one, -0.25; atol=1e-3))  && (isapprox(two, -0.2; atol=1e-3)|| isapprox(two, -0.25; atol=1e-3))
              base_color = RGB(0,1,0) # green
            # .5, -.5
            elseif (isapprox(one, .5; atol=1e-3) || isapprox(one, .55; atol=1e-3) || isapprox(one, .45; atol=1e-3)) && (isapprox(two, -0.45; atol=1e-3) || isapprox(two, -0.5; atol=1e-3) || isapprox(two, -0.55; atol=1e-3))
                base_color = RGB(1,0,0)  # red

            # .3, -.65
            elseif (isapprox(one, .3; atol=1e-3) || isapprox(one, .35; atol=1e-3)) && (isapprox(two, -0.65; atol=1e-3) || isapprox(two, -.60; atol=1e-3))
                base_color = RGB(1, 1, 0) # yellow
            else 
              println(pair_val)
            end
          else 
            println("converged to 3 or more numbers") 
          end

          # Blend with white based on time (faster = stronger color)
          faded_color = RGB(
            base_color.r * time_val,
            base_color.g  * time_val,
            base_color.b *  time_val
          )
          heatmap_data[i, j] = faded_color
      end
  end

  x_range = maximum(x_vals) - minimum(x_vals)
  y_range = maximum(y_vals) - minimum(y_vals)
  println(minimum(y_vals))
  println("max time")
  println(max_time)
  aspect_ratio_value = x_range / y_range
  # Step 6: Plot heatmap with color scaling based on `heatmap_data`
  p = plot(     
          aspect_ratio=aspect_ratio_value,

         xlims=(minimum(x_vals), maximum(x_vals)),
          ylims=(minimum(y_vals), maximum(y_vals))
  )
          heatmap!(p, 
          x_vals, 
          y_vals,
          heatmap_data, xlabel="β value", ylabel="fraction of nodes with γ = 1, p = 1", 
          title="Convergence Heatmap", 
          yflip = false,
          
          # ratio=1,
          # size=(600, 600),
          colorbar=false,
          )

  savefig(p, filename * ".png")
end

