using Graphs
using GraphPlot
using Plots
using Compose
using LaTeXStrings

# include("graphs.jl")
include("models.jl")
# include("runs.jl")
# include("dist.jl")
include("plots.jl")


function star_graph(n::Int)
  g = SimpleGraph(n)  # n nodes
  for i in 2:n
      add_edge!(g, 1, i)  # connect node 1 to all others
  end
  return g
end



function star_graph_experiment()

  g = star_graph(5)

  # p = gplot(g, layout=spring_layout)
  # draw(SVG("star_graph.svg", 600px, 600px), p) 


  # plot 1 has b = 1 and x1(t) = -1 

  # plot 2 has b = 2.5 and x1(t) = -.5

  # for each plot, we will do it with both versions of mine, and with doing b = 1 g = 0, g = 1 b = 0, b = 1 g = 1

  # for each, we vary x2,3,4,5 from -1 to 1 

  other_x_vals = range(-1., 1., length=402) |> collect 

  x_vals = filter(x -> x != 0., other_x_vals)

  x1_val = -.5

  y_val_mine_full_g1 = []

  # 1, 1 
  for other_x_val in x_vals
    # bs 
    # all 1 
    bs = fill(1., 5)
    g_ps = fill((1.0, 1.0), 5)  
    selfs = fill(0., 5)
    ops = [x1_val, other_x_val, other_x_val, other_x_val, other_x_val]
  
    all_ops = mine_sim(1, g, ops, selfs,bs, g_ps, true )

    push!(y_val_mine_full_g1, all_ops[2][1])
  end

  y_val_mine_full_g5 = []

  # 1, 1 
  for other_x_val in x_vals
    # bs 
    # all 1 
    bs = fill(1., 5)
    g_ps = fill((0.5, 1.0), 5)  
    selfs = fill(0., 5)
    ops = [x1_val, other_x_val, other_x_val, other_x_val, other_x_val]
  
    all_ops = mine_sim(1, g, ops, selfs,bs, g_ps, true )

    push!(y_val_mine_full_g5, all_ops[2][1])
  end

  y_val_mine_g1 = []

  # 1, 1 
  for other_x_val in x_vals
    # bs 
    # all 1 
    bs = fill(1., 5)
    g_ps = fill((1.0, 1.0), 5)  
    selfs = fill(0., 5)
    ops = [x1_val, other_x_val, other_x_val, other_x_val, other_x_val]
  
    all_ops = mine_sim(1, g, ops, selfs,bs, g_ps, false )

    push!(y_val_mine_g1, all_ops[2][1])
  end

  y_val_mine_g5 = []

  # 1, 1 
  for other_x_val in x_vals
    # bs 
    # all 1 
    bs = fill(1., 5)
    g_ps = fill((0.5, 1.0), 5)  
    selfs = fill(0., 5)
    ops = [x1_val, other_x_val, other_x_val, other_x_val, other_x_val]
  
    all_ops = mine_sim(1, g, ops, selfs,bs, g_ps, false )

    push!(y_val_mine_g5, all_ops[2][1])
  end


  y_val_mine_full_10 = []

  # 1, 0 
  for other_x_val in x_vals
    # bs 
    # all 1 
    bs = fill(1., 5)
    g_ps = fill((0.0, 0.0), 5)  
    selfs = fill(0., 5)
    ops = [x1_val, other_x_val, other_x_val, other_x_val, other_x_val]
  
    all_ops = mine_sim(1, g, ops, selfs,bs, g_ps, true )

    push!(y_val_mine_full_10, all_ops[2][1])
  end

  y_val_mine_full_30 = []

  # 1, 0 
  for other_x_val in x_vals
    # bs 
    # all 1 
    bs = fill(3., 5)
    g_ps = fill((0.0, 0.0), 5)  
    selfs = fill(0., 5)
    ops = [x1_val, other_x_val, other_x_val, other_x_val, other_x_val]
  
    all_ops = mine_sim(1, g, ops, selfs,bs, g_ps, true )

    push!(y_val_mine_full_30, all_ops[2][1])
  end


  y_val_beba_1 = []

  # beba
  for other_x_val in x_vals
    # bs 
    # all 1 
    bs = fill(1., 5)
    # g_ps = fill((1.0, 1.0), 5)  
    selfs = fill(0., 5)
    ops = [x1_val, other_x_val, other_x_val, other_x_val, other_x_val]
  
    all_ops = beba_sim(1, g, ops, selfs,bs)

    push!(y_val_beba_1, all_ops[2][1])
  end

  y_val_beba_3 = []

  # beba
  for other_x_val in x_vals
    # bs 
    # all 1 
    bs = fill(3., 5)
    # g_ps = fill((1.0, 1.0), 5)  
    selfs = fill(0., 5)
    ops = [x1_val, other_x_val, other_x_val, other_x_val, other_x_val]
  
    all_ops = beba_sim(1, g, ops, selfs,bs)

    push!(y_val_beba_3, all_ops[2][1])
  end



  labels = [
  "CBB/CBM β = 1, γ = 0",

  "CBB/CBM β = 3, γ = 0",

  "CBB β = 1/0, γ = 1",

  "CBM β = 1/0, γ = 1",

  "BEBA β = 1", 

  "BEBA β = 3", 
  ]

  ys = [y_val_mine_full_10, y_val_mine_full_30, y_val_mine_full_g1, y_val_mine_g1, y_val_beba_1, y_val_beba_3
  ]

  line_graph_mult_y(x_vals, ys, labels, "star graph 1s", "Central Node Update on Start Graph, " * L"y_1 = -0.5")










end

star_graph_experiment()