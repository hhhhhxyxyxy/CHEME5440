# Example for constructing a streamplot in Julia for phase portraits
# CHEME5440/7770 - Spring 2021
# The example makes use of Makie.jl and AbstractPlotting.jl for plotting
# To install the packages, use the following commands inside the Julia REPL:
# using Pkg
# Pkg.add("Makie")
# Pkg.add("AbstractPlotting")

using Makie
using AbstractPlotting
using AbstractPlotting.MakieLayout
AbstractPlotting.inline!(true)
using DifferentialEquations
using Plots                     # Include Plots.jl for plotting
gr(show = true)  # Use the gr backend for plotting and show plots



# Model for precise adaptation
# D1: D 1
# D2: D 2
function phase_protrait(D_1, D_2)

    v = 0.00001     #gammaD/ gammaN, inhibition rate constant
    a =  D_1^2/(0.1+D_1^2) #f(D_1)
    b = D_2^2/(0.1+D_2^2)  #f(D_2)
    c = 1/(1+10*a^2)  #g(f(D_1))
    d = 1/(1+10*b^2)  #g(f(D_2))

    u = (d-D_1)*v #dD_1/dtau
    w = (c-D_2)*v #dD_2/dtau


    return Point(u,w)
end

function estimate_steady_state()
  steady_state_prob = SteadyStateProblem(phase_protrait(10,10),0)
  steady_state_soln  = solve(steady_state_prob,SSRootfind())
  return abs.(steady_state_soln)
  end
XSS = estimate_steady_state()

plot(XSS[:,1],XSS[:,2], legend=false)


# Construct the streamplot
plt1 = Scene(resolution =(1000,1000))
streamplot!(plt1, phase_protrait, 0..1, 0..1, colormap = :plasma,
    gridsize= (32,32), arrow_size = 0.01)

# Display the plot
display(plt1)

# Save the plot
save("phase_protrait.png", plt1)
