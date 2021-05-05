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


# Model for precise adaptation
# N1: N 1
# N2: N 2
function phase_protrait(N_1, N_2)

    #a=(1/(1+10(N_2)^2))^2
    u = (1/(1+10(N_2)^2))^2/(0.1+(1/(1+10(N_2)^2))^2) - N_1

    #b=(1/(1+10(N_1)^2))^2
    v = (1/(1+10(N_1)^2))^2/(0.1+(1/(1+10(N_1)^2))^2) - N_2

    return Point(u,v)
end

# Construct the streamplot
plt1 = Scene(resolution =(1000,1000))
streamplot!(plt1, phase_protrait, 0..1, 0..1,
    gridsize= (32,32), arrow_size = 0.01)

# Display the plot
display(plt1)

# Save the plot
save("Notch_phase_protrait.png", plt1)
