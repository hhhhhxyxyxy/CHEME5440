# Example for use of DifferentialEquations.jl
# CHEME5440/7770 - Spring 2021
# To install the package, use the following commands inside the Julia REPL:
# using Pkg
# Pkg.add("DifferentialEquations")
# For this example, you will also need Plots; To add:
# Pkg.add("Plots")


using DifferentialEquations     # Include DifferentialEquations.jl
using Plots                     # Include Plots.jl for plotting
gr(show = true)  # Use the gr backend for plotting and show plots


# Model parameters

l=0.0001
a_b = 1
d_b = 1
k_b = 0
a_bp = 0.1
d_bp = 0.01
k_bp = 1
k_plus = 1
k_minus = 1
alpha_1_minus = l/(1+l)
alpha_1_plus = l/(1+l)
beta_1 = 2.5l/(1+l)
E_1 = 0.01  #nM
B = 0.002   #uM




# du: Diffrerential equations
# u: Time-dependent variables, [E_1_star, B, B_p, E_1_star_mult_B, E_1_star_mult_B_p]
# p: Additional model parameters
# t: time
# Note "!" point after function name is a Julia convention that indicates
# that the function will modify values in one or more of the function arguments.  In this case,
# twostate will modify values in the input vector, du
function twostate!(du,u,p,t)
 du[1] = alpha_1_plus*E_1 -alpha_1_minus*u[1] - a_b * u[1]*u[2]-a_bp*u[1]*u[3]+d_b*u[4]+d_bp*u[5]
 du[2] = d_b*u[4]+k_b*u[4]+beta_1*u[4]+k_minus*u[3]-a_b*u[1]*u[2]-a_bp*u[1]*u[3]-k_plus*u[2]
 du[3] = d_bp*u[5]+k_bp*u[5]+beta_1*u[5]-a_bp*u[1]*u[3]
 du[4] = a_b*u[1]*u[2]-d_b*u[4]-beta_1*u[1]*u[2]
 du[5] = a_bp*u[1]*u[3]-d_bp*u[5]-beta_1*u[1]*u[3]
end

#get the initial conditions, ie their steady states
function estimate_steady_state()
  steady_state_prob = SteadyStateProblem(twostate!,u0)
  steady_state_soln  = solve(steady_state_prob,SSRootfind())
  return abs.(steady_state_soln)
  end


u0 = estimate_steady_state()
u0[1] = B                               #initial conc of B = 0.1uM
tspan = (0.0,1000.0)                     #time interval (start time, end time)
prob = ODEProblem(twostate!,u0,tspan)     #Create an ODE problem for the twostate fxn
sol = solve(prob)                       #Solve the system


#Plot the results; Output vs time
plt1 = plot(sol, xaxis="time", yaxis = "Output",legend=false,title="l=0.0001")
display(plt1)



#Save the three plots as PNG files
savefig(plt1, "./l0.0001.png")
