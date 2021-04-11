clearconsole()
#import model files
include("isolatedpressure.jl")
include("soluteperm.jl")
include("isolatedmodel.jl")
using TypedTables,DelimitedFiles
function testfunction(Lpt,rs)
    n_nodes=2000
    N=100
    #Known metabolic parameters
    Kt = 2.9712e-7  ##hydraulic conductivity control value taken from ACS Nano paper
    Svt = 200;      #tumor vascular density
    D = 1.375e-07;  #solute diffusion coefficient
    Perm,sigma=soluteperm(Lpt,rs) #Vascular permeability and solute reflection coefficient respectively
    R=1.
    Pv=25.
    Pvv=1.
    r= (range(0,stop=R,length=N))./R
    dr=1/(N-1)
    kd=1278*60
    att=R*sqrt(Lpt*Svt/Kt)
    at = [att for i=1:N]
    #Call pressure profile
    P= Isolated_Pressure(N,R,Lpt,Svt,Kt,Pvv,at)
    #Generate Velocity Profile
    v=zeros(N,1)
    v[2:N-1]=[-.5*(P[i+1]-P[i-1])*dr^-1 for i=2:N-1]
    #velocity at the R=1 for isolated tumor
    v[N]=-.5*(dr^-1)*(-1*((at[N])^2+P[N-1]*((-1/r[N])*dr^-1 + dr^-2))/((1/r[N])*dr^-1 + dr^-2) - P[N-1])
    pe=Lpt*Pv*(Pvv-P[end])*(1-sigma)/Perm
    sim_error = false
    sol=Isolated_Model(N,Kt,Lpt,Svt,D,sigma,Perm,R,Pv,Pvv,kd,n_nodes,at)
    try length(sol.t) == 1 ? Error("Numerical Instability in Integration") : NaN
        sim_error=false
    catch e
        sim_error = true
    end
    return pe, sim_error
    #Call module to solve model for concentration profile

end
function gentable()
    lp = range(5e-7,stop=2e-6,length=100)
    ps = range(16,stop=26,length=100)
    LPT=[]
    Particle_Radius=[]
    Peclet=[]
    Simulation_Error=[]
    for pr in ps #iterate over range of radius in range of (5e-7,2e-6) cm
        for lptval in lp #iterate over range of rhydraulic conductivity values in range of (5e-7,2e-6) cm
            pe,se = testfunction(lptval,pr)
            append!(LPT,lptval)
            append!(Particle_Radius,pr)
            append!(Peclet,pe)
            append!(Simulation_Error,se)
        end
    end
    println(" ")
    #println(LPT)
    #println(Particle_Radius)
    #println(Peclet)
    #println(Simulation_Error)
    t = Table(LPT=LPT,Particle_Radius=Particle_Radius,Peclet=Peclet,Simulation_Error=Simulation_Error)
    return t
end
mean(vec)=sum(vec)/length(vec)
function accumprofile(sol,n_time,n_nodes)
    t_nodes = (60/n_time)*[i for i=1:n_time]
    c_model=zeros(n_time)
    for j = 1:n_time
        c_model[j] = mean(sol[:,j*Int64(round(n_nodes/n_time))]);  #drug concentration at j hours
    end
     t = t_nodes
     return t,c_model
 end

#t=gentable()
#using CSV
#CSV.write("output2.csv", t)
#Specify time and spatial n values
n_nodes=2000
N=100
#Known metabolic parameters
Kt = 2.9712e-7  ##hydraulic conductivity control value taken from ACS Nano paper
Lpt =8.6717e-7 #vascular hydraulic conductivity control value taken from ACS nano paper
Svt = 200;      #tumor vascular density
D = 1.375e-07;  #solute diffusion coefficient
rs = 16;      #partical radius (nm)

#highvalues
#rs=20.74
#Lpt=1.85e-6
Perm,sigma=soluteperm(Lpt,rs) #Vascular permeability and solute reflection coefficient respectively
R=1.
Pv=25.
Pvv=1.
r= (range(0,stop=R,length=N))./R
dr=1/(N-1)
kd=1278*60
att=R*sqrt(Lpt*Svt/Kt)
at = [att for i=1:N]
#Call pressure profile
P= Isolated_Pressure(N,R,Lpt,Svt,Kt,Pvv,at)
#Generate Velocity Profile
v=zeros(N,1)
v[2:N-1]=[-.5*(P[i+1]-P[i-1])*dr^-1 for i=2:N-1]
#velocity at the R=1 for isolated tumor
v[N]=-.5*(dr^-1)*(-1*((at[N])^2+P[N-1]*((-1/r[N])*dr^-1 + dr^-2))/((1/r[N])*dr^-1 + dr^-2) - P[N-1])

#Call module to solve model for concentration profile
sol =Isolated_Model(N,Kt,Lpt,Svt,D,sigma,Perm,R,Pv,Pvv,kd,n_nodes,at)
sol2=Isolated_Model_HP(N,Kt,Lpt,Svt,D,sigma,Perm,R,Pv,Pvv,kd,n_nodes,at)
#Plot the data
t,c_model_lowp = accumprofile(sol,100,n_nodes)
t,c_model_highp=accumprofile(sol2,100,n_nodes)

using Plots
#plot pressure profile

#plot concentration profile
ntime=(size(sol.t))[1]
c10min= sol[:,Int64(round(ntime/6))]
c30min=sol[:,ntime÷2]
c60min=sol[:,ntime]
c10min2= sol2[:,Int64(round(ntime/6))]
c30min2=sol2[:,ntime÷2]
c60min2=sol2[:,ntime]
concentrationplot= plot(r,[c10min,c30min,c60min,c10min2,c30min2,c60min2],
xlabel="Dimensionless Distance",
ylabel="Dimensionless Concentration",
title="Concentration vs Radial Postion",
label=["10min-LP" "30min-LP" "1hr-LP" "10min-HP" "30min-HP" "1hr-HP"],
legend=:topleft,
reuse=false,
xlims=(0,.7),
ylims=(0,.04)
)
#savefig(concentrationplot, "LP HP concentration vs position - high particle size")
display(concentrationplot)
#set up and plot error of low peclet vs high pecelt
c10minerror=broadcast(abs,(c10min .- c10min2)./(c10min))
c30minerror=broadcast(abs,(c30min .- c30min2)./(c30min))
c60minerror=broadcast(abs,(c60min .- c60min2)./(c60min))
errorplot= plot(r,[c10minerror,c30minerror,c60minerror],
xlabel="Dimensionless Distance",
ylabel="ϵ",
title="ϵ vs Radial Postion",
label=["10min-error" "30min-error" "1hr-error" ],
legend=:topleft,
reuse=false
)
#savefig(errorplot, "error vs position - high particle size")
display(errorplot)

accumplot= plot(t,[c_model_lowp,c_model_highp],
xlabel="Time(min)",
ylabel="Concentration accum",
title="Accumulation vs Time",
label=["Low-P" "High-P"],
legend=:topleft,
reuse=false
)
#savefig(accumplot, "C:\\Users\\sammo\\OneDrive\\Desktop\\Cancer project Julia\\accumulation vs time - normal.png")
display(accumplot)
caccum_error=broadcast(abs,(c_model_lowp .- c_model_highp)./(c_model_lowp))
accumerrorplot= plot(t,caccum_error,
xlabel="Time(min)",
ylabel="ϵ",
title="Accumulation Error vs Time",
label=["ϵ"],
legend=:topleft,
reuse=false
)
#savefig(accumerrorplot, "C:\\Users\\sammo\\OneDrive\\Desktop\\Cancer project Julia\\accumulation error vs time - normal.png")
display(accumerrorplot)
