include("isolatedpressure.jl")
include("soluteperm.jl")
include("isolatedmodel.jl")

mean(vec)=sum(vec)/length(vec)
function accumprofile(sol,n_nodes,n_time)
    t_nodes = (60/n_time)*[i for i=1:n_time]
    c_model=zeros(n_time)
    for j = 1:n_time
        c_model[j] = mean(sol[:,j])  #drug concentration at j hours
    end
     t = t_nodes
     return t,c_model
 end

function MSTHighPeclet!(c,p,t)
    #dc = preallocated vector corresponding to f(c,t)= dC/dt
    #c = concentration vector at all spatial points
    #p = vector of parameters as define in modelanalysis
    #t = current time
    #Unpackage parameters
    dc = zeros(length(c))
    P,N,sigma,Perm,Lpt,Svt,Kt,Pvv,Pv,D,r,dr,kd=p
    c0=1 #initial concentration
    tspan=3600 #timespan
    cv=c0*exp(-t*tspan/kd) #time-dependent exponential decay term
    Pe=Lpt*Pv*(Pvv-P[1])*(1-sigma)/Perm #Initial Peclet value
    #interior boundary condition
    dc[1]=tspan*(2*D*(c[2]-c[1])/dr^2 +Lpt*Svt*Pv*(Pvv-P[1])*(1-sigma)*cv)
    for j=2:N-1
        dc[j]=tspan*(((2*D/r[j])*((c[j+1]-c[j])/dr))+(D*(c[j+1]-2*c[j]+c[j-1])/dr^2)+ Kt*((P[j+1]-P[j-1])/(2*dr))*((c[j+1]-c[j])/dr)+ Lpt*Svt*Pv*(Pvv-P[j])*(1-sigma)*cv)
    end
    dc[N]=0 # concentration at the outtermost boundary
    return dc
end

function linstepsizeanalysis(p,n)
    dc = zeros(n)
    P,N,sigma,Perm,Lpt,Svt,Kt,Pvv,Pv,D,r,dr,kd=p
    c0=1 #initial concentration
    tspan=3600 #timespan
    #interior boundary condition
    dc[1]=tspan*(2*D*(-1/dr^2))
    for j=2:N-1
        dc[j]=tspan*((2*D/r[j])*(-1/dr)+(D*(-2)/dr^2)+ Kt*((P[j+1]-P[j-1])/(2*dr))*(-1/dr))
    end
    dc[N]=0 # concentration at the outtermost boundary
    dc = broadcast(abs,dc)
    return lambda
end
function ExplicitEuler(t0,c0,dt,p)
    n_time=Int64(1/dt)
    n_space = length(c0)
    c_out = zeros(n_space,n_time)
    t_out=zeros(n_time,1)
    c_out[:,1]=c0
    t_out[1]=t0
    t=t0
    for j = 2:n_time
        f = MSTHighPeclet!(c_out[:,j-1],p,t)
        c_out[:,j]= c_out[:,j-1] + dt*f
        t=t0+dt
        t_out[j]=t
    end
    return c_out,t_out
end
#n_nodes=2000
#n_nodes = n_time
n_nodes = 50
#N = n spatial
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
r= (1/R)*(range(0,stop=R,length=N))
dr=1/(N-1)
#Get Pressure profile
P = Isolated_Pressure(N,R,Lpt,Svt,Kt,Pvv,at)
#define initial concentration values
c0=zeros(N,1)
c0[N]=0.
#define timespan
#package parameters
p=P,N,sigma,Perm,Lpt,Svt,Kt,Pvv,Pv,D,r,dr,kd
#Define and solve system of ODEs using time constant of 1/n_nodes and the QNDF solver for stiff systems
dt=(1/(n_nodes))
t0=0
c_euler,t_euler = ExplicitEuler(t0,c0,dt,p)
lambda = linstepsizeanalysis(p,100)
h=2/lambda
println("lambda max")
println(lambda)
println("h")
println(h)
#handle EE solution
ntime=n_nodes
c10min= c_euler[:,Int64(round(ntime/6))]
c30min=c_euler[:,ntime÷2]
c60min=c_euler[:,ntime]

#call and handle julia solver solution
sol = Isolated_Model_HP(N,Kt,Lpt,Svt,D,sigma,Perm,R,Pv,Pvv,kd,n_nodes,at)
ntime2=(size(sol.t))[1]
c10min2= sol[:,Int64(round(ntime2/6))]
c30min2=sol[:,ntime2÷2]
c60min2=sol[:,ntime2]

#fetch matlab solution from csv
#c= readlm("c_matlab.csv")
using Plots
#plot
concentrationplot= plot(r,[c10min,c30min,c60min,c10min2,c30min2,c60min2],
xlabel="Dimensionless Distance",
ylabel="Dimensionless Concentration",
title="Concentration vs Radial Postion",
label=["10min-EE" "30min-EE" "1hr-EE" "10min-J" "30min-J" "1hr-J"],
legend=:topleft,
reuse=false
)
display(concentrationplot)

#set up and plot error of low peclet vs high pecelt
c10minerror=broadcast(abs,(c10min .- c10min2)./(c10min2))
c30minerror=broadcast(abs,(c30min .- c30min2)./(c30min2))
c60minerror=broadcast(abs,(c60min .- c60min2)./(c60min2))
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

t,c_model_EE = accumprofile(c_euler,100,n_nodes)
t,c_model_J=accumprofile(sol,100,n_nodes)
accumplot= plot(t,[c_model_EE,c_model_J],
xlabel="Time(min)",
ylabel="Concentration accum",
title="Accumulation vs Time",
label=["Explicit Euler" "Julia Solver"],
legend=:topleft,
reuse=false
)
#savefig(accumplot, "C:\\Users\\sammo\\OneDrive\\Desktop\\Cancer project Julia\\accumulation vs time - normal.png")
display(accumplot)
caccum_error=broadcast(abs,(c_model_J .- c_model_EE)./(c_model_J))
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
