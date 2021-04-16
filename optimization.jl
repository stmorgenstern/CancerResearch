using JuMP,Ipopt,DifferentialEquations,LinearAlgebra
clearconsole()
#import model files
include("isolatedpressure.jl")
include("soluteperm.jl")
include("isolatedmodel.jl")
#Specify time and spatial n values
n_nodes=2000
N=100
#Known metabolic parameters
Kt = 2.9712e-7  ##hydraulic conductivity control value taken from ACS Nano paper
Lpt = 8.6717e-7 #vascular hydraulic conductivity control value taken from ACS nano paper
Svt = 200;      #tumor vascular density
D = 1.375e-07;  #solute diffusion coefficient
rs = 16;      #partical radius (nm)
Perm,sigma=soluteperm(Lpt,rs) #Vascular permeability and solute reflection coefficient respectively
R=1.
Pv=25.
Pvv=1.
r= (range(0,stop=R,length=N))./R
dr=1/(N-1)
kd=1278*60
att=R*sqrt(Lpt*Svt/Kt)
at = [att for i=1:N]
t=(3600/100)*[i for i=1:n_nodes]



co=1
#Call pressure profile
P= Isolated_Pressure(N,R,Lpt,Svt,Kt,Pvv,at)
#Generate Velocity Profile
v=zeros(N,1)
v[2:N-1]=[-.5*(P[i+1]-P[i-1])*dr^-1 for i=2:N-1]
#velocity at the R=1 for isolated tumor
v[N]=-.5*(dr^-1)*(-1*((at[N])^2+P[N-1]*((-1/r[N])*dr^-1 + dr^-2))/((1/r[N])*dr^-1 + dr^-2) - P[N-1])
#package parameters and call explicit euler
c0=zeros(N,1)
c0[N]=0.
#define timespan
#package parameters
p=P,N,sigma,Perm,Lpt,Svt,Kt,Pvv,Pv,D,r,dr,kd
dt=(1/(n_nodes))
t0=0
sol,t_euler = ExplicitEulerMSTHighPeclet!(t0,c0,dt,p)
#Call module to  model for concentration profile
#sol =Isolated_Model_HP(N,Kt,Lpt,Svt,D,sigma,Perm,R,Pv,Pvv,kd,n_nodes,at)

t_euler,c_model = Accumulation_Model(sol,100,n_nodes)


c_data(Peff) = (co*Peff*Svt*kd/(1-Peff*Svt*kd)).*((broadcast(exp,-Peff*Svt.*t_euler)-broadcast(exp,(-1/kd).*t_euler)))
objfunc(Peff)= dot(c_data(Peff)- c_model,c_data(Peff)- c_model)
model = Model(Ipopt.Optimizer)
# @variable(model,0<=x<=1, start=0.5) start can be a vector
@variable(model,0<=Peff <= 1e-3)
register(model, :objfunc, 1, objfunc, autodiff=true)
@NLobjective(model,Min,objfunc(Peff))
JuMP.optimize!(model)
result = value(Peff)
@show result
println("All done")
