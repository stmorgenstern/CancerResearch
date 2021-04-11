clearconsole()
#import model files
include("isolatedvelocity.jl")
include("isolatedpressure.jl")
include("soluteperm.jl")
include("isolatedmodel.jl")
include("makeplots.jl")
#Specify time and spatial n values
n_time = 2000
n_spatial = 100
#Known metabolic parameters
Kt = 2.9712e-7  ##hydraulic conductivity control value taken from ACS Nano paper
Lpt = 8.6717e-7 #vascular hydraulic conductivity control value taken from ACS nano paper
Svt = 200;      #tumor vascular density
D = 1.375e-07;  #solute diffusion coefficient
rs = 16;      #partical radius (nm)
Perm,sigma=soluteperm(Lpt,rs) #Vascular permeability and solute reflection coefficient respectively
R=1. # Tumor Radius
Pv=25. #
Pvv=1.
#Define Spatial domain
r= (range(0,stop=R,length=n_spatial))./R
dr=1/(n_spatial-1)

kd=1278*60 # blood circulation time of drug in hours;
att=R*sqrt(Lpt*Svt/Kt) #parameter alpha for tumor, ignoring lymphatics
at = [att for i=1:n_spatial]

#Call pressure profile
P= Isolated_Pressure(n_spatial,R,Lpt,Svt,Kt,Pvv,at)
#Generate Velocity Profile
v =Isolated_Velocity(n_spatial,P,r, dr)
#Call module to solve model for concentration profile
sol =Isolated_Model(n_spatial,Kt,Lpt,Svt,D,sigma,Perm,R,Pv,Pvv,kd,n_time,at)

#Plot the data
using Plots
#plot pressure profile
pplot = pressureplot(r,P)
display(pplot)
#plot velocity profile
vplot = velocityplot(r,v)
display(vplot)
#plot concentration profile at 3 time nodes
cplot = singleconcplot(r,sol)
display(cplot)

#plot spatially integrated accumulation
aplot = singleaccumplot(sol)
display(aplot)
