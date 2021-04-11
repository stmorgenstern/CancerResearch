include("isolatedpressure.jl")
include("soluteperm.jl")
include("isolatedmodel.jl")
include("makeplots.jl")
n_nodes=2000
#n_nodes = n_time
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

using DelimitedFiles
#fetch matlab solution from csv
c_mat= readdlm("c_matlab_hp.csv",',')
c_mat = transpose(c_mat)
print(size(c_mat))

#handle EE solution
ntime=n_nodes
c10min= c_mat[:,Int64(round(ntime/6))]
c30min=c_mat[:,ntime÷2]
c60min=c_mat[:,ntime]

#call and handle julia solver solution
sol = Isolated_Model_HP(N,Kt,Lpt,Svt,D,sigma,Perm,R,Pv,Pvv,kd,n_nodes,at)
println(size(sol))

ntime2=(size(sol.t))[1]
c10min2= sol[:,Int64(round(ntime2/6))]
c30min2=sol[:,ntime2÷2]
c60min2=sol[:,ntime2]

concentrationplot=doubleconcplot(r,sol,c_mat,["10min-M" "30min-M" "1hr-M" "10min-J" "30min-J" "1hr-J"])
display(concentrationplot)

errorplot=errorconcplot(r,sol,c_mat)
display(errorplot)

accumplot=doubleaccumplot(sol,c_mat,["MATLAB Solver" "Julia Solver"])
display(accumplot)

accumerrorplot=Accumerrorplot(sol,c_mat)
display(accumerrorplot)
