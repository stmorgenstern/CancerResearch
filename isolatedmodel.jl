include("isolatedpressure.jl")
include("soluteperm.jl")
using DifferentialEquations
#Define the discretized system of time dependent ODEs
function MST!(dc,c,p,t)
    #dc = preallocated vector corresponding to f(c,t)= dC/dt
    #c = concentration vector at all spatial points
    #p = vector of parameters as define in modelanalysis
    #t = current time
    #Unpackage parameters
    P,N,sigma,Perm,Lpt,Svt,Kt,Pvv,Pv,D,r,dr,kd=p
    c0=1 #initial concentration
    tspan=3600 #timespan
    cv=c0*exp(-t*tspan/kd) #time-dependent exponential decay term
    Pe=Lpt*Pv*(Pvv-P[1])*(1-sigma)/Perm #Initial Peclet value
    dc[1]= tspan*(2*D*(c[2]-c[1])*dr^-2 + Lpt*Svt*Pv*(Pvv-P[1])*(1-sigma)*(cv*exp(Pe)-c[1])/(exp(Pe)-1)) # R=0 boundary condition
    for j=2:N-1
        #Define interior nodes
        Pe=Lpt*Pv*(Pvv-P[j])*(1-sigma)/Perm #Peclet number
        dc[j]=tspan*(((2*D/r[j])*((c[j+1]-c[j])/dr))+(D*(c[j+1]-2*c[j]+c[j-1])/dr^2)+ Kt*((P[j+1]-P[j-1])/(2*dr))*((c[j+1]-c[j])/dr)+ Lpt*Svt*Pv*(Pvv-P[j])*(1-sigma)*(cv*exp(Pe)-c[j])/(exp(Pe)-1))
    end
    dc[N]=0 # concentration at the outtermost boundary
end
function MSTHighPeclet!(dc,c,p,t)
    #dc = preallocated vector corresponding to f(c,t)= dC/dt
    #c = concentration vector at all spatial points
    #p = vector of parameters as define in modelanalysis
    #t = current time
    #Unpackage parameters
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
end


function Isolated_Model(N,Kt,Lpt,Svt,D,sigma,Perm,R,Pv,Pvv,kd,n_nodes,at)
    #define radial position space
    r= (1/R)*(range(0,stop=R,length=N))
    dr=1/(N-1)
    #Get Pressure profile
    P = Isolated_Pressure(N,R,Lpt,Svt,Kt,Pvv,at)
    #define initial concentration values
    c0=zeros(N,1)
    c0[N]=0.
    #define timespan
    time_end =1.
    tspan=(0.,time_end)
    #package parameters
    Pe=Lpt*Pv*(Pvv-P[end])*(1-sigma)/Perm
    p=P,N,sigma,Perm,Lpt,Svt,Kt,Pvv,Pv,D,r,dr,kd
    #Define and solve system of ODEs using time constant of 1/n_nodes and the QNDF solver for stiff systems
    dt=(1/(n_nodes))
    prob=ODEProblem(MST!,c0,tspan,p)
    sol=solve(prob,QNDF(),saveat=dt)
    return sol
end
function Isolated_Model_HP(N,Kt,Lpt,Svt,D,sigma,Perm,R,Pv,Pvv,kd,n_nodes,at)
    #define radial position space
    r= (1/R)*(range(0,stop=R,length=N))
    dr=1/(N-1)
    #Get Pressure profile
    P = Isolated_Pressure(N,R,Lpt,Svt,Kt,Pvv,at)
    #define initial concentration values
    c0=zeros(N,1)
    c0[N]=0.
    #define timespan
    time_end =1.
    tspan=(0.,time_end)
    #package parameters
    p=P,N,sigma,Perm,Lpt,Svt,Kt,Pvv,Pv,D,r,dr,kd
    #Define and solve system of ODEs using time constant of 1/n_nodes and the QNDF solver for stiff systems
    dt=(1/(n_nodes))
    prob=ODEProblem(MSTHighPeclet!,c0,tspan,p)
    sol= solve(prob,Rosenbrock23(),saveat=dt)
    return sol
end

mean(vec)=sum(vec)/length(vec)
function Accumulation_Model(sol,n_spatial,n_time)
    #solution a_ij in sol represents the concentration at radius r_i and time t_j
    if n_time <=100 #if less than 100 time nodes preserve original dimension
        t = (60/n_time)*[i for i=1:n_time]
        c_model=zeros(n_time)
        for j = 1:n_time
            c_model[j] = mean(sol[:,j])
        end
    else#if greater than than 100 time nodes, sample the time domain 100 times
        t = (60/100)*[i for i=1:100]
        c_model=zeros(100)
        spacingfactor = n_time ÷ 100
        for j = 1:100
            c_model[j] = mean(sol[:,spacingfactor*j])
        end
    end
    return t, c_model
end
function MSTHighPecletEE!(c,p,t)
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
function ExplicitEulerMSTHighPeclet!(t0,c0,dt,p)
    n_time=Int64(1/dt)
    n_space = length(c0)
    c_out = zeros(n_space,n_time)
    t_out=zeros(n_time,1)
    c_out[:,1]=c0
    t_out[1]=t0
    t=t0
    for j = 2:n_time
        f = MSTHighPecletEE!(c_out[:,j-1],p,t)
        c_out[:,j]= c_out[:,j-1] + dt*f
        t=t0+dt
        t_out[j]=t
    end
    return c_out,t_out
end
