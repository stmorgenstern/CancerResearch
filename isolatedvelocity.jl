function Isolated_Velocity(n_spatial,pressure,r, dr)
    #initialize vector
    v=zeros(n_spatial,1)
    #Set interior nodes
    v[2:n_spatial-1]=[-.5*(P[i+1]-P[i-1])*dr^-1 for i=2:n_spatial-1]
    #set outtermost boundary nodes
    v[n_spatial]=-.5*(dr^-1)*(-1*((at[n_spatial])^2+P[n_spatial-1]*((-1/r[n_spatial])*dr^-1 + dr^-2))/((1/r[n_spatial])*dr^-1 + dr^-2) - P[n_spatial-1])
    return v
end
