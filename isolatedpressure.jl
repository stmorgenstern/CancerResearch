function Isolated_Pressure(N,R,Lpt,Svt,Kt,Pvv,at)
  #define position space
  r= (1/R)*(range(0,stop=R,length=N))
  dr=1/(N-1)
  #initialize A matrix and F vector
  A=zeros(N,N)
  F=zeros(N,1)
  #Fill in the A matrix interior values
  M=N
  for i=2:(M-1)
    A[i,i-1]= -1*(r[i]^-1)*dr^-1 + dr^-2
    A[i,i]=-2*(dr^-2)-(at[i]^2)
    A[i,i+1]=(r[i]^-1)*dr^-1 + dr^-2
    F[i]= -(at[i]^2)*Pvv
  end
  #Apply boundary conditions
  A[1,1]=-2*dr^-2 -(at[1]^2)
  A[1,2]=2*dr^-2
  F[1]=-(at[1]^2)*Pvv
  A[N,N]=1
  # Solve the linear system for Pressure
  P=A\F
  return P
end
