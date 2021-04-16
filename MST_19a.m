function [f] = MST_19a(t,c,P,N,sigma,Perm,Lpt,Svt,Kt,Pvv,Pv,D,r,dr,kd)
% Concentration Profile
co = 1; % dimensionless drug concentration
% kd = 24*3600; % blood circulation time of drug in seconds

tspan = 3600;


f = zeros(N,1);   
cv = co*exp(-t*tspan/kd);  % vascular concentration of the drug following exponential decay

%low peclect regime 
%Pe = Lpt*Pv*(Pvv-P(1))*(1-sigma)/Perm;
%f(1) = tspan*(2*D*(c(2)-c(1))/dr^2 + ... 
%    Lpt*Svt*Pv*(Pvv-P(1))*(1-sigma)*(cv*exp(Pe)-c(1))/(exp(Pe)-1));


%for j = 2:N-1
%    Pe = Lpt*Pv*(Pvv-P(j))*(1-sigma)/Perm;
%      f(j) = tspan*(((2*D/r(j))*((c(j+1)-c(j))/dr))+(D*(c(j+1)-2*c(j)...
%      + c(j-1))/dr^2)+ Kt*((P(j+1)-P(j-1))/(2*dr))*((c(j+1)-c(j))/dr)...
 %     + Lpt*Svt*Pv*(Pvv-P(j))*(1-sigma)*(cv*exp(Pe)-c(j))/(exp(Pe)-1));
%end

%f(N) = 0; 


 %high peclet regime 
 Pe = Lpt*Pv*(Pvv-P(1))*(1-sigma)/Perm;
 f(1) = tspan*(2*D*(c(2)-c(1))/dr^2 + ... 
      Lpt*Svt*Pv*(Pvv-P(1))*(1-sigma)*cv);
 
 for j = 2:N-1    
       f(j) = tspan*(((2*D/r(j))*((c(j+1)-c(j))/dr))+(D*(c(j+1)-2*c(j)...
       + c(j-1))/dr^2)+ Kt*((P(j+1)-P(j-1))/(2*dr))*((c(j+1)-c(j))/dr)...
       + Lpt*Svt*Pv*(Pvv-P(j))*(1-sigma)*cv);
 end
 
f(N) = 0;



end