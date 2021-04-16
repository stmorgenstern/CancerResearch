function [time,c] = Isolated_Model_19a(N,Kt,Lpt,Svt,D,sigma,Perm,R,Pv,Pvv,kd,n_nodes)
% Solution of Baxter's model
% We ignore lymphatics and osmotic pressures


r = linspace(0,R,N);        % vector of spatial grid points
r = r./R;
dr = 1./(N-1);              % distance between grid points

P = Isolated_Pressure_19a(N,R,Lpt,Svt,Kt,Pvv);

% solution of transient mass transport model

%initial solute concentration
c_0 = zeros(N,1);
c_0(N) = 0;

time_end = 1;         % length of simulation (seconds)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ODE15s %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options = odeset('Stats','off');

t_span = linspace(0,time_end,n_nodes); % discretization w.r.t. time    

[time,c] = ode15s(@MST_19a,t_span,c_0,options,P,N,sigma,Perm,...
    Lpt,Svt,Kt,Pvv,Pv,D,r,dr,kd);

end