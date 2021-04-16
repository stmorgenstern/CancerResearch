function [Perm,sigma] = solutePerm_19a(Lpt,rs)
% calculate diffusion coefficient from Stoke's Einstein
kB = 1.380648*10^(-23);       % Boltzmann Constant (J/K)
T = 310.15;                   % Temperature K
eta = 3*10^(-5);              % viscosity of blood (mmHg-sec)
conv_mu  = 133.322365;        % (Pascal/mmHg)
etac = eta*conv_mu;            % Pascal-sec
pore_conv = 10^(-9);          % (m/nm)
r_partc = rs*pore_conv;       % radius (m)
D0 = kB*T/(6*pi*etac*r_partc)*1e4; % Diffusivity (cm^2/s)

% Bungay and Brenner
a = [-73/60,77293/50400,-22.5083,-5.6117,-0.3363,-1.216,1.647];
b = [7/60;-2227/50400;4.0180;-3.9788;-1.9215;4.392;5.006];
% Calculate the pore size
gamma = 1e-3;
eta = 3e-5;  % Blood viscosity (mmHg/s)
L = 5e-4;    % Vessel wall thickness (cm)
r_pore = sqrt(8*eta*L*Lpt/gamma)*1e7; % nm
% rs = 30; % solute radius (nm) 30nm for FITC
lambda = rs/r_pore;
t1 = 0;
t2 = 0;
p1 = 0;
p2 = 0;
for i = 1:7
    if i<3
        t1 = t1 + a(i)*(1-lambda)^i;
        p1 = p1 + b(i)*(1-lambda)^i;
    else
        t2 = t2 + a(i)*lambda^(i-3);
        p2 = p2 + b(i)*lambda^(i-3);
    end
end
Kt = 9/4*pi^2*sqrt(2)*(1-lambda)^(-5/2)*(1+t1)+t2;
Ks = 9/4*pi^2*sqrt(2)*(1-lambda)^(-5/2)*(1+p1)+p2;
Phi = (1-lambda)^2;
H = 6*pi*Phi/Kt;
W = Phi*(2-Phi)*Ks/(2*Kt);
Perm = gamma*H*D0/L;
sigma = 1 - W;
end