function P = Isolated_Pressure_19a(N,R,Lpt,Svt,Kt,Pvv)
% Pressure Profile
att = R*sqrt(Lpt*Svt/Kt); % parameter alpha for tumor, ignoring lymphatics
% an = R*sqrt(Lpn*Svn/Kn);  % parameter alpha for normal, ignoring lymphatics

% theta=(Kt/Kn)*(att*cosh(att)-sinh(att));
% phi=(1+an)*sinh(att);

r = linspace(0,R,N);        % vector of spatial grid points
r = r./R;
dr = 1./(N-1);              % distance between grid points

A = zeros(N,N);
F = zeros(N,1);

at = zeros(N,1);

M = N;   % for isolated tumor

% allows for potential spatial variation in tumor alpha
for i=1:M
    at(i) = att;%  + 2*rand; 
end

% % Analytical solution from Baxter's article
% for i = 2:M
%     po(i)=  1 - sinh(att*r(i))/r(i)/sinh(att); 
%     vo(i) = (att*r(i)*cosh(att*r(i)) - sinh(att*r(i)))/sinh(att)/r(i)^2;
% end

for i = 2:M-1
    A(i,i-1) = -1./r(i)/dr + 1./dr^2;
    A(i,i) = -2./dr^2 - (at(i)^2);
    A(i,i+1) = 1./r(i)/dr + 1./dr^2;
    F(i) = -(at(i)^2)*Pvv;
end

% BOUNDARY CONDITINS FOR ISOLATED TUMOR
A(1,1) = -2/dr^2 - (at(1)^2);
A(1,2) = 2./dr^2;
F(1) = -at(1)^2*Pvv;

A(N,N) = 1.;

% Solution of linear problem for pressure distribution
P = A\F;
end