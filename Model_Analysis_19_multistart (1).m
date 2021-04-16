function Model_Analysis_19_multistart
clc
% Input Peff in data set
% idx = 1 -- control
% idx = 2 -- 3mg/kg
% idx = 3 -- 30mg/kg
idx = 1;

n_nodes = 2000;
n_time = 100;
N = 100;
%tspan = 3600;


num = 10;
Kt_set = linspace(1.0e-6,1.0e-8,num);
Lpt_set = linspace(1.0e-7,1.0e-5,num);
Peff_set = zeros(length(Kt_set),length(Lpt_set)); 
fval_set = zeros(length(Kt_set),length(Lpt_set));

tic
for a = 1:length(Kt_set)
   for b = 1:length(Lpt_set) %1:length(Lpt_set)
w = [a,b];
% Input known metabolic parameters
% Kt = 9.0e-7
% Lpt = 1.5e-6
Kt = Kt_set(a);                  % Hydraulic conductivity of tumor 
Lpt = Lpt_set(b) ;               % hydraulic conductivity of tumor vessels
Svt = 200;                         % tumor vascular density
D = 1.375e-07;                     % solute diffusion coefficient
% sigma = sigma_set(idx);          % solute reflection coefficient
% Perm = Perm_set(idx);            % solute vascular permeability 
rs = 32/2;                         % partical radius (nm)
[Perm,sigma] = solutePerm_19(Lpt,rs); % pore theory

R = 1.;                            % tumor radius (cm)
Pv = 25;                           % vascular pressure (mmHg)
Pvv = 1.;                          % vascular pressure dimensionless

r = linspace(0,R,N);               % vector of spatial grid points
r = r./R;
dr = 1./(N-1);                     % distance between grid points

kd = 1278*60;                      % blood circulation time of drug in hours;

co = 1;
%Peff = Peff_set(idx);

%% CALCULATE PRESSURE PROFILE %%%%%%
P = Isolated_Pressure_19(N,R,Lpt,Svt,Kt,Pvv);
% P(1)
% 
% figure
% plot(r,P);
% xlabel('Dimensionless radial Position');
% ylabel('Dimensionless Pressure');
% title('Pressure')

%% CALCULATE VELOCITY PROFILE %%%%%%
v = zeros(N,1);
for i=2:N-1
    v(i)=-(P(i+1)-P(i-1))/2/dr;
end

att = R*sqrt(Lpt*Svt/Kt); % parameter alpha for tumor, ignoring lymphatics
for i=1:N
    at(i) = att;%  + 2*rand; 
end

v(N)= -(-(at(N)*at(N) + P(N-1)*(-1./r(N)/dr + 1./dr^2))...
    /(1./r(N)/dr+1./dr^2) - P(N-1))/2/dr; % VELOCITY AT R=1 FOR ISOLATED TUMOR

%% RUN CALCULATION %%%%%%
[time,c] = Isolated_Model_19(N,Kt,Lpt,Svt,D,sigma,Perm,R,Pv,Pvv,kd,n_nodes);
time
c;
  
length(time)
%% Calculate Peff: Optimzation Problem  

Peff0_list = linspace(1e-6, 1e-8, 1000);

t_nodes = (10/n_time).*[1:n_time]; % min

for j = 1:n_time
    c_model(j) = mean(c(j*round(n_nodes/n_time),:));  %drug concentration at j hours
end

Lpt ;
Kt ;

c_model;

% Interpolates simulation data
t = (600/n_time).*[1:n_time]; % seconds
% c_data = (co*Peff*Svt*kd/(1-Peff*Svt*kd))*...
%     (exp(-Peff*Svt*t)-exp(-t/kd));

%Peff0_matrix = linspace(1e-8, 1e-6, 100); 
 
%for multi = 1:10 

%Peff0min = 1e-8
%Peff0max = 1e-6
%Peff0 = random('Peff0_matrix')

%initialPeff0 = randi([1 5]);
%Peff0 = initialPeff0*10^-8

f=@(Peff)obj(Peff,c_model,t,co,Svt,kd);              %declare object function 
lb = [0] ;  %lower bound on Peff
ub =  [1e-3] ; %upper bound on Peff 
A = [] ;                    %no linear inequality constraint LHS 
B = [];                     %no linear inequality constraint RHS 
Aeq = [];                   %no linear equality constraint LHS 
beq = [];                   %no linear equality constraint RHS 
p = [] ;            %no parameters for nonlinear inequality constraints 
%nonlcon = [];         %no nonlinear constraints (model) 

%Peff0 = [5e-8] ;  %reasonable initial guess
options = optimoptions(@fmincon, 'Algorithm', 'sqp');  %use SQP algorithm 

nonlcon = [];
    %call the SQP solver and solve the problem locally 

for multi = 1:10 
    
    Peff0 = Peff0_list(randi(numel(Peff0_list)));
    
    [Peff,fval,exitflag,output,lambda] = fmincon(f,Peff0,A,B,Aeq,beq,lb,ub,nonlcon,options);

    Peff_list(multi,1) = Peff;
    fval_list(multi,1) = fval;
    exitflag;

end 

Peff_list;
fval_list;

[fval,index] = min(fval_list) ;
Peff = Peff_list(index,1);

%  c_data = (co*Peff*Svt*kd/(1-Peff*Svt*kd))*...
%      (exp(-Peff*Svt*t)-exp(-t/kd));
% 
% figure
% plot(t_nodes,c_data,'r.')
% hold on
% plot(t_nodes,c_model)
% hold on
% title('Average Concentration Accumulation Over Time Domain')
% xlabel('Time (min)')
% ylabel('Concentration Accumulation')
% legend('Data','Model')

%Peff_out(multi) = Peff;
%fval_out(multi) = fval;

%end 

Peff_set(a,b) = Peff;
fval_set(a,b) = fval;

Lpt;
Kt;

if b == 10
    if a == 1
        c_data = (co*Peff*Svt*kd/(1-Peff*Svt*kd))*...
             (exp(-Peff*Svt*t)-exp(-t/kd));
        figure
        plot(t_nodes,c_data,'r.')
        hold on
        plot(t_nodes,c_model)
        hold on
        title('Average Concentration Accumulation Over Time Domain')
        xlabel('Time (min)')
        ylabel('Concentration Accumulation')
        legend('Data','Model')
    end 
end
         
% 


%  c_data = (co*Peff*Svt*kd/(1-Peff*Svt*kd))*...
%      (exp(-Peff*Svt*t)-exp(-t/kd));
% 
% figure
% plot(t_nodes,c_data,'r.')
% hold on
% plot(t_nodes,c_model)
% hold on
% title('Average Concentration Accumulation Over Time Domain')
% xlabel('Time (min)')
% ylabel('Concentration Accumulation')
% legend('Data','Model')

    
    end   % end of lpt loop 
end   %end of kt loop


Peff_set
fval_set

toc 

%rows: each row is same Lpt value
%columns: each column is the same kt value 

%fval is highest when Lpt = 1e-6, so the lowest Lpt value 
%fval is overall highest when Lpt = 1e-6 (lowest) and Kt = 9e-7 (highest) 
%% Plot Heat Map

figure 
h = heatmap(Lpt_set, Kt_set, Peff_set)
h.Title = 'High Pe Effective Permeability Heat Map' ;
h.XLabel = 'Vascular Hydraulic Conductivity (Lp)';
h.YLabel = 'Interstitial Hydraulic Conductivity (K)';
h.Colormap = parula;
h.GridVisible = 'off';
h.ColorMethod = 'none';


end 

%% Objective Function 
function f = obj(Peff,c_model,t,co,Svt,kd)

% c_data = (co*Peff*Svt*kd/(1-Peff*Svt*kd))*...
%     (exp(-Peff*Svt*t)-exp(-t/kd));

sse_it = zeros(length(t),1);
%sum of squared errors 
for m = 1:length(t)
    sse_it(m) = (((co*Peff*Svt*kd/(1-Peff*Svt*kd))*...
        (exp(-Peff*Svt*t(m))-exp(-t(m)/kd))) - c_model(m))^2;
end 
Peff ;
f = sum(sse_it);

end 




%% RESULTS ANALYSIS %%%%%%
% Derive source term
% t = (3600/n_nodes).*[1:n_nodes]'; % seconds
% cv = co.*exp(-t/kd);
% 
% S_conv = zeros(n_nodes,N);
% S_diff = zeros(n_nodes,N);
% Pe = Lpt.*Pv.*(Pvv-P).*(1-sigma)./Perm;
% for j = 1:N
%     S_diff(:,j) = Perm*Svt.*(cv-c(:,j)).*(Pe(j)/(exp(Pe(j))-1));
%     S_conv(:,j) = Lpt*Svt*(1-sigma).*(Pvv - P(j)).*cv;
% end
% 
% % Derive convective and diffusive term
% conv = zeros(n_nodes,N);
% diff = zeros(n_nodes,N);
% diff(:,1) = 2*D.*(c(:,2)-c(:,1))/dr^2;
% 
% for j = 2:N-1
%     conv(:,j) = Kt.*((P(j+1)-P(j-1))./(2*dr)).*((c(:,j+1)-c(:,j))./dr);
%     diff(:,j) = ((2*D/r(j)).*((c(:,j+1) - c(:,j))./dr))+...
%         (D.*(c(:,j+1)-2*c(:,j)+c(:,j-1))./dr^2);
% end
% 
% % derive c(N+1)
% S = Lpt*Svt*Pv*(Pvv-P(N))*(1-sigma)*(cv*exp(Pe(N))-c(j))/(exp(Pe(N))-1);
% c_end = ((D/dr^2).*(c(:,N)-c(:,N-1))-S)./(2*D/r(N)/dr+D/dr^2-Kt*v(N)/dr)+c(:,N);
% 
% conv(:,N) = -Kt.*v(N).*((c_end-c(:,N))./dr);
% diff(:,N) = ((2*D/r(j)).*((c_end - c(:,N))./dr))+...
%         (D.*(c_end-2*c(:,N)+c(:,N-1))./dr^2);

% %% Plot pressure and velocity profie 
% figure
% plot(r,P);
% xlabel('Dimensionless radial Position');
% ylabel('Dimensionless Pressure');
% title('Pressure')
% 
% figure
% plot(r,v);
% xlabel('Dimensionless radial Position');
% ylabel('Dimensionless Velocity');
% title('Velocity')
% 
% %% Plot concentration profile
% c_10min = c(round(n_nodes/6),:);
% c_30min = c(n_nodes/2,:);     %drug concentration at 0.5 hour
% c_1hr = c(n_nodes,:);  %drug concentration at 1 hour
% figure
% plot(r,c_10min,'r',r,c_30min,'g',r,c_1hr,'b')
% xlabel('Dimensionless Distance');
% ylabel('Dimensionless Concentration');
% legend('10min','30min','1hr')
% title('Concentration')

% %%%%%% Plot solute flow over space domain %%%%%%
% % 10 min
% dcdt_10min = MST(time(round(n_nodes/6)),c(round(n_nodes/6),:),...
%     P,N,sigma,Perm,Lpt,Svt,Kt,Pvv,Pv,D,r,dr,kd);
% figure
% plot(r,dcdt_10min,'k')


% %% Plot convection-diffusion over space domain 
% % 10 min
% conv_10min = conv(round(n_nodes/6),:);
% diff_10min = diff(round(n_nodes/6),:);
% S_conv_10min = S_conv(round(n_nodes/6),:);
% S_diff_10min = S_diff(round(n_nodes/6),:);
% figure
% plot(r,conv_10min,'b',r,diff_10min,'r',r,S_conv_10min,'g',r,S_diff_10min,'m')
% legend('Convection','Diffusion','S_c_o_n_v','S_d_i_f_f')
% title('Solute Transport at 10min Post-injection')
% 
% % plot the enlarged view
% conv_enlarge = conv_10min(1,[90:100]);
% diff_enlarge = diff_10min(1,[90:100]);
% S_conv_enlarge = S_conv(1,[90:100]);
% S_diff_enlarge = S_diff(1,[90:100]);
% r_enlarge = r(1,[90:100]);
% figure
% plot(r_enlarge,conv_enlarge,'b',r_enlarge,diff_enlarge,'r',...
%     r_enlarge,S_conv_enlarge,'g',r_enlarge,S_diff_enlarge,'m')
% axis([0.9 1 -inf inf])
% %legend('Convection','Diffusion','S_c_o_n_v','S_d_i_f_f')
% 
% % 30 min
% conv_30min = conv(n_nodes/2,:);
% diff_30min = diff(n_nodes/2,:);
% S_conv_30min = S_conv(n_nodes/2,:);
% S_diff_30min = S_diff(n_nodes/2,:);
% figure
% plot(r,conv_30min,'b',r,diff_30min,'r',r,S_conv_30min,'g',r,S_diff_30min,'m')
% legend('Convection','Diffusion','S_c_o_n_v','S_d_i_f_f')
% title('Solute Transport at 30min Post-injection')
% 
% % 1 hr
% conv_1hr = conv(n_nodes,:);
% diff_1hr = diff(n_nodes,:);
% S_conv_1hr = S_conv(n_nodes,:);
% S_diff_1hr = S_diff(n_nodes,:);
% figure
% plot(r,conv_1hr,'b',r,diff_1hr,'r',r,S_conv_1hr,'g',r,S_diff_1hr,'m')
% legend('Convection','Diffusion','S_c_o_n_v','S_d_i_f_f')
% % title('Solute Transport at 1hr Post-injection')
% %%%%%% Plot average convection-diffusion over the whole time horizon %%%%%%
% for j = 1:n_time
%     conv_span(j) = mean(conv(j*round(n_nodes/n_time),:));
%     diff_span(j) = mean(diff(j*round(n_nodes/n_time),:));
%     S_conv_span(j) = mean(S_conv(j*round(n_nodes/n_time),:));
%     S_diff_span(j) = mean(S_diff(j*round(n_nodes/n_time),:));
% end
% 
% t_nodes = (60/n_time).*[1:n_time]; % min
% figure
% plot(t_nodes,conv_span,'b',t_nodes,diff_span,'r',...
%     t_nodes,S_conv_span,'g',t_nodes,S_diff_span,'m')
% legend('Convection','Diffusion','S_c_o_n_v','S_d_i_f_f')
% xlabel('time (min)')
% title('Solute Transport for The Whole Time Horizon')
% 
%% Plot Model vs Data

% t_nodes = (60/n_time).*[1:n_time]; % min
% 
% for j = 1:n_time
%     c_model(j) = mean(c(j*round(n_nodes/n_time),:));  %drug concentration at j hours
% end
% 
% % Interpolates simulation data
% t = (3600/n_time).*[1:n_time]; % seconds
% c_data = (co*Peff*Svt*kd/(1-Peff*Svt*kd))*...
%     (exp(-Peff*Svt*t)-exp(-t/kd));
% 
% figure
% plot(t_nodes,c_data,'r.')
% hold on
% plot(t_nodes,c_model)
% hold on
% xlabel('time (min)')
% ylabel('data concentration vs model concentration')



% %% Plot Convection Diffusion Ratio
% for idx = 1:3
%     Kt = Kt_set(idx);         % Hydraulic conductivity of tumor 
%     Lpt = Lpt_set(idx);             % hydraulic conductivity of tumor vessels
%     sigma = sigma_set(idx);   % solute reflection coefficient
%     Perm = Perm_set(idx);            % solute vascular permeability 
%     Peff = Peff_set(idx);
%     P = Isolated_Pressure(N,R,Lpt,Svt,Kt,Pvv);
%     [time,c] = Isolated_Model(N,Kt,Lpt,Svt,D,sigma,Perm,R,Pv,Pvv,kd,n_nodes);
%     Pe = Lpt.*Pv.*(Pvv-P).*(1-sigma)./Perm;
%     for j = 1:N
%         S_diff(:,j) = Perm*Svt.*(cv-c(:,j)).*(Pe(j)/(exp(Pe(j))-1));
%         S_conv(:,j) = Lpt*Svt*(1-sigma).*(Pvv - P(j)).*cv;
%     end
%     for j = 1:n_time
%         S_conv_span(j) = mean(S_conv(j*round(n_nodes/n_time),:));
%         S_diff_span(j) = mean(S_diff(j*round(n_nodes/n_time),:));
%     end
%     Avg_S_conv(idx) = mean(S_conv_span);
%     Avg_S_diff(idx) = mean(S_diff_span);
% end
% C_D_ratio = [Avg_S_diff;Avg_S_conv]'./sum([Avg_S_diff;Avg_S_conv]',2);
% figure
% bar(C_D_ratio,0.4,'stacked')
% legend('Diffusion','Convection')
% xticklabels({'Control','3 mg/kg','30 mg/kg'})
% end