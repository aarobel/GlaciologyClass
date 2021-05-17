%% Ice stream box model (Robel et al., JGR, 2013)
% This is the driver for the ice stream box model described in 
% Robel et al., JGR, 2013.
% In this script, model parameters and settings are specified and the model
% is integrated using ODE45 (though other MATLAB integrators may be used as
% well)
%
% See Ice_Stream_Box_Model_RHS for implementation of model equations.
% 
% Code written by Alex Robel, last update August, 2016
% Debugging help from Christian Schoof, Eli Tziperman, Eric DeGiuli and
% Elisa Mantelli

%% Set parameters
global o
p.year = 3600*24*365;   %seconds in a year
t_final = 1e4;          %total time of integration

p.L=500e3;              %ice stream length
p.W=40e3;               %ice stream width
p.n=3;                  %nye-glen law exponent
p.q_g = 0.07;           %geothermal heat flux
p.htill_init = 1;       %initial till thickness

p.T_s=21;               %ice surface temperature

p.rho_i = 917;          %density of glacial ice
p.L_f = 3.335e5;        %latent heat of fusion for ice
p.K_i = 2.1;            %thermal diffusivity of ice
p.A_f = 5e-25;          %nye-glen law rate factor
p.g = 9.81;             %acceleration due to gravity

p.e_c = 0.3;            %till consolidation void ratio
p.tau0 = 9.44e8;        %empirical till coefficient
p.c = 21.7;             %empirical till exponent
p.C_ice = 1.94e6;       %specific heat capacity of ice

p.A = p.L*p.W;              %ice stream area

p.eta_b = 10;           %basal ice layer thickness
p.h_t_min = 1e-3;       %minimum till thickness (prevents unfrozen till thickness from going to zero)

p.a = 0.1./p.year;      %accumulation rate

%% Set initial conditions and intergration time

p.ic=[650 .6 p.htill_init 0];   %Initial conditions
p.tspan=[0,p.year*t_final];     %time steps

%% Run Model

options = odeset('RelTol',1e-6,'AbsTol',1e-6);      %set ode integration settings
[time,T] = ode45(@(t,X) Ice_Stream_Box_Model_RHS(t,X,p),p.tspan,p.ic,options);  %integrate box model

%% Make some plots
figure(1);set(1,'units','pixels','position',[0 0 702 1202])
subplot(4,1,4);
plot((o.t(1:o.nt)./p.year),o.e(1:o.nt),'k','linewidth',5);hold on;
ylabel('Till Water Content (m)','fontsize',20);
set(gca,'fontsize',24)
% xlim([0 1870])
xlabel('Time (yr)','fontsize',24)

subplot(4,1,2);
plot((o.t(1:o.nt)./p.year),o.h(1:o.nt),'k','linewidth',5);hold on;
ylabel('Ice Thickness (m)','fontsize',20);
set(gca,'fontsize',24)
% xlim([0 1870]);
ylim([400 700])
xlabel('Time (yr)','fontsize',24)

subplot(4,1,1);
plot((o.t(1:o.nt)./p.year),p.year.*o.u(1:o.nt),'k','linewidth',5);hold on;
ylabel('Sliding Velocity (m/yr)','fontsize',20);
set(gca,'fontsize',24)
% xlim([0 1870]);
ylim([0 300])
xlabel('Time (yr)','fontsize',24)

subplot(4,1,3);
plot((o.t(1:o.nt)./p.year),1e3.*o.qg(1:o.nt),'k','linewidth',5);hold on;
plot((o.t(1:o.nt)./p.year),1e3.*o.vh(1:o.nt),'k--','linewidth',5);hold on;
plot((o.t(1:o.nt)./p.year),1e3.*o.fh(1:o.nt),'k:','linewidth',5);hold on;
legend('Geothermal','Vertical Conduction','Frictional')
ylabel('Heat Flux (mW/m^2)','fontsize',20);
set(gca,'fontsize',24)
% xlim([0 1870]);
ylim([0 80])
xlabel('Time (yr)','fontsize',24)

% subplot(5,1,4);
% plot((o.t(1:o.nt)./p.year)-4690,o.htill(1:o.nt),'k','linewidth',2);hold on;
% ylabel('Till Thickness (m)','fontsize',16);
% set(gca,'fontsize',16)
% xlim([0 1870])
% 
% subplot(5,1,5);
% plot((o.t(1:o.nt)./p.year)-4690,o.tb(1:o.nt),'k','linewidth',2);hold on;
% ylabel('Basal Temp (K dep p-m)','fontsize',16);
% set(gca,'fontsize',16)
% xlim([0 1870])