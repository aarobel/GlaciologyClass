function rhs=Ice_Stream_Box_Model_RHS(t,X,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ice stream box model (Robel et al., JGR, 2013)
% This is the dynamics function for the ice stream box model described in 
% Robel et al., JGR, 2013. 
% In this script, model parameters and settings are specified and the model
% is integrated using ODE45 (though other MATLAB integrators may be used as
% well)
%
% See Ice_Stream_Box_Model.m for model driver.
% 
% Code written by Alex Robel, last update August, 2016
% Debugging help from Christian Schoof, Eli Tziperman, Eric DeGiuli and
% Elisa Mantelli

persistent nt
global o

%% Read in variables
h = X(1);
e = X(2); 
h_till = X(3);
T_b = X(4);

%% Set parameters
n=p.n;
q_g = p.q_g;
T_s=p.T_s;
rho_i = p.rho_i;
L_f = p.L_f;
K_i = p.K_i;
A_f = p.A_f;
g = p.g;

e_c = p.e_c;
tau0 = p.tau0;
c = p.c;
C_ice = p.C_ice;

L = p.L;
W = p.W;
A = p.A;

eta_b = p.eta_b;
h_t_min = p.h_t_min;

a = p.a;

%% Diagnostic model equations

deltaT = T_s-T_b;               %difference between basal and surface ice temperature
e = max([e,e_c]);               %make sure till void ratio doesn't go below threshold
h_till = max([h_till,h_t_min]);h_till=min([h_till p.htill_init]);   %make sure unfrozen till thickness stays within bounds
T_b = max([T_b 0]);             %make sure basal ice temperature never goes above zero

tau_d = rho_i*g*(h^2/L);        %calculate driving stress
tau_f = tau0*exp(-c*e);         %calculate basal shear stress

U = (A_f/256)*(W^(n+1))*((max([(tau_d-tau_f)/h,0]))^n); %calculate centerline ice stream velocity

%% Prognostic model equations (see Robel et al., JGR, 2013 for details)

if((h_till==h_t_min && T_b==0 && ((tau_f*U) + q_g - (K_i*deltaT/h))<0) || (h_till==h_t_min && T_b>0))   %till is frozen
    U=0;
    dedt = 0;
    dhtdt = 0;
    dTbdt = -1*(1/(eta_b*C_ice))*((tau_f*U) + q_g - (K_i*deltaT/h));
    
    else if((e == e_c && h_till==p.htill_init && ((tau_f*U) + q_g - (K_i*deltaT/h))<0) ||  (e == e_c && h_till<p.htill_init))   %basal ice is temperate, till is partially frozen
        U=0;
        dedt = 0;
        dhtdt = ((tau_f*U) + q_g - (K_i*deltaT/h))/(L_f*rho_i);
        dTbdt = 0;
        
        else                             %basal ice is temperate, till is thawed, void ratio varies
            dedt = ((tau_f*U) + q_g - (K_i*deltaT/h))/(h_till*L_f*rho_i);
            dhtdt = 0;
            dTbdt = 0;
    end
end

dhdt = ((A*a)/(L*W))-(h*U/L);

disp(['% of integration finished: ' num2str(100*(t/p.tspan(2))) '%']);

%% save diagnostics to output variable:
% first, advance time step used for saving diagnostics.
if isempty(nt); nt=1; else; nt=nt+1; end
if nt==1; init=1; else; init=0; end

max_nt=10000;

%% matlab sometimes makes a trial step, and then goes back to get
%% more accurate results, take care of this by not saving trial
%% steps in output arrays:
%%fprintf(1,'%g, %d; ',t, nt);
while nt>2 && t<=o.t(nt-1);
  nt=nt-1;
  %%fprintf(1,'\n !! \n');
end
o.nt=nt;

if init; ZZ=zeros(max_nt,1)*NaN; end

if init; o.h=ZZ; end;   o.h(nt) = h;
if init; o.e=ZZ; end;   o.e(nt) = e;
if init; o.u=ZZ; end;   o.u(nt) = U;
if init; o.htill=ZZ; end;   o.htill(nt) = h_till;
if init; o.tb=ZZ; end;   o.tb(nt) = T_b;
if init; o.ts=ZZ; end;   o.ts(nt) = T_s;

if init; o.tauf=ZZ; end;   o.tauf(nt) = tau_f;
if init; o.fh=ZZ; end;   o.fh(nt) = tau_f*U;
if init; o.qg=ZZ; end;   o.qg(nt) = q_g;
if init; o.vh=ZZ; end;   o.vh(nt) = K_i*deltaT/h;

if init; o.dedt=ZZ; end;   o.dedt(nt) = dedt;
if init; o.dhtdt=ZZ; end;   o.dhtdt(nt) = dhtdt;
if init; o.dTbdt=ZZ; end;   o.dTbdt(nt) = dTbdt;

if init; o.t=ZZ; end;   o.t(nt)=t;

rhs= [dhdt;
    dedt
    dhtdt
    dTbdt];