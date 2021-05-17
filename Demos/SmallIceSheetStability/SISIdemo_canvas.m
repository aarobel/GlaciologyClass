%% Set parameters

sigma = 2;              %ice sheet shape parameter
a_0   = -1;             %SL SMB
beta  = 5e-3;           %vertical SMB gradient

%% Run Simple Model
dt = 1000;              %time step length
nt = 20;               %number of time steps

L = nan.*ones(nt,1);    %pre-allocate
L(1) = 100e3;           %set initial condition

for t = 1:nt            %loop for forward Euler
    dLdt = (a_0/sigma)*(L(t)^(1/2)) + (2/3)*beta*L(t);  %calculate LHS of model ODE
    
    L(t+1) = L(t) + dLdt*dt; %advance forward euler
end

figure;plot(dt*(0:nt)./1e3,L./1e3,'linewidth',3)
xlabel('Time (kyr)','fontsize',20)
ylabel('Ice Sheet Extent (km)','fontsize',20)
set(gca,'fontsize',20)
