clear

%% Parameters

b_max = 1;          % max SMB at top of glacier
beta = 5e-3;        % vertical SMB gradient
bx = 1e-1;          % bed slope
gamma = beta*bx;    % horizontal SMB gradient
H = 100;            % glacier thickness (constant)

dt = 1;
nt = 1000;

%% Run Mountain Glacier model to steady state

L = nan.*ones(nt+1,1);

L(1) = 5e3;

for t = 1:nt
    dLdt = (b_max*L(t) - 0.5*gamma*(L(t)^2))/H;
    
    L(t+1) = L(t) + dLdt*dt;        %first-order forward euler method for solving ODE
end

%% Run Mountain Glacier model with step change in forcing

b_t_init = b_max - gamma*L(end);

L(1) = L(end);
b_max_init = b_max;
b_max = 0.9;

for t = 1:nt
    dLdt = (b_max*L(t) - 0.5*gamma*(L(t)^2))/H;
    
    L(t+1) = L(t) + dLdt*dt;        %first-order forward euler method for solving ODE
end

time = linspace(0,nt*dt,nt+1);
figure;
plot(time,L,'b','linewidth',4);hold on
xlabel('Time (yr)','fontsize',30)
ylabel('Glacier Length (m)','fontsize',30)
set(gca,'fontsize',30)
%% Calculate linear response parameters
tau = -H/b_t_init;                                  % linear response time scale
Lprime = (b_max-b_max_init)*(-L(1)/b_t_init);       % linear response sensitivity

L_lin = L(1)+Lprime*(1-exp(-time/tau));             % predicted linear response
plot(time,L_lin,'k--','linewidth',4)