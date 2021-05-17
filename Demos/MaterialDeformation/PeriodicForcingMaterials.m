%% Ice parameters

E = 1e10;      %elastic modulus (in Pa)
mu = 2e14;     %linearized viscosity of ice 

%% Periodic Stress Forcing for Ice

ts = 3600*(0:1:24*4);             %time steps for simulations  
period = 3600*24;                 %period of forcing
sigma = 1e4*sin(2*pi*ts/period);  %a few days of diurnal tidal forcing

nt = length(ts);
dt = ts(2)-ts(1);                 %time step

%% Elastic Material
strain_elas = 0;                  %initial condition

for t=1:nt
    strain_elas(t) = sigma(t)/E;  %elastic material
end

figure(1);set(1,'units','normalized','position',[0 0.1 0.5 0.4]);
[ax,h1,h2] = plotyy(ts./(3600*24),sigma,ts./(3600*24),strain_elas);hold on
set(h1,'linewidth',3,'Linestyle','--','Color','k');
set(h2,'linewidth',3,'Linestyle','-','Color','r');
xlabel('time (days)','fontsize',20);ylabel(ax(1),'Applied Stress');ylabel(ax(2),'Strain (unitless)','fontsize',20)
set(ax,'fontsize',20)
ylim(ax(2),[-2e-6 2e-6])

legend('Stress Forcing','Elastic Deformation')
%% Viscous Material
strain_visc = -1e-6;

for t=1:nt-1
    strainrate_visc = sigma(t)/mu;
    
    strain_visc(t+1) = strain_visc(t) + strainrate_visc*dt;
end

hold(ax(2))
plot(ax(2),ts./(3600*24),strain_visc,'b','linewidth',3);hold on

legend('Stress Forcing','Elastic Deformation','Viscous Deformation')
%% Viscoelastic Material
strain_ve = 0;
dsigma_dt = diff(sigma)./diff(ts);

for t=1:nt-1
    strainrate_ve = (sigma(t)/mu) + (dsigma_dt(t)/E);
    
    strain_ve(t+1) = strain_ve(t) + strainrate_ve*dt;
end


plot(ax(2),ts./(3600*24),strain_ve,'m','linewidth',3);hold on
legend('Stress Forcing','Elastic Deformation','Viscous Deformation','Viscoelastic Deformation')

