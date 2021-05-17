load Milankovitch.dat                           %Load Insolation Data
LR04 = dlmread('LR04stack.txt', '\t', 5, 0);    %Load LR04 stack for comparison

time = Milankovitch(:,1); 
ins = Milankovitch(:,3);
ins_dm = (ins-mean(ins));                       %de-mean insolation input data
ins_dm_norm = ins_dm./(max(ins_dm)-min(ins_dm));%normalize by range
ins_dm_vnorm = ins_dm./20;

%% Calder Model
dt = 1e3;                       %time step length
nt = 2001;                      %number of time steps (to match with insolation input data)
V = nan.*ones(nt,1);            %pre-allocate volume vector
V(1) = 0;                       %initial condition for volume

%parameters
i0 = 502;       
km = 0.55e-5;
ka = 0.22e-5;

for t = 1:nt
    if(ins(t)>i0)
        k = km;
    else
        k = ka;
    end
    
    dVdt = -k*(ins(t)-i0);      %calculate RHS of model
    
    V(t+1) = V(t) + dVdt*dt;    %advance model one step using a forward-euler ODE approximation
end

figure(1);set(1,'units','normalized','position',[0 0.1 0.5 0.38]);
[ax,h1,h2] = plotyy(time./1e3,ins,time./1e3,V(1:end-1)-mean(V(1:end-1))); %plot insolation on left axis and model volume output on right axis
hold(ax(2));
plot(ax(2),-LR04(:,1),LR04(:,2)-mean(LR04(:,2)),'Color',[0 0 1],'linewidth',3) %add LR04 to right axis for comparison
xlim([-500 0])
set(h1,'linewidth',3,'Color',[0.7 0.7 0.7])
set(h2,'linewidth',3)
xlabel('time (kyr)','fontsize',20)
ylabel(ax(1),'Insolation @ 65N (W/m^2)','fontsize',20,'Color',[0.7 0.7 0.7])
ylabel(ax(2),'Ice "Volume"','fontsize',20)
set(ax,'fontsize',20,'XLim',[-500 0])

%% Imbrie and Imbrie Model
dt = 1e3;
nt = 2001;
V = nan.*ones(nt,1);
V(1) = 0;

taum = 42e3;
taua = 10e3;

for t = 1:nt
    if(V(t)>ins_dm_norm(t))
        tau=taum;
    else
        tau=taua;
    end
    
    dVdt = (ins_dm_norm(t)-V(t))/tau;
    
    V(t+1) = V(t) + dVdt*dt;
end

figure(2);set(2,'units','normalized','position',[0 0.1 0.5 0.38]);
[ax,h1,h2] = plotyy(time./1e3,ins,time./1e3,V(1:end-1)-mean(V(1:end-1)));
hold(ax(2));
plot(ax(2),-LR04(:,1),(LR04(:,2)-mean(LR04(:,2)))./10,'Color',[0 0 1],'linewidth',3)
xlim([-500 0])
set(h1,'linewidth',3,'Color',[0.7 0.7 0.7])
set(h2,'linewidth',3)
xlabel('time (kyr)','fontsize',20)
ylabel(ax(1),'Insolation @ 65N (W/m^2)','fontsize',20,'Color',[0.7 0.7 0.7])
ylabel(ax(2),'Ice "Volume"','fontsize',20)
set(ax,'fontsize',20,'XLim',[-500 0])

%% Paillard Model
ins_dm_vnorm2 = ins_dm_vnorm(1009:end);
time2 = time(1009:end);

dt = 1e3;
nt = length(ins_dm_vnorm2);
V = nan.*ones(nt,1);
V(1) = 0;

taug = 50e3;
tauf = 25e3;
taui = 10e3;
i0 = -0.75;
i1 = 0;
Vmax = 1;

st = 3; %st=1 is full glacial, st=2 is mild glacial, st=3 is interglacial
for t = 1:nt
    if(ins_dm_vnorm2(t)<i0 && st==3)
        Vr = 1;
        taur = taug;
        st = 1;
    else if(V(t)>Vmax && st==1)
        Vr = 1;
        taur=taug;
        st=2;
        else if(ins_dm_vnorm2(t)>i1 && st==2)
            Vr = 0;
            taur=taui;
            st=3;
            end;end;end
    
    dVdt = (Vr-V(t))/taur - ins_dm_vnorm2(t)/tauf;
    
    V(t+1) = V(t) + dVdt*dt;
end

figure(3);set(3,'units','normalized','position',[0 0.1 0.5 0.38]);
[ax,h1,h2] = plotyy(time2./1e3,ins_dm_vnorm2,time2./1e3,V(1:end-1)-mean(V(1:end-1)));
hold(ax(2));
plot(ax(2),-LR04(:,1),(LR04(:,2)-mean(LR04(:,2)))./2,'Color',[0 0 1],'linewidth',3)
xlim([-500 0])
set(h1,'linewidth',3,'Color',[0.7 0.7 0.7])
set(h2,'linewidth',3)
xlabel('time (kyr)','fontsize',20)
ylabel(ax(1),'Insolation @ 65N (W/m^2)','fontsize',20,'Color',[0.7 0.7 0.7])
ylabel(ax(2),'Ice "Volume"','fontsize',20)
set(ax,'fontsize',20,'XLim',[-500 0])