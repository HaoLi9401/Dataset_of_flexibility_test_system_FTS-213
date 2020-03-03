% The flexibility resource planning codes for the FTS-213 test system.
% Author: Hao Li, email: h-l16@mails.tsinghua.edu.cn
% updated: 2020_March_02

% Yalmip and Gurobi packages are required to run the codes.
% Yalmip package: https://yalmip.github.io/
% Gurobi package: https://www.gurobi.com/

load('data_20200214.mat');
load('month_week');
load('DCtrans');
load('TieLine2');
load('Cluster_day_number');

a = case213_20200213;
%%%%%%%%%%%%%%%%%%
a.RE(:,7) = a.RE(:,7)*3;
%%%%%%%%%%%%%%%%%%

%% Parameter setting
Horizon = 24;
TimeScale = 1;
CoeffReseve_load = 0.03;% reserve coefficient for load
CoeffReserve_VRE = 0.05;% reserve coefficient for VRE
efficiency_TES = 0.98;% efficiency of CSP TES(Thermal energy storage)
efficiency_ele_ther = 0.415;% efficiency of electric power output over input thermal power of CSP power block
Capacity_TES_CSP = 2.4;% TES capacity relative to the capacity of power block of CSP
ratio_TES_t0 = 0.5;% initial energy ratio of TES

LT_ES_Battery = 10; % life time of battery 
LT_ES_HydroPump = 30;% life time of pumped hydro
LT_ES_CAES = 25;% life time of CAES
LT = 40;% life time of others

w = 0.08;
txl = w*((1+w)^LT)/(((1+w)^LT)-1);
txl_ES_Battery = w*((1+w)^LT_ES_Battery)/(((1+w)^LT_ES_Battery)-1);
txl_ES_HydroPump = w*((1+w)^LT_ES_HydroPump)/(((1+w)^LT_ES_HydroPump)-1);
txl_ES_CAES = w*((1+w)^LT_ES_CAES)/(((1+w)^LT_ES_CAES)-1);
theta_Load = 3/7*1000/1e6;% M$/MWh 
theta_RE = 0.42/7*1000/1e6;% M$/MWh 
cost_thermal_coal_peak_depth = 3e4/1e6; % M$/MW
cost_thermal_gas_peak_depth = 2e5/1e6; % M$/MW

cost_ES_Battery = 0.828; % M$/MW
cost_ES_HydroPump = 0.165*16; % M$/MW
cost_ES_CAES = 0.105*16;% M$/MW
cost_DR1_Capacity = 80000/7/1e6; % M$/MW, load shifting
cost_DR2_Capacity = 100000/7/1e6; % M$/MW, load shedding
cost_DR1_Energy = 350/7/1e6;%compensate cost，M$/MWh, load shifting
cost_DR2_Energy = 450/7/1e6;%compensate cost，M$/MWh, load shedding

% ES efficiency 
efficiency_ES_Battery = 0.86;
efficiency_ES_HydroPump = 0.8;
efficiency_ES_CAES = 0.52;

E_ES0 = 0.4; %initial energy rate of ES
ratio_DR1 = 0.05; % potential of load shifting DR
ratio_DR2 = 0.05; % potential of load shedding DR
M = 1e6;


% Minimum stable level
% Since the on/off status is neglected in the planning model,
% decrease the minimum stable level as compensation.
a.gen(:,10) = a.gen(:,10)*0.8;
a.CSP(:,10) = a.CSP(:,10)*0;

N = length(a.bus(:,1));
Nbr = length(a.branch(:,1));
Ngen = length(a.gen(:,1));
NRE = length(a.RE(:,1));
NCSP = length(a.CSP(:,1));
NStorage = size(a.Storage,1);
NPumpHydro = size(a.HydroPump,1);

Area_bus = [];
for b = 1:N-1
    if a.bus(b,11)-a.bus(b+1,11)==-1
        Area_bus = [Area_bus;b];
    end
end
Area_bus = [Area_bus;N];

M_bus_G = zeros(N,Ngen);
for row = 1:N
    if abs(find(a.gen(:,1)==row))>0
        M_bus_G(row,find(a.gen(:,1)==row)) = 1;
    end
end

M_bus_RE = zeros(N,NRE);
for row = 1:N
    if abs(find(a.RE(:,1)==row))>0
        M_bus_RE(row,find(a.RE(:,1)==row)) = 1;
    end
end

M_bus_CSP = zeros(N,NCSP);
for row = 1:N
    if abs(find(a.CSP(:,1)==row))>0
        M_bus_CSP(row,find(a.CSP(:,1)==row)) = 1;
    end
end

M_bus_BatteryAndCAES = zeros(N,NStorage);
for row = 1:N
    if abs(find(a.Storage(:,1)==row))>0
        M_bus_BatteryAndCAES(row,find(a.Storage(:,1)==row)) = 1;
    end
end
M_bus_PumpHydro = zeros(N,NPumpHydro);
for row = 1:N
    if abs(find(a.HydroPump(:,1)==row))>0
        M_bus_PumpHydro(row,find(a.HydroPump(:,1)==row)) = 1;
    end
end

% ------ Hydro power ------
PHydro_max = zeros(365,1);
PHydro_min = zeros(365,1);
for i = 1:365
    if i<=31
        PHydro_max(i) = PHydro_all(1,3);
        PHydro_min(i) = PHydro_all(3,3);
    else
        for j = 2:12
            if (i>=month_week(j-1,3)+1)&&(i<=month_week(j,3))
                PHydro_max(i) = PHydro_all(1,j+2);
                PHydro_min(i) = PHydro_all(3,j+2);
                break;
            end
        end
    end
end
EHydro_month = zeros(12,1);
for i = 1:12
    EHydro_month(i) = month_week(i,2)*24*PHydro_all(2,i+2);
end
EHydro_day = zeros(5,365); 
mean_PLoad = zeros(365,5);
for Area = 1:5
    mean_PLoad(:,Area) = mean(PLoad_all{Area},2);% 没有算DC外送负荷
end
for i = 1:5
    for j = 1:365
        if j<=31
            EHydro_day(i,j) = EHydro_month(1)/sum(mean_PLoad(1:31,i))*mean_PLoad(j,i);
        else
            for k = 2:12
                if (j>=month_week(k-1,3)+1)&&(j<=month_week(k,3))
                    EHydro_day(i,j) = EHydro_month(k)/sum(mean_PLoad(month_week(k-1,3)+1:month_week(k,3),i))*mean_PLoad(j,i);
                    break;
                end
            end
        end
    end
end
% ------ Hydro power ------

gen_total = [];
for row = 1:size(a.gen,1)
    for j = 1:a.gen(row,32)
        gen_total = [gen_total;a.gen(row,:)];
    end
end
CSP_total = [];
for row = 1:size(a.CSP,1)
    for j = 1:a.CSP(row,32)
        CSP_total = [CSP_total;a.CSP(row,:)];
    end
end  
M_gen_gentotal = zeros(Ngen,size(gen_total,1));
for row1 = 1:Ngen
    for row2 = 1:size(gen_total,1)
        if a.gen(row1,[1,22]) == gen_total(row2,[1,22])
            M_gen_gentotal(row1,row2) = 1;
        end
    end
end
M_CSP_CSPtotal = zeros(NCSP,size(CSP_total,1));
for row1 = 1:NCSP
    for row2 = 1:size(CSP_total,1)
        if a.CSP(row1,[1,22]) == CSP_total(row2,[1,22])
            M_CSP_CSPtotal(row1,row2) = 1;
        end
    end
end											   


%% Typical days
day = (Cluster_day_number(:,1))'
N_scenario = length(day);

% ------ Load ------
% --- Local load ---
PLoad_day0 = zeros(N,Horizon,N_scenario);% 实际负荷
for sce = 1:N_scenario
    for row = 1:N
        if a.bus(row,11) == 1
            PLoad_day0(row,:,sce) = a.bus(row,3)*PLoad_all{6}((day(sce)-1)*Horizon+1:day(sce)*Horizon);
        elseif a.bus(row,11) == 2
            PLoad_day0(row,:,sce) = a.bus(row,3)*PLoad_all{7}((day(sce)-1)*Horizon+1:day(sce)*Horizon);
        elseif a.bus(row,11) == 3
            PLoad_day0(row,:,sce) = a.bus(row,3)*PLoad_all{8}((day(sce)-1)*Horizon+1:day(sce)*Horizon);
        elseif a.bus(row,11) == 4
            PLoad_day0(row,:,sce) = a.bus(row,3)*PLoad_all{9}((day(sce)-1)*Horizon+1:day(sce)*Horizon);
        else
            PLoad_day0(row,:,sce) = a.bus(row,3)*PLoad_all{10}((day(sce)-1)*Horizon+1:day(sce)*Horizon);
        end
    end
end
PLoad_day = PLoad_day0;
% --- Local load ---
% --- External load ---
for sce = 1:N_scenario
    if day(sce) <= 31
        for b = 1:size(DCtrans,1)
            PLoad_day(find(a.bus(:,1)==DCtrans(b,1)),:,sce) = PLoad_day(find(a.bus(:,1)==DCtrans(b,1)),:,sce)+DCtrans(b,2)*0.6;
        end
    else
        for month = 2:12
            if day(sce)>month_week(month-1,3) && day(sce)<=month_week(month,3)
                for b = 1:size(DCtrans,1)
                    PLoad_day(find(a.bus(:,1)==DCtrans(b,1)),:,sce) = PLoad_day(find(a.bus(:,1)==DCtrans(b,1)),:,sce)+DCtrans(b,month+1)*0.6;
                end
            end
        end
    end
end
% --- External load ---
% ------ Load ------

% ------ Wind, PV, CSP ------
PRE_day = zeros(NRE,Horizon,N_scenario);
for sce = 1:N_scenario
    for row = 1:NRE
        if a.RE(row,3) == 3
            if a.RE(row,2) == 1
                PRE_day(row,:,sce) = a.RE(row,7)*PSolar_all{6}((day(sce)-1)*Horizon+1:day(sce)*Horizon);
            elseif a.RE(row,2) == 2
                PRE_day(row,:,sce) = a.RE(row,7)*PSolar_all{7}((day(sce)-1)*Horizon+1:day(sce)*Horizon);
            elseif a.RE(row,2) == 3
                PRE_day(row,:,sce) = a.RE(row,7)*PSolar_all{8}((day(sce)-1)*Horizon+1:day(sce)*Horizon);
            elseif a.RE(row,2) == 4
                PRE_day(row,:,sce) = a.RE(row,7)*PSolar_all{9}((day(sce)-1)*Horizon+1:day(sce)*Horizon);
            else
                PRE_day(row,:,sce) = a.RE(row,7)*PSolar_all{10}((day(sce)-1)*Horizon+1:day(sce)*Horizon);
            end
        else
            if a.RE(row,2) == 1
                PRE_day(row,:,sce) = a.RE(row,7)*PWind_all{6}((day(sce)-1)*Horizon+1:day(sce)*Horizon);
            elseif a.RE(row,2) == 2
                PRE_day(row,:,sce) = a.RE(row,7)*PWind_all{7}((day(sce)-1)*Horizon+1:day(sce)*Horizon);
            elseif a.RE(row,2) == 3
                PRE_day(row,:,sce) = a.RE(row,7)*PWind_all{8}((day(sce)-1)*Horizon+1:day(sce)*Horizon);
            elseif a.RE(row,2) == 4
                PRE_day(row,:,sce) = a.RE(row,7)*PWind_all{9}((day(sce)-1)*Horizon+1:day(sce)*Horizon);
            else
                PRE_day(row,:,sce) = a.RE(row,7)*PWind_all{10}((day(sce)-1)*Horizon+1:day(sce)*Horizon);
            end
        end
    end
end
PCSP_day = zeros(NCSP,Horizon,N_scenario);% 实际CSP可用热功率
for sce = 1:N_scenario
    for row = 1:NCSP
            if a.CSP(row,31) == 1
                PCSP_day(row,:,sce) = a.CSP(row,9)*PCSP_all{6}((day(sce)-1)*Horizon+1:day(sce)*Horizon);
            elseif a.CSP(row,31) == 2
                PCSP_day(row,:,sce) = a.CSP(row,9)*PCSP_all{7}((day(sce)-1)*Horizon+1:day(sce)*Horizon);
            elseif a.CSP(row,31) == 3
                PCSP_day(row,:,sce) = a.CSP(row,9)*PCSP_all{8}((day(sce)-1)*Horizon+1:day(sce)*Horizon);
            elseif a.CSP(row,31) == 4
                PCSP_day(row,:,sce) = a.CSP(row,9)*PCSP_all{9}((day(sce)-1)*Horizon+1:day(sce)*Horizon);
            else
                PCSP_day(row,:,sce) = a.CSP(row,9)*PCSP_all{10}((day(sce)-1)*Horizon+1:day(sce)*Horizon);
            end
    end
end
% ------ Wind, PV, CSP ------

M_min_gen_thermal = zeros(size(gen_total(find(gen_total(:,22)~=1),:),1),Horizon,N_scenario);
M_max_gen_thermal = zeros(size(gen_total(find(gen_total(:,22)~=1),:),1),Horizon,N_scenario);
M_min_gen_CSP = zeros(size(CSP_total,1),Horizon,N_scenario);
M_max_gen_CSP = zeros(size(CSP_total,1),Horizon,N_scenario);
M_max_gen = zeros(size(gen_total,1),Horizon,N_scenario);
M_max_gen2 = zeros(size(gen_total,1),1,N_scenario);
M_min_gen2 = zeros(size(gen_total,1),1,N_scenario);
M_max_CSP2 = zeros(size(CSP_total,1),1,N_scenario);
M_min_CSP2 = zeros(size(CSP_total,1),1,N_scenario);
M_CostStartUp = zeros(size(gen_total,1),Horizon-1,N_scenario);
M_CostGeneration = zeros(size(gen_total,1),Horizon,N_scenario);
M_RampUp = zeros(size(gen_total,1),1,N_scenario);
M_RampDown = zeros(size(gen_total,1),1,N_scenario);
M_RampUp_CSP = zeros(size(CSP_total,1),1,N_scenario);
M_RampDown_CSP = zeros(size(CSP_total,1),1,N_scenario);
M_days_gentotal = zeros(size(gen_total,1),Horizon,N_scenario);
M_days_bus = zeros(N,Horizon,N_scenario);
M_days_RE = zeros(size(a.RE,1),Horizon,N_scenario);
for sce = 1:N_scenario
    M_min_gen_thermal(:,:,sce) = (gen_total(find(gen_total(:,22)~=1),10)*ones(1,Horizon));
    M_max_gen_thermal(:,:,sce) = (gen_total(find(gen_total(:,22)~=1),9)*ones(1,Horizon));
	M_min_gen_CSP(:,:,sce) = (CSP_total(:,10)*ones(1,Horizon));
    M_max_gen_CSP(:,:,sce) = (CSP_total(:,9)*ones(1,Horizon));
    M_max_gen(:,:,sce) = (gen_total(:,9)*ones(1,Horizon));
    M_CostStartUp(:,:,sce) = gen_total(:,28)*ones(1,Horizon-1);
    M_CostGeneration(:,:,sce) = gen_total(:,30)*ones(1,Horizon);
    M_RampUp(:,:,sce) = gen_total(:,23);
    M_RampDown(:,:,sce) = gen_total(:,24);
	M_RampUp_CSP(:,:,sce) = CSP_total(:,23);
	M_RampDown_CSP(:,:,sce) = CSP_total(:,24);
    M_max_gen2(:,:,sce) = gen_total(:,9);
    M_min_gen2(:,:,sce) = gen_total(:,10);
	M_max_CSP2(:,:,sce) = CSP_total(:,9);
    M_min_CSP2(:,:,sce) = CSP_total(:,10);
    M_days_gentotal(:,:,sce) = Cluster_day_number(sce,2)*ones(size(gen_total,1),Horizon);
    M_days_bus(:,:,sce) = Cluster_day_number(sce,2)*ones(N,Horizon);
    M_days_RE(:,:,sce) = Cluster_day_number(sce,2)*ones(NRE,Horizon);
end


%% Optimization
% Operation variables
PG_G = sdpvar(size(gen_total,1),Horizon,N_scenario,'full');% Hydro and thermal power
PG_RE = sdpvar(NRE,Horizon,N_scenario,'full');% Wind and PV power
PG_CSP = sdpvar(size(CSP_total,1),Horizon,N_scenario,'full');% CSP power
Et_TES = sdpvar(size(CSP_total,1),Horizon,N_scenario,'full');% TES energy of CSP
PC_Load = sdpvar(N,Horizon,N_scenario,'full');% Load shedding
onoff = ones(size(gen_total,1),Horizon,N_scenario);% on/off status of hydro and thermal units
onoff_CSP = ones(size(CSP_total,1),Horizon,N_scenario);% on/off status of CSP units
PG_G_bus = sdpvar(N,Horizon,N_scenario,'full');% injection power of hydro and thermal units into each bus
PG_RE_bus = sdpvar(N,Horizon,N_scenario,'full');% injection power of wind and PV units into each bus
PG_CSP_bus = sdpvar(N,Horizon,N_scenario,'full');% injection power of CSP units into each bus
PG_ES_BatteryAndCAES_bus = sdpvar(N,Horizon,N_scenario,'full');% injection power of battery and CAES into each bus
PG_ES_PumpHydro_bus = sdpvar(N,Horizon,N_scenario,'full');% injection power of pumped hydro into each bus
PG_ES_BatteryAndCAES = sdpvar(NStorage,Horizon,N_scenario,'full');% output power of battery and CAES 
E_ES_BatteryAndCAES = sdpvar(NStorage,Horizon,N_scenario,'full');% energy of battery and CAES 
PG_ES_PumpHydro = sdpvar(NPumpHydro,Horizon,N_scenario,'full');% output power of pumped hydro
E_ES_PumpHydro = sdpvar(NPumpHydro,Horizon,N_scenario,'full');% energy of pumped hydro 
PG_DR2 = sdpvar(N,Horizon,N_scenario,'full');% load changing power of load shedding DR 
PG_DR1 = sdpvar(N,Horizon,N_scenario,'full');% load changing power of load shifting DR 
abs_PG_DR1 = sdpvar(N,Horizon,N_scenario,'full');% Auxiliary variables to indicate the absolute value of load shifting DR
Branch = sdpvar(size(a.branch,1)+size(a.candi_AC,1)+size(a.candi_DC,1),Horizon,N_scenario,'full');% power flow
theta = sdpvar(N,Horizon,N_scenario,'full');% phases angle 
% Investment variables
invest_line_AC = binvar(size(a.candi_AC,1),1); % AC branch
invest_line_DC = binvar(size(a.candi_DC,1),1); % DC branch
invest_gen_peak_depth = sdpvar(size(gen_total,1),1,'full');% flexibility retrofit of power plants
invest_ES_BatteryAndCAES = sdpvar(NStorage,1,'full');% battery and CAES
invest_DR1 = sdpvar(N,1,'full');% load shifting DR
invest_DR2 = sdpvar(N,1,'full');% load shedding DR
invest_ES_PumpHydro = binvar(NPumpHydro,1,'full');% pumped hydro 

% Constraints
F = [];
% Upper and lower bounds of flexibility retrofit 
F = [F, gen_total(find(gen_total(:,22)==2),10)-gen_total(find(gen_total(:,22)==2),9)*0.3 >= invest_gen_peak_depth(find(gen_total(:,22)==2)) >= 0];
F = [F, gen_total(find(gen_total(:,22)==5),10)-gen_total(find(gen_total(:,22)==5),9)*0.3 >= invest_gen_peak_depth(find(gen_total(:,22)==5)) >= 0];
F = [F, invest_gen_peak_depth(find(gen_total(:,22)==1)) == 0];
% Upper and lower bounds of battery and CAES
F = [F, 500 >= invest_ES_BatteryAndCAES >= 0];
% Upper and lower bounds of DR
F = [F, 0 <= invest_DR1 <= a.bus(:,3)*ratio_DR1];
F = [F, 0 <= invest_DR2 <= a.bus(:,3)*ratio_DR2];

for k = 1:Horizon
    % Ramp rate constraints 
    if k >= 2
        F = [F, PG_G(:,k,:)-PG_G(:,k-1,:) <= onoff(:,k-1,:).*M_RampUp*TimeScale*60+(onoff(:,k,:)-onoff(:,k-1,:)).*M_min_gen2+(1-onoff(:,k,:)).*M_max_gen2];
        F = [F, -PG_G(:,k,:)+PG_G(:,k-1,:) <= onoff(:,k,:).*M_RampDown*TimeScale*60+(-onoff(:,k,:)+onoff(:,k-1,:)).*M_min_gen2+(1-onoff(:,k-1,:)).*M_max_gen2];
		F = [F, PG_CSP(:,k,:)-PG_CSP(:,k-1,:) <= onoff_CSP(:,k-1,:).*M_RampUp_CSP*TimeScale*60+(onoff_CSP(:,k,:)-onoff_CSP(:,k-1,:)).*M_min_CSP2+(1-onoff_CSP(:,k,:)).*M_max_CSP2];
        F = [F, -PG_CSP(:,k,:)+PG_CSP(:,k-1,:) <= onoff_CSP(:,k,:).*M_RampDown_CSP*TimeScale*60+(-onoff_CSP(:,k,:)+onoff_CSP(:,k-1,:)).*M_min_CSP2+(1-onoff_CSP(:,k-1,:)).*M_max_CSP2];
    end
end
% Relation between power and energy of ES
F = [F, PG_ES_BatteryAndCAES(:,1:Horizon-1,:) == E_ES_BatteryAndCAES(:,2:Horizon,:)-E_ES_BatteryAndCAES(:,1:Horizon-1,:)];
F = [F, PG_ES_PumpHydro(:,1:Horizon-1,:) == E_ES_PumpHydro(:,2:Horizon,:)-E_ES_PumpHydro(:,1:Horizon-1,:)];

for sce = 1:N_scenario
    for k = 1:Horizon
        for i = 1:N
            F = [F, PG_G_bus(i,k,sce)+PG_RE_bus(i,k,sce)+PG_CSP_bus(i,k,sce)-(PLoad_day(i,k,sce)-PC_Load(i,k,sce)-PG_DR1(i,k,sce)-PG_DR2(i,k,sce))+PG_ES_BatteryAndCAES_bus(i,k,sce)+PG_ES_PumpHydro_bus(i,k,sce) == -sum(Branch(find(a.branch(:,2)==i),k,sce))-sum(Branch(size(a.branch,1)+find(a.candi_AC(:,2)==i),k,sce))-sum(Branch(size(a.branch,1)+size(a.candi_AC,1)+find(a.candi_DC(:,2)==i),k,sce))+sum(Branch(find(a.branch(:,1)==i),k,sce))+sum(Branch(size(a.branch,1)+find(a.candi_AC(:,1)==i),k,sce))+sum(Branch(size(a.branch,1)+size(a.candi_AC,1)+find(a.candi_DC(:,1)==i),k,sce))];%+PG_CSP_bus(i,k,sce)
        end
    end
end
for sce = 1:N_scenario
    F = [F, PG_G_bus(:,:,sce) == M_bus_G*M_gen_gentotal*PG_G(:,:,sce)];
    F = [F, PG_RE_bus(:,:,sce) == M_bus_RE*PG_RE(:,:,sce)];
	F = [F, PG_CSP_bus(:,:,sce) == M_bus_CSP*M_CSP_CSPtotal*PG_CSP(:,:,sce)];
    F = [F, PG_ES_BatteryAndCAES_bus(:,:,sce) == M_bus_BatteryAndCAES*PG_ES_BatteryAndCAES(:,:,sce)];
    F = [F, PG_ES_PumpHydro_bus(:,:,sce) == M_bus_PumpHydro*PG_ES_PumpHydro(:,:,sce)];
    % Uppwer and lower bounds of power output for hydro, thermal and CSP units
    F = [F, onoff(find(gen_total(:,22)~=1),:,sce).*(M_min_gen_thermal(:,:,sce)-invest_gen_peak_depth(find(gen_total(:,22)~=1))*ones(1,Horizon)) <= PG_G(find(gen_total(:,22)~=1),:,sce) <= onoff(find(gen_total(:,22)~=1),:,sce).*(M_max_gen_thermal(:,:,sce))];
    F = [F, onoff_CSP(:,:,sce).*(M_min_gen_CSP(:,:,sce)) <= PG_CSP(:,:,sce) <= onoff_CSP(:,:,sce).*(M_max_gen_CSP(:,:,sce))];
    F = [F, onoff(find(gen_total(:,22)==1),:,sce)*PHydro_min(day(sce)).*(gen_total(find(gen_total(:,22)==1),9)*ones(1,Horizon)) <= PG_G(find(gen_total(:,22)==1),:,sce) <= onoff(find(gen_total(:,22)==1),:,sce)*PHydro_max(day(sce)).*(gen_total(find(gen_total(:,22)==1),9)*ones(1,Horizon))];
    % Power flow of existing branches
    for line = 1:size(a.branch,1)
        F = [F, Branch(line,:,sce)/100 == (theta(a.branch(line,1),:,sce)-theta(a.branch(line,2),:,sce))/a.branch(line,4)];% 可以简化
        F = [F, -a.branch(line,6)/100*ones(1,Horizon) <= Branch(line,:,sce)/100 <= a.branch(line,6)/100*ones(1,Horizon)];
    end
    % Power flow of candidate AC branches
    for line = size(a.branch,1)+1:size(a.branch,1)+size(a.candi_AC,1)
        F = [F, -M*(1-invest_line_AC(line-size(a.branch,1)))*ones(1,Horizon) <= Branch(line,:,sce)/100-(theta(a.candi_AC(line-size(a.branch,1),1),:,sce)-theta(a.candi_AC(line-size(a.branch,1),2),:,sce))/a.candi_AC(line-size(a.branch,1),4) <= M*(1-invest_line_AC(line-size(a.branch,1)))*ones(1,Horizon) ];
        F = [F, -a.candi_AC(line-size(a.branch,1),6)*invest_line_AC(line-size(a.branch,1))*ones(1,Horizon) <= Branch(line,:,sce) <= a.candi_AC(line-size(a.branch,1),6)*invest_line_AC(line-size(a.branch,1))*ones(1,Horizon)];
    end
    % Power flow of candidate DC branches
    for line = size(a.branch,1)+size(a.candi_AC,1)+1:size(a.branch,1)+size(a.candi_AC,1)+size(a.candi_DC,1)
        F = [F, -a.candi_DC(line-size(a.branch,1)-size(a.candi_AC,1),3)*invest_line_DC(line-size(a.branch,1)-size(a.candi_AC,1))*ones(1,Horizon) <= Branch(line,:,sce) <= a.candi_DC(line-size(a.branch,1)-size(a.candi_AC,1),3)*invest_line_DC(line-size(a.branch,1)-size(a.candi_AC,1))*ones(1,Horizon)];
    end
    % Energy constraints for hydro units
    for Area = 1:5
        F = [F, sum(PG_G((intersect(find(gen_total(:,22)==1),find(gen_total(:,31)==Area))),1:Horizon,sce),2) <= gen_total(intersect(find(gen_total(:,22)==1),find(gen_total(:,31)==Area)),9)*EHydro_day(Area,day(sce))/24*Horizon];
    end
    % ES power output constraints 
    F = [F, -invest_ES_BatteryAndCAES*ones(1,Horizon) <= PG_ES_BatteryAndCAES(:,:,sce) <= invest_ES_BatteryAndCAES*ones(1,Horizon)];
    F = [F, -a.HydroPump(:,3).*invest_ES_PumpHydro*ones(1,Horizon) <= PG_ES_PumpHydro(:,:,sce) <= a.HydroPump(:,3).*invest_ES_PumpHydro*ones(1,Horizon)];
    % ES energy constraints 
    F = [F, 0 <= E_ES_BatteryAndCAES(:,:,sce) <= (invest_ES_BatteryAndCAES.*a.Storage(:,12))*ones(1,Horizon)];% 电量是容量的XX倍
    F = [F, 0 <= E_ES_PumpHydro(:,:,sce) <= (invest_ES_PumpHydro.*a.HydroPump(:,12)).*a.HydroPump(:,3)*ones(1,Horizon)];% 电量是容量的XX倍  
    F = [F, E_ES_BatteryAndCAES(:,1,sce) == (invest_ES_BatteryAndCAES.*a.Storage(:,12))*E_ES0];%初始电量状态
    F = [F, E_ES_PumpHydro(:,1,sce) == (invest_ES_PumpHydro.*a.HydroPump(:,12)).*a.HydroPump(:,3)*E_ES0];%初始电量状态  
    % Load shifting DR constraints 
    F = [F, -invest_DR1*ones(1,Horizon) <= PG_DR1(:,:,sce) <= invest_DR1*ones(1,Horizon)];
    F = [F, PG_DR1(:,:,sce) <= abs_PG_DR1(:,:,sce)];
    F = [F, -PG_DR1(:,:,sce) <= abs_PG_DR1(:,:,sce)];
    % Load shedding DR constraints 
    F = [F, invest_DR2*ones(1,Horizon) >= PG_DR2(:,:,sce) >= 0];
	% CSP internal constraints 
    F = [F, PG_CSP(:,:,sce)/efficiency_ele_ther+Et_TES(:,1:Horizon,sce)-[Capacity_TES_CSP*CSP_total(:,9)*ratio_TES_t0,Et_TES(:,1:Horizon-1,sce)] <= PCSP_day(:,:,sce)];%热功率平衡
    F = [F, 0 <= Et_TES(:,:,sce) <= Capacity_TES_CSP*CSP_total(:,9)*ones(1,Horizon)];%储热能量约束
    F = [F, Et_TES(:,Horizon,sce) == Capacity_TES_CSP*CSP_total(:,9)*ratio_TES_t0];%储热初始能量,储热始末状态一致
    
end
% Wind and PV output constraints
F = [F, 0 <= PG_RE <= PRE_day];
% Load shedding constraints
F = [F, 0 <= PC_Load <= PLoad_day0];
% Reserve constraints
F = [F, sum(onoff.*(M_max_gen)-PG_G,1) >= sum(CoeffReseve_load*PLoad_day0,1)+sum(CoeffReserve_VRE*PG_RE,1)];
% ES initial and final states constraints 
F = [F, E_ES_BatteryAndCAES(:,1,:) == E_ES_BatteryAndCAES(:,Horizon,:)];
F = [F, E_ES_PumpHydro(:,1,:) == E_ES_PumpHydro(:,Horizon,:)];
% Initial and final states constraints Load shifting DR
F = [F, sum(PG_DR1,2) == zeros(N,1,N_scenario)];

obj = (sum(sum(sum(M_CostGeneration.*PG_G.*M_days_gentotal)))+theta_Load*sum(sum(sum(PC_Load.*M_days_bus)))+theta_RE*(sum(sum(sum((PRE_day-PG_RE).*M_days_RE)))))/sum(Cluster_day_number(1:N_scenario,2))*365;% generation cost
obj = obj+(sum(invest_line_AC.*a.candi_AC(:,14))+sum(invest_line_DC.*a.candi_DC(:,4)))*txl;% branch investment cost
obj = obj+sum(invest_gen_peak_depth(find(gen_total(:,22)==2)))*cost_thermal_coal_peak_depth*txl;% retrofit cost of coal-fired thermal units
obj = obj+sum(invest_gen_peak_depth(find(gen_total(:,22)==5)))*cost_thermal_gas_peak_depth*txl;% retrofit cost of gas-fired thermal units
obj = obj+sum(invest_ES_BatteryAndCAES(find(a.Storage(:,2)==2),1).*a.Storage(find(a.Storage(:,2)==2),10).*a.Storage(find(a.Storage(:,2)==2),12)*txl_ES_Battery);% battery investment cost
obj = obj+sum(invest_ES_BatteryAndCAES(find(a.Storage(:,2)==3),1).*a.Storage(find(a.Storage(:,2)==3),10).*a.Storage(find(a.Storage(:,2)==3),12)*txl_ES_CAES);% CAES investment cost
obj = obj+sum(invest_ES_PumpHydro.*a.HydroPump(:,10).*a.HydroPump(:,12).*a.HydroPump(:,3)*txl_ES_HydroPump);% Pumped hydro investment cost
obj = obj+sum(invest_DR1)*cost_DR1_Capacity*txl+sum(sum(sum(abs_PG_DR1.*M_days_bus)))/2*cost_DR1_Energy/sum(Cluster_day_number(1:N_scenario,2))*365;% investment cost and compensation cost of load shifting DR
obj = obj+sum(invest_DR2)*cost_DR2_Capacity*txl+sum(sum(sum(PG_DR2.*M_days_bus)))*cost_DR2_Energy/sum(Cluster_day_number(1:N_scenario,2))*365;% investment cost and compensation cost of load shedding DR
ops = sdpsettings('solver','gurobi','verbose',0,'gurobi.MIPGap',5e-3);
% Run optimization
ans = optimize(F,obj,ops)

% Save results
z_obj = value(obj);
z_cost_toal = value(obj)+sum(a.gen(:,29).*a.gen(:,9))+sum(a.RE(:,7).*a.RE(:,11));% total cost
z_cost_ope = value((sum(sum(sum(M_CostGeneration.*PG_G.*M_days_gentotal)))+theta_Load*sum(sum(sum(PC_Load.*M_days_bus)))+theta_RE*(sum(sum(sum((PRE_day-PG_RE).*M_days_RE)))))/sum(Cluster_day_number(1:N_scenario,2))*365);% generation cost
z_cost_invest_line_AC = value(sum(invest_line_AC.*a.candi_AC(:,14)))*txl;% AC branch investment cost
z_cost_invest_line_DC = value(sum(invest_line_DC.*a.candi_DC(:,4)))*txl;% DC branch investment cost
z_cost_invest_peak_coal = value(sum(invest_gen_peak_depth(find(gen_total(:,22)==2)))*cost_thermal_coal_peak_depth*txl);% retrofit cost of coal-fired thermal units
z_cost_invest_peak_gas = value(sum(invest_gen_peak_depth(find(gen_total(:,22)==5)))*cost_thermal_gas_peak_depth*txl);% retrofit cost of gas-fired thermal units
z_cost_invest_ES_Battery = value(sum(invest_ES_BatteryAndCAES(find(a.Storage(:,2)==2),1).*a.Storage(find(a.Storage(:,2)==2),10).*a.Storage(find(a.Storage(:,2)==2),12)*txl_ES_Battery));% battery investment cost
z_cost_invest_ES_CAES = value(sum(invest_ES_BatteryAndCAES(find(a.Storage(:,2)==3),1).*a.Storage(find(a.Storage(:,2)==3),10).*a.Storage(find(a.Storage(:,2)==3),12)*txl_ES_CAES));% CAES investment cost
z_cost_invest_ES_PumpHydro = value(sum(invest_ES_PumpHydro.*a.HydroPump(:,10).*a.HydroPump(:,12).*a.HydroPump(:,3)*txl_ES_HydroPump));% Pumped hydro investment cost
z_cost_invest_DR1 = value(sum(invest_DR1)*cost_DR1_Capacity*txl);% investment cost and of load shifting DR
z_cost_invest_DR1_compen = value(sum(sum(sum(abs_PG_DR1.*M_days_bus)))/2*cost_DR1_Energy/sum(Cluster_day_number(1:N_scenario,2))*365);% compensation cost of load shifting DR
z_cost_invest_DR2 = value(sum(invest_DR2)*cost_DR2_Capacity*txl);% investment cost of load shedding DR
z_cost_invest_DR2_compen = value(sum(sum(sum(PG_DR2.*M_days_bus)))*cost_DR2_Energy/sum(Cluster_day_number(1:N_scenario,2))*365);% compensation cost of load shedding DR
z_cost_invest_sum = z_cost_invest_line_AC+z_cost_invest_line_DC+z_cost_invest_peak_coal+z_cost_invest_peak_gas+z_cost_invest_ES_Battery+z_cost_invest_ES_CAES+z_cost_invest_ES_PumpHydro+z_cost_invest_DR1+z_cost_invest_DR2; %total annual investment cost
z_PG_G = value(PG_G);
z_PG_RE = value(PG_RE);
z_PC_Load = value(PC_Load);
z_PG_ES = value(PG_ES_BatteryAndCAES);
z_E_ES = value(E_ES_BatteryAndCAES);
z_PG_DR1 = value(PG_DR1);
z_PG_DR2 = value(PG_DR2);
z_PG_CSP = value(PG_CSP);
z_RCR = (sum(sum(sum(PRE_day)))-sum(sum(sum(z_PG_RE))))/sum(sum(sum(PRE_day)));
z_PG_Thermal = z_PG_G(find(gen_total(:,22)~=1),:,:);
z_PG_Hydro = z_PG_G(find(gen_total(:,22)==1),:,:);
z_PG_Wind = z_PG_RE(find(a.RE(:,3)==4),:,:);
z_PG_Solar = z_PG_RE(find(a.RE(:,3)==3),:,:);
z_Branch = value(Branch);

z_invest_line_AC = value(invest_line_AC);
z_invest_line_DC = value(invest_line_DC);
z_invest_gen_peak_depth = value(invest_gen_peak_depth);
z_invest_sum_peak = sum(z_invest_gen_peak_depth);
z_invest_ES_BatteryAndCAES = value(invest_ES_BatteryAndCAES);
z_invest_ES_PumpHydro = value(invest_ES_PumpHydro);
z_invest_DR1 = value(invest_DR1);
z_invest_DR2 = value(invest_DR2);

% Investment capacity 
z_Capacity_invest_gen_peak_depth_Coal = sum(z_invest_gen_peak_depth(find(gen_total(:,22)==2)));
z_Capacity_invest_gen_peak_depth_Gas = sum(z_invest_gen_peak_depth(find(gen_total(:,22)==5)));
z_Capacity_invest_line_AC_CollectionGrid = sum(z_invest_line_AC(find(a.candi_AC(:,15)==2)).*a.candi_AC(find(a.candi_AC(:,15)==2),6));
z_Capacity_invest_line_AC_LocalTransGrid = sum(z_invest_line_AC(find(a.candi_AC(:,15)==3)).*a.candi_AC(find(a.candi_AC(:,15)==3),6));
z_Capacity_invest_line_AC_Tieline = sum(z_invest_line_AC(find(a.candi_AC(:,15)==1)).*a.candi_AC(find(a.candi_AC(:,15)==1),6));
z_Capacity_invest_line_DC_CollectionGrid = sum(z_invest_line_DC(find(a.candi_DC(:,5)==2)).*a.candi_DC(find(a.candi_DC(:,5)==2),3));
z_Capacity_invest_line_DC_LocalTransGrid = sum(z_invest_line_DC(find(a.candi_DC(:,5)==3)).*a.candi_DC(find(a.candi_DC(:,5)==3),3));
z_Capacity_invest_line_DC_Tieline = sum(z_invest_line_DC(find(a.candi_DC(:,5)==1)).*a.candi_DC(find(a.candi_DC(:,5)==1),3));
z_Capacity_invest_ES_Battery = sum(z_invest_ES_BatteryAndCAES(find(a.Storage(:,2)==2)));
z_Capacity_invest_ES_CAES = sum(z_invest_ES_BatteryAndCAES(find(a.Storage(:,2)==3)));
z_Capacity_invest_ES_HydroPump = sum(z_invest_ES_PumpHydro.*a.HydroPump(:,3));
z_Capacity_invest_DR1 = sum(z_invest_DR1);
z_Capacity_invest_DR2 = sum(z_invest_DR2);

% Investment cost
z_cost_invest_line_AC_CollectionGrid = sum(z_invest_line_AC(find(a.candi_AC(:,15)==2)).*a.candi_AC(find(a.candi_AC(:,15)==2),14))*txl;
z_cost_invest_line_AC_LocalTransGrid = sum(z_invest_line_AC(find(a.candi_AC(:,15)==3)).*a.candi_AC(find(a.candi_AC(:,15)==3),14))*txl;
z_cost_invest_line_AC_Tieline = sum(z_invest_line_AC(find(a.candi_AC(:,15)==1)).*a.candi_AC(find(a.candi_AC(:,15)==1),14))*txl;
z_cost_invest_line_DC_CollectionGrid = sum(z_invest_line_DC(find(a.candi_DC(:,5)==2)).*a.candi_DC(find(a.candi_DC(:,5)==2),4))*txl;
z_cost_invest_line_DC_LocalTransGrid = sum(z_invest_line_DC(find(a.candi_DC(:,5)==3)).*a.candi_DC(find(a.candi_DC(:,5)==3),4))*txl;
z_cost_invest_line_DC_Tieline = sum(z_invest_line_DC(find(a.candi_DC(:,5)==1)).*a.candi_DC(find(a.candi_DC(:,5)==1),4))*txl;

str = 'Results_Flexibility_Resource_Planning';
save(str)