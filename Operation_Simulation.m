% The operation simualtion codes for the FTS-213 test system.
% Author: Hao Li, email: h-l16@mails.tsinghua.edu.cn
% updated: 2020_March_02

% Matpower, Yalmip and Gurobi packages are required to run the codes.
% Matpower package: R. D. Zimmerman, C. E. Murillo-Sanchez (2019). MATPOWER (Version 7.0) [Software]. Available: https://matpower.org
% Yalmip package: https://yalmip.github.io/
% Gurobi package: https://www.gurobi.com/

% Hydrological condition is set as dry year.

load('data_20200214.mat');
load('month_week');
load('DCtrans');
load('TieLine2');
load('Cluster_day_number');

%% Parameter setting
Horizon = 24;
TimeScale = 1;% 1h as the time scale
CoeffReseve_load = 0.03;% reserve coefficient for load
CoeffReserve_VRE = 0.05;% reserve coefficient for VRE
efficiency_TES = 0.98;% efficiency of CSP TES(Thermal energy storage)
efficiency_ele_ther = 0.415;% efficiency of electric power output over input thermal power of CSP power block
Capacity_TES_CSP = 2.4; % TES capacity relative to the capacity of power block of CSP
ratio_TES_t0 = 0.5;% initial energy ratio of TES
TES_hour = 10;% TES-hour
theta_Load = 3/7*1000;% $/MWh % load shedding penalty 

% read the dataset of FTS-213
a = case213_20200302;

N = length(a.bus(:,1));
Nbr = length(a.branch(:,1));
Ngen = length(a.gen(:,1));
NRE = length(a.RE(:,1));
NCSP = length(a.CSP(:,1));

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

GSDF = makePTDF(a);

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
    mean_PLoad(:,Area) = mean(PLoad_all{Area},2);
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

% ------ Inter-regional power exchange ------
TieLine_section = [2,1;2,4;2,3;5,2];
section = zeros(4,8760);
for s = 1:4
    for i = 1:365
        if i <= month_week(1,3)
            if TieLine2(s,1) > 0
                section(s,(i-1)*24+1:i*24) = TieLine2(s,1);
            else
                section(s,(i-1)*24+1:i*24) = TieLine2(s,1);
            end
        end
        for j = 2:12
            if (i>month_week(j-1,3))&&(i<=month_week(j,3))
                if TieLine2(s,j) > 0
                    section(s,(i-1)*24+1:i*24) = TieLine2(s,j);
                else
                    section(s,(i-1)*24+1:i*24) = TieLine2(s,j);
                end
                break;
            end
        end
    end
end
% ------ Inter-regional power exchange ------

N_day = size(Cluster_day_number,1);
z_obj = zeros(N_day,1);
z_PG_G = cell(N_day,1);
z_PG_RE = cell(N_day,1);
z_PC_Load = cell(N_day,1);
z_RCR = zeros(N_day,1);
z_PG_Thermal_Coal =cell(N_day,1);
z_PG_Hydro = cell(N_day,1);
z_PG_Wind = cell(N_day,1);
z_PG_Solar = cell(N_day,1);
z_sum_PC_load =cell(N_day,1);
z_Branch = cell(N_day,1);
z_op_result = zeros(N_day,1);
z_Hydro_CR_ = zeros(N_day,1);
z_onoff_gen = cell(N_day,1);
z_section = cell(N_day,1);

Count_day = 0;
%% Typical days
for day = (Cluster_day_number(:,1))'
    day
    Count_day = Count_day+1
    
    % ------ Load ------
    % --- Local load ---
    PLoad_day0 = zeros(N,Horizon);
    for row = 1:N
        if a.bus(row,11) == 1
            PLoad_day0(row,:) = a.bus(row,3)*PLoad_all{6}((day-1)*Horizon+1:day*Horizon);
        elseif a.bus(row,11) == 2
            PLoad_day0(row,:) = a.bus(row,3)*PLoad_all{7}((day-1)*Horizon+1:day*Horizon);
        elseif a.bus(row,11) == 3
            PLoad_day0(row,:) = a.bus(row,3)*PLoad_all{8}((day-1)*Horizon+1:day*Horizon);
        elseif a.bus(row,11) == 4
            PLoad_day0(row,:) = a.bus(row,3)*PLoad_all{9}((day-1)*Horizon+1:day*Horizon);
        else
            PLoad_day0(row,:) = a.bus(row,3)*PLoad_all{10}((day-1)*Horizon+1:day*Horizon);
        end
    end
    PLoad_day = PLoad_day0;
    % --- Local load ---
    % --- External load ---
    if day <= 31
        for b = 1:size(DCtrans,1)
            PLoad_day(find(a.bus(:,1)==DCtrans(b,1)),1:Horizon) = PLoad_day(find(a.bus(:,1)==DCtrans(b,1)),1:Horizon)+DCtrans(b,2)*0.6;
        end
    else
        for month = 2:12
            if day>month_week(month-1,3) && day<=month_week(month,3)
                for b = 1:size(DCtrans,1)
                    PLoad_day(find(a.bus(:,1)==DCtrans(b,1)),1:Horizon) = PLoad_day(find(a.bus(:,1)==DCtrans(b,1)),1:Horizon)+DCtrans(b,month+1)*0.6;
                end
            end
        end
    end
    % --- External load ---
    % ------ Load ------
    
    % ------ Wind, PV, CSP ------
    PRE_hour = zeros(NRE,Horizon);
    for row = 1:NRE
        if a.RE(row,3) == 3
            if a.RE(row,2) == 1
                PRE_hour(row,:) = a.RE(row,7)*PSolar_all{6}((day-1)*Horizon+1:day*Horizon);
            elseif a.RE(row,2) == 2
                PRE_hour(row,:) = a.RE(row,7)*PSolar_all{7}((day-1)*Horizon+1:day*Horizon);
            elseif a.RE(row,2) == 3
                PRE_hour(row,:) = a.RE(row,7)*PSolar_all{8}((day-1)*Horizon+1:day*Horizon);
            elseif a.RE(row,2) == 4
                PRE_hour(row,:) = a.RE(row,7)*PSolar_all{9}((day-1)*Horizon+1:day*Horizon);
            else
                PRE_hour(row,:) = a.RE(row,7)*PSolar_all{10}((day-1)*Horizon+1:day*Horizon);
            end
        else
            if a.RE(row,2) == 1
                PRE_hour(row,:) = a.RE(row,7)*PWind_all{6}((day-1)*Horizon+1:day*Horizon);
            elseif a.RE(row,2) == 2
                PRE_hour(row,:) = a.RE(row,7)*PWind_all{7}((day-1)*Horizon+1:day*Horizon);
            elseif a.RE(row,2) == 3
                PRE_hour(row,:) = a.RE(row,7)*PWind_all{8}((day-1)*Horizon+1:day*Horizon);
            elseif a.RE(row,2) == 4
                PRE_hour(row,:) = a.RE(row,7)*PWind_all{9}((day-1)*Horizon+1:day*Horizon);
            else
                PRE_hour(row,:) = a.RE(row,7)*PWind_all{10}((day-1)*Horizon+1:day*Horizon);
            end
        end
    end
    PtCSP_hour = zeros(NCSP,Horizon);
    for row = 1:NCSP
        if a.CSP(row,31) == 1
            PtCSP_hour(row,:) = a.CSP(row,9)*PCSP_all{6}((day-1)*Horizon+1:day*Horizon);
        elseif a.CSP(row,31) == 2
            PtCSP_hour(row,:) = a.CSP(row,9)*PCSP_all{7}((day-1)*Horizon+1:day*Horizon);
        elseif a.CSP(row,31) == 3
            PtCSP_hour(row,:) = a.CSP(row,9)*PCSP_all{8}((day-1)*Horizon+1:day*Horizon);
        elseif a.CSP(row,31) == 4
            PtCSP_hour(row,:) = a.CSP(row,9)*PCSP_all{9}((day-1)*Horizon+1:day*Horizon);
        else
            PtCSP_hour(row,:) = a.CSP(row,9)*PCSP_all{10}((day-1)*Horizon+1:day*Horizon);
        end
    end
    % ------ Wind, PV, CSP ------
    
    %% Optimization
    PG_G = sdpvar(size(gen_total,1),Horizon,'full');% Hydro and thermal power
    PG_RE = sdpvar(NRE,Horizon,'full');% Wind and PV power
    PG_CSP = sdpvar(size(CSP_total,1),Horizon,'full');% CSP power
    PC_Load = sdpvar(N,Horizon,'full');% Load shedding
    onoff_gen = binvar(size(gen_total,1),Horizon,'full');% on/off status of hydro and thermal units
    onoff_CSP = binvar(size(CSP_total,1),Horizon,'full');% on/off status of CSP units
    Branch = sdpvar(size(a.branch,1),Horizon,'full');% power flow
    Cost_StartUp = sdpvar(size(gen_total,1),Horizon-1,'full');% startup cost
    Pt_TES_charge = sdpvar(size(CSP_total,1),Horizon,'full');% input thermal power for TES of CSP
    Pt_TES_discharge= sdpvar(size(CSP_total,1),Horizon,'full');% output thermal power for TES of CSP
    Et_TES = sdpvar(size(CSP_total,1),Horizon,'full');% TES energy of CSP
    F = [];
    for k = 1:Horizon
        % Minimum on/off time constraints for thermal units
        if k >= 2
            for i = 1:size(gen_total,1)
                if (gen_total(i,22)==2)||(gen_total(i,22)==5)
                    for tao = k+1:min(k+gen_total(i,25)-1,Horizon)
                        F = [F, onoff_gen(i,k)-onoff_gen(i,k-1) <= onoff_gen(i,tao)];
                    end
                    for tao = k+1:min(k+gen_total(i,26)-1,Horizon)
                        F = [F, onoff_gen(i,k-1)-onoff_gen(i,k) <= 1-onoff_gen(i,tao)];
                    end
                end
            end
        end
        % Minimum on/off time constraints for CSP units
        if k >= 2
            for i = 1:size(CSP_total,1)
                for tao = k+1:min(k+CSP_total(i,25)-1,Horizon)
                    F = [F, onoff_CSP(i,k)-onoff_CSP(i,k-1) <= onoff_CSP(i,tao)];
                end
                for tao = k+1:min(k+CSP_total(i,26)-1,Horizon)
                    F = [F, onoff_CSP(i,k-1)-onoff_CSP(i,k) <= 1-onoff_CSP(i,tao)];
                end
            end
        end
        % Ramp rate constraints
        if k >= 2
            F = [F, PG_G(:,k)-PG_G(:,k-1) <= onoff_gen(:,k-1).*gen_total(:,23)*TimeScale*60+(onoff_gen(:,k)-onoff_gen(:,k-1)).*gen_total(:,10)+(1-onoff_gen(:,k)).*gen_total(:,9)];
            F = [F, -PG_G(:,k)+PG_G(:,k-1) <= onoff_gen(:,k).*gen_total(:,24)*TimeScale*60+(-onoff_gen(:,k)+onoff_gen(:,k-1)).*gen_total(:,10)+(1-onoff_gen(:,k-1)).*gen_total(:,9)];
            F = [F, PG_CSP(:,k)-PG_CSP(:,k-1) <= onoff_CSP(:,k-1).*CSP_total(:,23)*TimeScale*60+(onoff_CSP(:,k)-onoff_CSP(:,k-1)).*CSP_total(:,10)+(1-onoff_CSP(:,k)).*CSP_total(:,9)];
            F = [F, -PG_CSP(:,k)+PG_CSP(:,k-1) <= onoff_CSP(:,k).*CSP_total(:,24)*TimeScale*60+(-onoff_CSP(:,k)+onoff_CSP(:,k-1)).*CSP_total(:,10)+(1-onoff_CSP(:,k-1)).*CSP_total(:,9)];
        end
    end
    % Upper/lower bounds of power output for hydro, thermal and CSP units
    F = [F, onoff_gen(union(find(gen_total(:,22)==2),find(gen_total(:,22)==5)),:).*(gen_total(union(find(gen_total(:,22)==2),find(gen_total(:,22)==5)),10)*ones(1,Horizon)) <= PG_G(union(find(gen_total(:,22)==2),find(gen_total(:,22)==5)),:) <= onoff_gen(union(find(gen_total(:,22)==2),find(gen_total(:,22)==5)),:).*(gen_total(union(find(gen_total(:,22)==2),find(gen_total(:,22)==5)),9)*ones(1,Horizon))];
    F = [F, onoff_CSP.*(CSP_total(:,10)*ones(1,Horizon)) <= PG_CSP <= onoff_CSP.*(CSP_total(:,9)*ones(1,Horizon))];
    F = [F, onoff_gen(find(gen_total(:,22)==1),:)*PHydro_min(day).*(gen_total(find(gen_total(:,22)==1),9)*ones(1,Horizon)) <= PG_G(find(gen_total(:,22)==1),:) <= onoff_gen(find(gen_total(:,22)==1),:)*PHydro_max(day).*(gen_total(find(gen_total(:,22)==1),9)*ones(1,Horizon))];
    % Power balance constraints
    F = [F, sum(PG_G,1)+sum(PG_RE,1)+sum(PG_CSP,1) == sum(PLoad_day-PC_Load,1)];
    % Power flow constraints
    F = [F, Branch == GSDF*(M_bus_G*M_gen_gentotal*PG_G+M_bus_RE*PG_RE+M_bus_CSP*M_CSP_CSPtotal*PG_CSP-(PLoad_day-PC_Load))];
    F = [F, -a.branch(:,6)*ones(1,Horizon) <= GSDF*(M_bus_G*M_gen_gentotal*PG_G+M_bus_RE*PG_RE+M_bus_CSP*M_CSP_CSPtotal*PG_CSP-(PLoad_day-PC_Load)) <= a.branch(:,6)*ones(1,Horizon)];
    % Inter-regional power constraints
    F = [F, -section(1,(day-1)*24+1:(day-1)*24+Horizon) <= sum(Branch(intersect(find(a.branch(:,1)<=Area_bus(1)),intersect(find(a.branch(:,2)>Area_bus(1)),find(a.branch(:,2)<=Area_bus(2)))),:),1) <= -section(1,(day-1)*24+1:(day-1)*24+Horizon)];
    F = [F, section(2,(day-1)*24+1:(day-1)*24+Horizon) <= sum(Branch(intersect(intersect(find(a.branch(:,1)>Area_bus(1)),find(a.branch(:,1)<=Area_bus(2))),intersect(find(a.branch(:,2)>Area_bus(3)),find(a.branch(:,2)<=Area_bus(4)))),:),1) <= section(2,(day-1)*24+1:(day-1)*24+Horizon)];
    F = [F, section(3,(day-1)*24+1:(day-1)*24+Horizon) <= sum(Branch(intersect(intersect(find(a.branch(:,1)>Area_bus(1)),find(a.branch(:,1)<=Area_bus(2))),intersect(find(a.branch(:,2)>Area_bus(2)),find(a.branch(:,2)<=Area_bus(3)))),:),1) <= section(3,(day-1)*24+1:(day-1)*24+Horizon)];
    F = [F, -section(4,(day-1)*24+1:(day-1)*24+Horizon) <= sum(Branch(intersect(intersect(find(a.branch(:,1)>Area_bus(1)),find(a.branch(:,1)<=Area_bus(2))),intersect(find(a.branch(:,2)>Area_bus(4)),find(a.branch(:,2)<=Area_bus(5)))),:),1) <= -section(4,(day-1)*24+1:(day-1)*24+Horizon)];
    % Energy constraints for hydro units
    for Area = 1:5
        F = [F, sum(PG_G((intersect(find(gen_total(:,22)==1),find(gen_total(:,31)==Area))),1:Horizon),2) <= gen_total(intersect(find(gen_total(:,22)==1),find(gen_total(:,31)==Area)),9)*EHydro_day(Area,day)/24*Horizon];
    end
    % Hydro units are set as online
    F = [F, onoff_gen(find(gen_total(:,22)==1),:) == 1];
    % Wind and PV output constraints
    F = [F, 0 <= PG_RE <= PRE_hour];
    % Load shedding constraints
    F = [F, 0 <= PC_Load <= PLoad_day0];
    % Reserve constraints
    F = [F, sum(onoff_gen.*(gen_total(:,9)*ones(1,Horizon))-PG_G,1)+sum(onoff_CSP.*(CSP_total(:,9)*ones(1,Horizon))-PG_CSP,1) >= sum(CoeffReseve_load*PLoad_day0,1)+sum(CoeffReserve_VRE*PG_RE,1)];
    % Startup cost constraints
    F = [F, Cost_StartUp >= (onoff_gen(:,2:Horizon)-onoff_gen(:,1:Horizon-1)).*(gen_total(:,28)*ones(1,Horizon-1))];
    F = [F, Cost_StartUp >= 0];
    % CSP internal constraints
    F = [F, PG_CSP/efficiency_ele_ther+Pt_TES_charge-Pt_TES_discharge <= PtCSP_hour];% Thermal power balance
    F = [F, Et_TES(:,2:Horizon)-Et_TES(:,1:Horizon-1) == Pt_TES_charge(:,1:Horizon-1)*efficiency_TES-Pt_TES_discharge(:,1:Horizon-1)/efficiency_TES];
    F = [F, 0 <= [Pt_TES_charge;Pt_TES_discharge] <= Capacity_TES_CSP*ones(size(CSP_total,1)*2,Horizon)];
    F = [F, 0 <= Et_TES <= Capacity_TES_CSP*CSP_total(:,9)*ones(1,Horizon)];
    F = [F, Et_TES(:,1) == Capacity_TES_CSP*CSP_total(:,9)*ratio_TES_t0];
    F = [F, Et_TES(:,1) == Et_TES(:,Horizon)];
    
    % Objective function
    obj = sum(gen_total(:,30)'*PG_G)+sum(CSP_total(:,30)'*PG_CSP)+sum(sum(Cost_StartUp))+theta_Load*sum(sum(PC_Load));
    % Run optimization
    ops = sdpsettings('solver','gurobi','verbose',0,'gurobi.MIPGap',5e-3,'gurobi.TimeLimit',1.45e5);
    ans = optimize(F,obj,ops)
    
    % Save results
    z_obj(Count_day) = value(obj);
    z_PG_G{Count_day} = value(PG_G);
    z_PG_RE{Count_day} = value(PG_RE);
    z_PG_CSP{Count_day} = value(PG_CSP);
    z_PC_Load{Count_day} = value(PC_Load);
    z_RCR(Count_day) = (sum(sum(PRE_hour))-sum(sum(z_PG_RE{Count_day})))/sum(sum(PRE_hour));
    z_PG_Thermal_Coal{Count_day} = z_PG_G{Count_day}(find(gen_total(:,22)==2),:);
    z_PG_Thermal_Gas{Count_day} = z_PG_G{Count_day}(find(gen_total(:,22)==5),:);
    z_PG_Hydro{Count_day} = z_PG_G{Count_day}(find(gen_total(:,22)==1),:);
    z_PG_Wind{Count_day} = z_PG_RE{Count_day}(find(a.RE(:,3)==4),:);
    z_PG_Solar{Count_day} = z_PG_RE{Count_day}(find(a.RE(:,3)==3),:);
    z_Et_TES = value(Et_TES);
    z_Pt_TES_charge = value(Pt_TES_charge);
    z_Pt_TES_discharge = value(Pt_TES_discharge);
    z_onoff_gen{Count_day} = value(onoff_gen);
    z_onoff_CSP{Count_day} = value(onoff_CSP);
    z_Branch{Count_day} = GSDF*(M_bus_G*M_gen_gentotal*z_PG_G{Count_day}+M_bus_RE*z_PG_RE{Count_day}-(PLoad_day-z_PC_Load{Count_day}));
    z_section{Count_day} = [
        value(sum(z_Branch{Count_day}(intersect(find(a.branch(:,1)<=Area_bus(1)),intersect(find(a.branch(:,2)>Area_bus(1)),find(a.branch(:,2)<=Area_bus(2)))),:),1));
        sum(z_Branch{Count_day}(intersect(intersect(find(a.branch(:,1)>Area_bus(1)),find(a.branch(:,1)<=Area_bus(2))),intersect(find(a.branch(:,2)>Area_bus(3)),find(a.branch(:,2)<=Area_bus(4)))),:),1)
        sum(z_Branch{Count_day}(intersect(intersect(find(a.branch(:,1)>Area_bus(1)),find(a.branch(:,1)<=Area_bus(2))),intersect(find(a.branch(:,2)>Area_bus(2)),find(a.branch(:,2)<=Area_bus(3)))),:),1)
        sum(z_Branch{Count_day}(intersect(intersect(find(a.branch(:,1)>Area_bus(1)),find(a.branch(:,1)<=Area_bus(2))),intersect(find(a.branch(:,2)>Area_bus(4)),find(a.branch(:,2)<=Area_bus(5)))),:),1)];
    
     str = strcat('Result_Operation_Simulation',int2str(Count_day));
     save(str);
end
z_obj_year = z_obj'*Cluster_day_number(:,2); % annual operation cost
z_RCR_year = z_RCR'*Cluster_day_number(:,2)/365; % annual average VRE curtailment rate
Energy_RE = 0;
Energy_all = 0;
for i = 1:20
    Energy_RE = Energy_RE+sum(sum(z_PG_RE{i}))*Cluster_day_number(i,2);
    Energy_all = Energy_all+(sum(sum(z_PG_RE{i}))+sum(sum(z_PG_G{i})))*Cluster_day_number(i,2);
end
z_Energy_rate_RE = Energy_RE/Energy_all; % Energy share of VRE
 save(str);
