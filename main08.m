%%%%%%%%%%%%%% HW3 Part 1 %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%Define general properties of element%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Concrete Parameters
D = 4*12;           %Diameter of Column [in]
L = 24*12;          %height of column [in] (from base to center of mass of inertial block) 
cch = 2;            %concrete clear cover
ds = 0.625;         %hoop diameter [in]
Adb = 0.62;         % area of double hoops [in^2
dc = D-2*cch-ds     %Confined diameter of concrete [in] (from hoop center to hoop center) 
fyh = 60.4;         %Yield strength of hoop [ksi]
s = 6;              %center to center spacing of hoops [in]
sp = s - 2*ds;      %clear spacing of hoops [in] (subtract 2 bar diameters because the hoops are doubled)
Ec = 3250;          %Modulus of elasticity of concrete [ksi]

fpc = 6;           %%%%-- MOHAMED CHANGED THIS --%%%% %Compressive strength of the column [ksi]
ne = 1.6;           % non-linearity coefficient of pre-peak stress-strain curve
epc = -ne*fpc/Ec; %%%%-- MOHAMED CHANGED THIS --%%%%  %strain at the unconfined compressive strength 
% % cf (negative)   %%-- AUSTIN, WHAT IS THIS? --%%
Gfc = 4.2*((fpc*1000)^0.5)/1000
d_ppku = -2*(Gfc/fpc)
Gl = 20   %Gauge length
ecu = epc+(fpc/Ec)+(d_ppku/Gl);      %Crushing strength of unconfined concrete

%k1 = ????          %unsure what this variable is for

Kc = 4.1;           %%%%-- MOHAMED CHANGED THIS --%%%%
Ke = 1 - sp/dc;     %%%%-- MOHAMED CHANGED THIS --%%%% % Eq (2) on Pg 18 of Homework 3 r2.3.pdf
f_pl = 0.03*fpc;   %%%%-- MOHAMED CHANGED THIS --%%%%  %passive confining stress provided by the lateral reinforcement 
d_ppkccu = -2*(1+40*(Ke*f_pl/fpc))*Gfc/fpc
fpcc = fpc + Ke*Kc*f_pl; %6 %%-- MOHAMED CHANGED THIS --%%     %Compressive strength of confined concrete [ksi]
epcc = epc-(1/20)*Ke*f_pl/fpc  %%-0.0027;     %%-- MOHAMED CHANGED THIS --%% 
eccu = epc+(fpc/Ec)+(d_ppkccu/Gl);     %Ultimate strain of confined concrete
fcr = 0.440;        %Tensile strength of the concrete in [ksi]
ecr = 0.00014;      %Tensile strain corresponding to f_cr
Esec = fpc/epc;       %secant modulus unconfined conc
rc = 1.5 / (1 - (1/ne) ); 
Esecc = -fpcc/epcc;    %secant modulus confined conc 
rcc = 1 / (1 - (Esecc/Ec)); 

       % subject to change?
lambda_conc = 8/Gl

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Longitudinal Steel Parameters %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Es = 29000;                 %Young's Modulus of Steel [ksi]
fy = 69;                    %[ksi]
ey = fy/Es;                 %Yield Strain
esh = 0.01;                 %Strain at onset of strain hardening
esu = 0.107;                %Longitudinal Reinforcement Strain at tensile strength
fsu = 101.8;                %Longitudinal Reinforcement tensile strength magnitude
P = 3  ;                    %strain hardening power term
dl = 1.375;                 %diameter of longitudinal bars [in]
lp1 = max(0.5*D, 0.08*L);   %Length of plastic hinge 1 [in]
lp2 = 6*dl;               %Plastic hinge length 2 [in]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Define the Loads %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w_conc = 0.150/(12^3);          %Concrete unit weight [kip/in^3]
P_col = -pi*(D/2)^2*L*w_conc;   %Self-weight of column [kip] (Negative for compression)
P_struct = -522;                %weight of super_structure [kip] (Negative for compression)
P_total= P_col +P_struct;       %total weight [kip]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Define Geometry%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%Concrete%%%%%%%%
nc = 20; % number of concrete fibers
t_conc = D/nc;  %width of concrete fibers [in]

[yc,Ac,Acc] = areas_circle(D,dc,t_conc,nc);  
% yc is the distance from the bottom of the section to the centroid of each
% fiber (vector)
% Ac is the area of unconfined concrete of each fiber (vector)
% Acc is the area of confined concrete of each fiber (vector)

%%%%%%%%%Steel%%%%%%%%%%%
n_bars = 18; %actual rebar count
n_steel = 5; %number of discretized steel fibers
ysi = [3.4; 8.8; 24; 39.2; 44.7]; %centroid of each discretized fiber [in]
Asi = [4.68; 6.24; 6.24; 6.24; 4.68]; %Area of each discretized fiber [in^2]


%%%%%%%%%%Properties of Transformed Section %%%%%%%%%%
n=Es/Ec;
At= (pi*(D/2)^2-n_bars*pi*(dl/2)^2)+ n_bars*pi*(dl/2)^2*n;  %Area of Transformed Section [in]
yt= D/2;  %Transformed Controid [in]

It=0;  %initialie transformed moment of inertia [in^4] (solved for below in for loop)
I_conc = 1/4*pi*(D/2)^4;
for i=1:length(ysi)                       
    It = It + Asi(i)*(yt-ysi(i))^2*n +n*n_bars*pi*(dl/2)^4 - n_bars*pi*(dl/2)^4;  %Account for transformed steel. Subtract moment of inertia of steel from conc bc we took entire amoment of ienrtia of column
end           
It = I_conc + It
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Moment Curvature %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Inital Strain and Curvature
e(1) = P_total/(Ec*At);  %Initial strain at centroid
M(1) = 0;   %Initial moment taken from centroid(1) = 0   
phi(1) = M(1)/(Ec*It);  %initial curvature [1/in]

phi_delta = 2*ey/(D*20);  %Curvature change

tol = 0.01; %Arbitrarily chosen can change later
count = 2; %Initalize counter Starts at two bc everything is initialied above
es_max = 0; %Initialize Max steel strain 
ec_min = 0; %Initialize min concrete strain

es_lim = esu;     %Initialize steel strain limit
ec_lim = eccu;   %Initialize conc strain limit;  $$$check here if limit is weird

%%%% Debuging Initializations
%error_plot(1) = 1;
%count_vec(1)=1;
%count_test = 1;

%Initialize vairables for flagging points
fcr_flag = 0
fcr_flag_triggered = 0
yield_flag=0
yield_flag_triggered = 0
ACI_flag=0
ACI_flag_triggered = 0
nom_moment_flag=0
nom_moment_flag_triggered = 0
spall_flag=0
spall_flag_triggered = 0

M_max= M(1)
phi_max = 0
while es_max <= es_lim && ec_min >=ec_lim
    de = -0.5*phi_delta;        %arbitrarily chosen strain change 
    e(count) = e(count-1)+de; %define strain at current step
       
    
    phi(count) = phi(count-1) + phi_delta;  %Define curvature at current step
    
    error = 1;  
    while abs(error)>tol
        Fc = zeros(nc,1); %Initialie concrete force vector for each layer
        Fs = zeros(length(Asi),1);  %Initialie steel force vector for each layer
        test_count = 0;
        ec_j = zeros(length(yc),1);
        es_j = zeros(length(ysi),1);
        
        for j=1:nc   %Loop through concrete layers
              %test_count = test_count+1; %debugging
              F_layer=0;
              sigma_c= 0; %Confined stress initialiation
              sigma_u = 0; %unconfined stress initialiation
              sigma_t = 0;  %tensile stress initialiation
              ec_j(j) = phi(count)*(yt-yc(j))+e(count);  %strain at conc layer j
              
              if ec_j(j)<0 %conc in compression
                  sigma_c = stress_confined_conc(ec_j(j), epcc, fpcc, lambda_conc,eccu,  rcc);
                  sigma_u = stress_unconfined_conc(ec_j(j), epc, ecu, fpc, ne, lambda_conc, rc);
              else
                  sigma_t =tens_conc(Ec,fcr,ec_j(j),ecr);
              end
              if fcr_flag_triggered==0  %Flag function to find cracking moment and curvature to plot later 
                 if sigma_t >= fcr
                     fcr_flag=1
                     fcr_flag_triggered = 1
                 end
              end
               if yield_flag_triggered==0  %Flag function to find yield moment and curvature to plot later 
                 if ec_j(j) <= epc
                     yield_flag=1
                     yield_flag_triggered = 1
                     steel_yield = 1 %Initialize to know later on if it was steel or concrete that yeilded first 
                 end
               end
               if ACI_flag_triggered==0  %Flag function to find ACI limit (e=-.003) moment and curvature to plot later 
                 if ec_j(j) <= -0.003
                     ACI_flag=1
                     ACI_flag_triggered = 1
                 end
               end
               if nom_moment_flag_triggered==0  %Flag function to find nominal moment limit (e=-.004) moment and curvature to plot later 
                 if ec_j(j) <= -0.004
                     nom_moment_flag=1
                     nom_moment_flag_triggered = 1
                 end
               end
               if spall_flag_triggered==0  %Flag function to find spalling moment and curvature to plot later 
                 if ec_j(j) <= ecu
                     spall_flag=1
                     spall_flag_triggered = 1
                 end
               end
              Fc(j) = sigma_c*Acc(j) + sigma_u*Ac(j)+sigma_t*(Acc(j) +Ac(j));
        end
        for j=1:length(ysi) %loop through steel layers
             Fc_doublecount(1)=0; %we use this to subtract the steel area from the cocnrete, bc we used the entire concrete area to calc Fc
             es_j(j) = phi(count)*(yt-ysi(j))+e(count);   %find strain at layer j
             fs_j = steel_stress_strain(Es,fy,P,esh,esu,es_j(j),fsu);  %find stress at layer j
              if yield_flag_triggered==0  %Flag function to find cracking moment and curvature to plot later 
                 if fs_j >= fy
                     yield_flag=1
                     yield_flag_triggered = 1
                     steel_yield = 1 %Initialize to know later on if it was steel or concrete that yeilded first 
                 end
             end
             if es_j(j) <0     %Account for concrete that was accounted for cocnrete loop  but is occupied by rebar (must be removed)
                 Fc_doublecount(j) = Asi(j)*stress_unconfined_conc(es_j(j), epcc, fpcc, lambda_conc, ecu, rc);
             else 
                 Fc_doublecount(j) = Asi(j)*tens_conc(Ec,fcr,es_j(j),ecr);
             end
             Fs(j) = fs_j*Asi(j) - Fc_doublecount(j);
              
        end
        F_tot = sum(Fc)+sum(Fs); %total force [kip]
        error = (F_tot-P_total)/(sum(Asi)*fy);
        %count_test= count_test+1;  %debugging 
        %count_vec(count_test) = count_test   %debugging 
        %error_plot(count_test) = error   %debugging 
        %plot(count_vec,error_plot, 'bo')
        hold on
        [de,e(count)]=brute_force(error,tol,de,e(count));
    end

    M(count) = sum(Fc.*(yt-yc))+sum(Fs.*(yt-ysi));

    if M(count)>M_max
        M_max = M(count);
        phi_max = phi(count);
    end
    if fcr_flag==1        %If loop picks up on flag concrete cracking triggered and saves the moment and curvature at this step 
        M_cr = M(count);
        phi_cr = phi(count);
        fcr_flag=0;
    end
    if yield_flag==1    %If loop picks up on flag triggered for  yielding saves the moment and curvature at this step 
        M_1styield = M(count);
        phi_1styield = phi(count);
        if steel_yield == 1
            disp("steel yielded first");
        else
            disp("concrete yielded first");
        end
        yield_flag = 0;
    end
    if ACI_flag==1        %If loop picks up on flag ACI limit and saves the moment and curvature at this step 
        M_ACI = M(count);
        phi_ACI = phi(count);
        ACI_flag=0;
    end
    if nom_moment_flag==1        %If loop picks up on flag nominal moment limit and saves the moment and curvature at this step 
        M_nom_moment = M(count)
        phi_nom_moment = phi(count)
        nom_moment_flag=0
    end
    if spall_flag==1        %If loop picks up on flag spalling onset  and saves the moment and curvature at this step 
        M_spall = M(count);
        phi_spall = phi(count);
        spall_flag=0;
    end
    es_max = max(es_j);
    ec_min = min(ec_j);
    count = count + 1 
end

figure(1)
plot(phi,M)
hold on
plot(phi_cr, M_cr,'ro')
hold on
plot(phi_1styield,M_1styield,'bo')
hold on
plot(phi_ACI,M_ACI,'mo')
hold on
plot(phi_nom_moment, M_nom_moment,'go')
hold on 
plot(phi_spall,M_spall,'ko')
hold on
plot(phi_max,M_max,'cyano')
hold on
plot(phi(end),M(end),'color' , [0.8500, 0.3250, 0.0980], 'Marker','s')
hold on 
xlabel('Curvature [1/in]')
ylabel("Moment [kip*in]")


%%%%%%Part2.3%%%%%%%
%Find yield point 
x_vec1 = [0:0.000001:0.0002];
y_vec1 = x_vec1*M_1styield/phi_1styield;
value_found = 0

for i=1:length(y_vec1)
    if y_vec1(i)>M_nom_moment && value_found == 0
        phi_yield = x_vec1(i);
        M_yield = M_nom_moment;
        value_found = 1;
    end
end

%%%%%%Part2.4%%%%%%%
slope_cr = M_cr/phi_cr
flexural_rigidity_transformed = Ec*It
error_It = abs(flexural_rigidity_transformed-slope_cr)/slope_cr

%%%%%%Part2.5%%%%%%%
%add bi_linear and quad-linear piecewise 
x_bi = [0; phi_yield; phi_max; phi(end)];
y_bi = [0; M_yield; M_max; M(end)];
plot(x_bi, y_bi);
hold on
x_quad_lin = [0; phi_cr; phi_yield; phi_max; phi(end)];
y_quad_lin = [0; M_cr; M_yield; M_max; M(end)];
plot(x_quad_lin,y_quad_lin);

legend("M-Phi curve", "Concrete Cracking", "First Yield","ACI Strain Limit","Nominal Moment", "Onset of Spalling", "Overstrength Point", "Ultimate Point",'Bi-Linear Model',"Four-Line Piecewise Model", "location","SouthEast")
grid on

%%%%%%Part2.6%%%%%%%
effective_flexural_stiffness = M_1styield/phi_1styield;
stiffness_modifier = effective_flexural_stiffness/slope_cr;

%%%%%%Part2.7%%%%%%%
%scale M and phi so they are unitless
lambda = 2.25;
%M_Priestley= M/(D^3*fpc)
phi_Priestley= lambda*abs(ey)/D
phi_yield
%plot(phi_Priestley, M_Priestley)
error = abs(phi_Priestley-phi_yield)/phi_yield


            

