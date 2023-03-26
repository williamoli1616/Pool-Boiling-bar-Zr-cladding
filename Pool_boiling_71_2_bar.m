%addpath('./funzioni_esame')
close all
%Fuel assembley geometry
pitch = 14.37; %mm

%Cladding
clad_length = 4.09; %m
clad_length =3.588;
clad_diam_outer = 11.20/1000; %m
clad_thickness = 0.71/1000; %m
clad_diam_inner = clad_diam_outer-2*clad_thickness; %m
Q_core =  3500*1000000; %MWth
N_rods = 81; %Number of fuel rods
N_assemblies = 764;
Q_assemblies = Q_core / N_assemblies;
Q_rod = Q_assemblies/N_rods;
rho_Zr = 6511; %Kg/m3; %Density of the cladding 
%V_clad_1 = (pi*((clad_diam_outer/2).^2)-pi*((clad_diam_inner/2).^2))*(clad_length); %Volume of the cladding in m3
V_clad = pi*((clad_diam_outer/2).^2)*clad_length;
cs_Zr = 270; % Specific heat capacity of Zirconium J/KgK
Mass_Zr = V_clad*rho_Zr; %Mass of the cladding Kg
Cs = Mass_Zr*cs_Zr; %Capacity of Zirconium J/K
Wetted_area = clad_diam_outer*pi*clad_length; %Outer wetted area of the cladding
T_wall_0 = 287.1733; %Initial temperature of the wall of the cladding in %°C 
%k_Zr= 8.8527+7.0820*(10^-3)*T+2.5329*(10^-6)*(T^2)+2.9918*(10^3)*(T^(-1); %Formula for thermal conductivity of zirconium as a function of temperature

%Coolant characteristics
P_in = 71.4; %bar
T_in = 287.1732; % °C
mu = 0.0000906; % Viscosity of water in Pascal*s
rho_in = XSteam('rhoL_p',P_in); % Density of water Kg/ m3
g = 9.81; % Accelleration gravity m/s2
k_w = 0.68; % Therman conductivity in W/m°C or W/mK of water
cp= XSteam('CpL_p',P_in); % kJ/(kg°C) Specific heat capacity of my fluid 
sigma= XSteam('st_p',P_in); % Surface tension in N/m
hg = XSteam('hV_p',P_in); %Enthalpy of saturated vapour in kJ/kg
hl = XSteam('hL_p',P_in); %Enthalpy of saturated liquid in kJ/kg
hlg = hg-hl; %kJ/kg
vg = XSteam('vV_p',P_in); %Specif volume of the gas bubble m3/Kg
beta = 7.51E-04; %Thermal expansion coefficients

%Celsius quantities converted into Kelvin and quantities in kJ/kg-->J/kg;
%kJ/kgK---->J/kgK
T_in = T_in + 273.15; %Converting the intial temperuature of the fluid in Kelvin
T_wall_0=T_wall_0 + 273.15; %Converting the Twall initial in Kelvin
cp = cp*1000;
hg = hg*1000;
hl = hl*1000;
hlg = hlg*1000;


%Natural convection of the fluid with the zirconium cladding
C_nc = 0.76; % Multiplication constant for Nu=C*(Ra)^n
Lc = clad_length; %Characteristic length of natural convection
n = 0.25; %Exoponent of the Rayileigh correlation
Q_rod = Q_rod*0.065; %Watt power of the fuel rod;
Pr = (cp*mu/k_w);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Grouping of constant values for the SIMPLE MODEL OF NATURAL CONVECTION
konstant = (Wetted_area/Cs) * (k_w/Lc) * C_nc * Pr^n * beta^n * ((Lc^3*g*rho_in^2)/(mu^2))^n;

%Time necessary to achieve onset of nucleate boling;
%Calculation of the wall temperature in natural convection
f = @(t,T) Q_rod/(Cs) - konstant*((T-T_in).^(n+1));
t0=0;
t_max=5;
h=0.0001;
[time_nc,T_wall_nc]=eulero_avanti(f,t0,t_max,T_wall_0,h);

for ii =1:1:length(T_wall_nc)   
alfa_conv(ii) = (konstant*Cs/Wetted_area)*((T_wall_nc(ii)-T_in)).^n;
qs_nc(ii) = alfa_conv(ii).*(T_wall_nc(ii)-T_in); %Heat flux calculated in natural convection W/m2
k_Zr(ii)= 8.8527 + 7.0820.*(10.^-3).*T_wall_nc(ii)+ 2.5329.*(10^-6).*(T_wall_nc(ii).^2)+2.9918.*(10.^3).*(T_wall_nc(ii).^(-1)); %Thermal conductivity of Zirconium
Bi(ii) = alfa_conv(ii)*(V_clad/Wetted_area)/k_Zr(ii); %Biot Number check up for uniform temperature
Gr(ii) = (beta*rho_in^2 *g *Lc^3.*(T_wall_nc(ii)-T_in)./(mu.^2)).^(n);
Ra(ii) = Gr(ii)*Pr; %Rayleigh check up
end

%Equation that describes "Heat flux" as a function of "Wall Temperature" in the
%natural convection region

int = linspace(min(T_wall_nc),max(T_wall_nc),100);
cnc = polyfit(T_wall_nc, qs_nc, 2); %Coefficients for the Heat flux as a function of temperature 


%Solve for the DELTA T CRITICAL for the onset of nucleate boling
constant1 = ((C_nc*k_w)/Lc)*(((g * beta/mu *Lc^3*rho_in^2))*(cp*mu/k_w))^n;
constant2 = ((hlg*k_w)/(8*sigma*T_in*vg));
fun = @(T) constant1.*(T-T_in).^(n-1) - constant2.*(T-T_in).^2;
nmax = 100;
toll=1e-5;
dfun = @(T) constant1.*(n-1).*(T-T_in).^(n-2) - 2.*constant2.*(T-T_in);
T0 = 570.3232;
[onb_new,it_new] = newton(T0,nmax,toll,fun,dfun);
T_Davis_Henderson = onb_new(end) % °K;
T_Davis_Henderson_Celsius= onb_new(end)-273.15 %°C
rc_max = 5*10^-6;

%Critical Raidus for ONB
qs_onb_Davis_Henderson = abs(polyval(cnc,T_Davis_Henderson)); %Heat exhcanged at the onset of nucleate boiling, at temperature T_Davis_Hendeson
rc = sqrt((2*sigma*T_in*vg*k_w)/(hlg*qs_onb_Davis_Henderson)); %Critical radii for the the qs of Davis Henderson
constant3 = (2*sigma*T_in*vg)/(hlg*rc_max);   

if rc > rc_max    
          %Solve for Wall temperature such that all cavities are activated.
          fun_onb = @(T) T-T_in - (rc_max/k_w).*constant1.*(T-T_in).^(n+1) -constant3
          T_1 = 550;
          T_2 = 700;
          nmax=1000;
          toll=1e-6;
          [Temp_onb,xdif,fx,it]=bisez(T_1,T_2,nmax,toll,fun_onb);
          T_onb = Temp_onb(end)
          T_onb_celsius = T_onb-273.15
          q_ONB=10.1628
      
         else
          T_onb=T_Davis_Henderson
          T_onb_celsius = T_Davis_Henderson_Celsius
          q_ONB = qs_onb_Davis_Henderson
         
      end      

%Calculation of the Heat Flux at the onset of nucleate boiling since the
%Interpolated curve at values close to Twall_initial gives a wrong
%estimation of the heat exchanged


%Plotting the "Wall Temperature as a function of time"
time = linspace(0,3,100);
t_coeff = polyfit(time_nc, T_wall_nc-273.15, 10); %Coefficients for the Twall vs time
temp_finder = @(t) polyval(t_coeff,t)-T_onb_celsius; %Function for the Bisection method

%Check for the time when ONB starts;
%figure(2)
%hold on
%plot(time,polyval(t_coeff,time),'k')
%plot(time,T_onb_celsius.*ones(1,length(time)),'k')
%hold off

t_1 = 0.00001;
t_2 = 100;
nmax=1000;
toll=1e-6;
[Time_for_ONB,xdif,fx,it]=bisez(t_1,t_2,nmax,toll,temp_finder);
Time_for_ONB = Time_for_ONB(end) %Seconds

%%%%%%%%%%%%%% Plotting the evolution of the HEAT FLUX vs TEMPERATURE
figure(1)
hold on
plot(T_wall_nc, qs_nc,LineWidth=2)
xlabel('Temp [°K]')
ylabel('Heat flux [W/m2]')
%plot(int,polyval(cnc,int),'-')
plot(T_wall_nc,q_ONB*ones(1,length(T_wall_nc)))%Making sure that the Interpolating curve for qs_nc, is well rapresented by the Polyval
legend('Heat flux', 'Heat flux at Onset of Nucleate Boiling');
grid on
hold off
%%%%%%%%%%%%%% Plotting the evolution in TIME of the WALL TEMPERATURE
figure(2)
hold on
plot(time_nc, T_wall_nc-273.15,'b',LineWidth=2)
xlabel('Time [s]')
ylabel('T-wall-NC [°C]')
plot(time_nc,T_onb_celsius*ones(1,length(time_nc)))
%plot(time, polyval(t_coeff,time),'-')
xlabel('Time [s]')
ylabel('T-wall-NC [°C]')
legend('Temperature curve in Natural Convection', 'Temperature Onset of Nucleate Boiling');
grid on
hold off
%%%%%%%%%%%%%% Plotting the evolution in TIME of the Difference between
%%%%%%%%%%%%%% Heat in and Heat out
figure(3)
plot(time_nc, Q_rod/Wetted_area-qs_nc,LineWidth=2)
xlabel('Time [s]')
ylabel('Heat flux difference [°C]')
title('Heat flux difference in Natural convection')
grid on




%% Nucleate boiling region

%%%%%%%%%%%%%%%% 
% Rosenhow Model
%%%%%%%%%%%%%%%%
rho_in_L = XSteam('rhoL_p',P_in); %Saturated liquid density
rho_in_V =XSteam('rhoV_p', P_in); %Saturated vapour density
Csf = 0.013;
deltaRHO=rho_in_L-rho_in_V;
%Grouping the Rosenhow constants qs_nb = R_constant * DeltaT^3
R_constant = mu*hlg*sqrt((g*deltaRHO)/sigma)*(cp/(Csf*hlg*Pr))^3; %Rosenhow constant

%Evaluating the Critical Heat flux using Kutateladze model
Kut_coeff=0.118;
qchf_kut = Kut_coeff*sqrt(rho_in_V)*hlg*(sigma*g*deltaRHO)^(1/4);

%Determining the ONB temperature using the Rosenhow model
T_onb_Rosenhow = (q_ONB/R_constant)^(1/3) + T_in; %Temperature for ONB using Rosenhow model
T_onb_Rosenhow_celsius = T_onb_Rosenhow -273.15

%Solving the transient differential equation for the nucleate boiling
%region
m=3;
f_nb_ROSENHOW = @(t,T)  Q_rod/(Cs) - ((R_constant*Wetted_area)/Cs)*((T-T_in).^(m));
t0_nb_R=0;
t_max_nb_R=10;
h=0.001;
T_wall_0_nb = T_onb; %Temperature of the ONB in Kelvin
[time_nb_R, T_wall_nb_R]=eulero_avanti(f_nb_ROSENHOW,t0_nb_R,t_max_nb_R,T_wall_0_nb,h);
T_wall_nb_R_Celsius = T_wall_nb_R-273.15;

for ii =1:1:length(T_wall_nb_R)   
qs_nb_R(ii) = R_constant.*((T_wall_nb_R(ii)-T_in).^m); %Heat flux calculated in Nucleate boiling W/m2
alfa_nb_R(ii) = mu*hlg*sqrt(g*(rho_in_L-rho_in_V)/sigma).*((cp.*(T_wall_nb_R(ii)-T_in))./(Csf*hlg*Pr)).^2;
heat_difference_R(ii) = Q_rod - qs_nb_R(ii).*Wetted_area;
end

%%%%%%%%%%%%%%%%%
% Mostinski Model
%%%%%%%%%%%%%%%%%
pc=221.1; %Critical pressure of water in bar
pr = P_in/pc %Reduced pressure
A = 0.1011*pc^0.69;
Fpr = 1.8*pr^(0.17) + 4*pr^(1.2) + 10*pr^(10);
T_onb_MOST = (q_ONB^0.3)/(A*Fpr) + T_in; %Temperature of the onset of nucleate boiling using Mostinski model
T_onb_MOST_celsius = (q_ONB^0.3)/(A*Fpr) + T_in - 273.15 %Temperature of the onset of nucleate boiling using Mostinski model in Celsius
f_nb_MOSTINSKI = @(t,T)  Q_rod/(Cs) - (Wetted_area/Cs)*(A*Fpr*(T-T_in))^(1/0.3);
t0_nb_R=0;
t_max_nb_R=5;
h=0.001;
T_wall_0_nb = T_onb; %Temperature of the ONB in Kelvin
[time_nb_MOST, T_wall_nb_MOST]=eulero_avanti(f_nb_MOSTINSKI,t0_nb_R,t_max_nb_R,T_wall_0_nb,h);
T_wall_nb_MOST_Celsius = T_wall_nb_MOST-273.15;

most_const=(A*Fpr)^(10/3);
for ii =1:1:length(T_wall_nb_MOST)   
qs_nb_MOST(ii) = most_const.*(T_wall_nb_MOST(ii)-T_in).^(10/3); %Heat flux calculated in Nucleate boiling W/m2
heat_difference_MOST(ii) = Q_rod - qs_nb_MOST(ii).*Wetted_area;
end

figure(4)
plot(time_nb_MOST, heat_difference_MOST,'LineWidth',2)
title('Time evaluation for steady state condition Heat Difference Mostinski Model')
xlabel('Time [s]')
ylabel('Heat flux variation [W/m2]')
grid on

figure(5)
plot(time_nb_MOST, T_wall_nb_MOST_Celsius,'LineWidth',2)
%title('Time evaluation for steady state condition Temp Mostinski Model')
xlabel('Time [s]')
ylabel('T-wall-NB-Mostinski [°C]')
grid on

figure(6)
plot(time_nb_R, heat_difference_R,'LineWidth',2)
%title('Time evaluation for steady state condition Heat Difference Rosenhow Model')
xlabel('Time [s]')
ylabel('Heat flux variation [W/m2]')
grid on

figure(7)
plot(time_nb_R, T_wall_nb_R_Celsius,'LineWidth',2)
%title('Time evaluation for steady state condition Temp Rosenhow Model')
xlabel('Time [s]')
ylabel('T-wall-NB-Rosenhow [°C]')
grid on



%% The Limitation of using a method where the temperature profile is uniform over time.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation that predicts when the Rosenhow model shoud stop working at time 't'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stop=1;
ii=0;
while stop>0
ii=ii+1;
qs_nb_R(ii) = R_constant.*(T_wall_nb_R(ii)-T_in).^3; %Heat flux calculated in Nucleate boiling W/m2
alfa_nb_R(ii) = mu*hlg*sqrt(g*(rho_in_L-rho_in_V)/sigma).*((cp*(T_wall_nb_R(ii)-T_in))./(Csf*hlg*Pr)).^(2);
heat_Bi_R(ii) = Q_rod - qs_nb_R(ii).*Wetted_area;
k_Zr(ii)= 8.8527 + 7.0820.*(10.^-3).*T_wall_nb_R(ii)+ 2.5329.*(10^-6).*(T_wall_nb_R(ii).^2)+2.9918.*(10.^3).*(T_wall_nb_R(ii).^(-1));
Bi_R(ii) = alfa_nb_R(ii)*(V_clad/Wetted_area)/k_Zr(ii);
    if Bi_R(ii)>0.1
       stop=-1;
       last_ii = ii;
   
    end
end

t_stop_71 = time_nb_R(ii)
for kk=1:1:ii
    time_Bi_R(kk) = time_nb_R(kk);
    T_wall_Bi_R(kk) = T_wall_nb_R(ii);
end
figure(8)
plot(time_Bi_R, heat_Bi_R)
title('Time evaluation for when Bi < 0.1,Rosenhow model')
grid on
figure(9)
plot(time_Bi_R, T_wall_Bi_R, LineWidth=4)
title('Time evaluation for when Bi < 0.1 for Wall Temp, Rosenhow model')
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation that predicts when the Mostinski model shoud stop working at time 't'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stop=1;
ii=0;
while stop>0
ii=ii+1;
qs_nb_MOST(ii) = most_const.*(T_wall_nb_MOST(ii)-T_in).^(10/3); %Heat flux calculated in Nucleate boiling W/m2
alfa_nb_MOST(ii) = most_const.*(T_wall_nb_MOST(ii)-T_in).^(7/3);
heat_Bi_MOST(ii) = Q_rod - qs_nb_MOST(ii).*Wetted_area;
k_Zr(ii)= 8.8527 + 7.0820.*(10.^-3).*T_wall_nb_R(ii)+ 2.5329.*(10^-6).*(T_wall_nb_R(ii).^2)+2.9918.*(10.^3).*(T_wall_nb_R(ii).^(-1));
Bi_MOST(ii) = alfa_nb_MOST(ii)*(V_clad/Wetted_area)/k_Zr(ii);
    if Bi_MOST(ii)>0.1
       stop=-1;
       last_ii = ii;
   
    end
end

t_stop_71_MOST = time_nb_R(ii)
for kk=1:1:ii
    time_Bi_MOST(kk) = time_nb_R(kk);
    T_wall_Bi_MOST(kk) = T_wall_nb_R(ii);
end
figure(10)
plot(time_Bi_MOST, heat_Bi_MOST)
title('Time evaluation for when Bi < 0.1, Mostinski model')
grid on
figure(11)
plot(time_Bi_MOST, T_wall_Bi_MOST, LineWidth=4)
title('Time evaluation for when Bi < 0.1 for Wall Temp, Mostinski model')
grid on



%% Transient calculation from the Leidenfrost point onwards
rho_in_V = XSteam('rhoV_p', P_in);
P2 = 0.09; %Berenson coefficient to calculate the Leidenfrost point
cpg = XSteam('CpV_p', P_in); %Specific heat capacity of the vapour phase
cpg = cpg*1000;
C2 = 0.68;
C1 = 0.4;
kv = 0.035; %Thermal conductivity of water
kv = kv^3;
mug = 9.0600E-06;
Lc=clad_diam_outer/2; %According to the Thermal Hydraulics of a BWR
film_const = (kv*g*rho_in_V*(rho_in_L-rho_in_V)/(Lc*mug));
%Calculation of the q_min Leidenfrost point
q_leiden = P2*rho_in_V*hlg*(g*sigma*(rho_in_L-rho_in_V)/(rho_in_L+rho_in_V)^2)^(1/4)
% Twall-Tsat at the LeidenFrost point

kv_1 = 0.2509; %Corrected the thermal conductivity
Berenson_Temp=(1/kv_1)*(0.127*rho_in_V*hlg).*((g*(rho_in_L-rho_in_V)./(rho_in_L-rho_in_V)).^(2/3))*((sigma./(g*(rho_in_L-rho_in_V))).^(1/2)).*(mug./(g*(rho_in_L-rho_in_V)))^(1/3);
T_wall_0_film = T_in+Berenson_Temp;
T_wall_0_film_Celsius = T_wall_0_film-273.15
f_film = @(t,T) Q_rod/Cs- C1*(Wetted_area/Cs)*(T-T_in)*((hlg*film_const+C2*cpg*film_const*(T-T_in))/(T-T_in))^(1/4);
t0_film=0;
t_max_film = 1500;
h=0.01;   
[time_film,T_wall_film]=eulero_avanti(f_film,t0_film,t_max_film,T_wall_0_film,h);

figure(13)
plot(time_film, T_wall_film-273.15,'LineWidth',2)
title('Time vs Wall superheating in Film Boiling')
xlabel('Time [s]')
ylabel('T-Wall-Film [°C]')


