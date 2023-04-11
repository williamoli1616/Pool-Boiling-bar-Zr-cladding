%The MATLAB CODE solves the heat equation in cylindrical coordinates
%for a zirconium rod with fix heat generation and heat removal at the
%border by convection. The heat equation has only a dependency on "r"



%Data
k=18.9738; %Zirconium cladding Thermal conductivity; [W/(m K)]
h = 4.0172e+03; %Heat transfer coefficient for convection [W/m2K]
R = 5.6e-3; %Rod length [m]
cs = 270; %Specific heat of Zirconium [[J/Kg K]
rho = 6511; %Density of the Zirconium cladding [Kg/m3]
q =7.0052e+6; %Heat density at the center [W/m3]
T_in = 560.3232; %Saturated liquid Temperature [K]
T0=560.5582; %Initial temperature of the rod [K]

Nr= 300; %Discretization for the radius
Nt = 300; %Discretization for the time
dr = R./Nr;
L = 1000; %seconds
r = linspace(0,R,Nr);
t = linspace(0,L,Nt);


%Simplify calculation
yo = 2*dr*h./k; % 
ma=rho*cs/k; %
na = q/k; %

%Initial condition
Ti = ones(1,Nr).*T0; %Initial temperature when everything stops.

%Solve ordinary differential equation
 
[t, T] = ode15s(@fpde,t,Ti,[],T_in,Nr,ma,na,yo,dr,r)

%recal boundary conditons
T(:,1) = (4*T(2) - T(3))/(3);
T(:,end) = (yo.*T_in + 4*T(end-1) -T(end-2))./(3+yo)

%Plotting
figure(1)
imagesc(r,t,T)
colormap jet
colorbar
title('Temperature profile in degrees K')
xlabel('radius [m]')
ylabel('time [s]')
grid on

figure(2)
contourf(r,t,T,20)
colormap jet
title('temperature profile in degrees K')
xlabel('Radius')
ylabel('time')
grid on


%Function
function dTdt = fpde(t,T,T_in,Nr,ma,na,yo,dr,r)
dTdt = zeros(Nr,1);
T(1) = (4*T(2)-T(3))/3; %Bc 1 %Heat flux at the center is null
T(end) = (yo.*T_in + 4*T(end-1) -T(end-2))./(3+yo); %Bc2
%na./ma = Heat input/rho*cs*k
for i=2:Nr-1
    d2Tdr2(i) = (T(i+1)-2*T(i)+T(i-1))./dr.^2; %Rewriting the second partial derivative
    dTdr(i) = (T(i+1)-T(i-1))./(2.*dr); %Rewriting the the first partial derivative
    dTdt(i) = (1/ma).*d2Tdr2(i) + (1/ma).*(1./r(i)).*dTdr(i) + na./(ma); %Cylindrical Heat equation
end
end





