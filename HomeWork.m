% ---------------------- FORMULAS -------------------------

% X = Thrust - Drag * cos(alpha) + Lift *sin(alpha)
% For small alpha : 
% X = Thrust - Drag + Lift*alpha = 1/2 * rho*S*V^2* (Ct-Cd+Cl*alpha)

%Z = -Drag(D) * sin(alpha) - Lift(L) * cos(alpha)
% For small alpha:
% Z = -(L+D*alpha) = -1/2 * rho*S*V^2* (Cl+Cd*alpha)  

% Roll Moment (L) = 1/2 * rho*S*V^2*b*Cl 
%Pitch Moment (M) = 1/2 * rho*S*V^2*c_bar*Cm
% Yaw Moment (N) = 1/2 * rho*S*V^2*b*Cn

% Xu = -[(q_bar*S)/(m*U0)]*(2*(Cd)0 + CDu) with q_bar = 1/2*rho*U0^2
% Xw = [(q_bar*S)/(m*U0)]*(CL0*(1-(2/pi*e*A)*CLalpha) with q_bar = 1/2*rho*U0^2
% Zu = -[(q_bar*S)/(2*CL0 + CLu) with q_bar = 1/2*rho*U0^2
% Zw = -[(q_bar*S)/((CD)0 + CLalpha) with q_bar = 1/2*rho*U0^2

% Mu = [(q_bar*S*C_bar)/(Iyy* U0)]*(CMu)
% Mw = [(q_bar*S*C_bar)/(Iyy* U0)]*(CMalpha)
% -Mw_dot = [(q_bar*S*C_bar^2)/(2*Iyy* U0^2)]*(CM_alpha_dot)
% -Mq = [(q_bar*S*C_bar)/(Iyy* U0)]*(CMu)


% GEOMETRIC DATA
S = 34.84; % WingArea in m^2
b = 11.8;% Wingspan in m
A = 4; % Aspect Ratio
c_bar = 3.29; % Mean aerodynamic chord in m
e = 0.85;

% MASS &  INERTIAL DATA
m = 9926.7; % Aircraft mass in kg
Ixx = 18486.6; % kgm^2
Iyy = 68965; % kgm^2
Izz = 91599; %kgm^2
Ixz = 3976.6; %kgm^2
g = 9.81;

% FLIGHT CONDITIONS
H = 0; %Altitude in m
rho = 1.225; % Air Density kg/m^3
U0 = 85; % Flight Speed in m/s
Mach = 0.25;
TETA0 = 11.2*(pi/180); % Initial Pitch Angle

% INITIAL STEADY STATE COEFFICIENTS
Cl0 = 0.62; % Lift Coefficient
Cd0 = 0.072; % Drag Coefficient

% LONGITUDINAL STABILITY DERIVATIVES
CDu = 0.0105;
CDalpha = 1.52;
CLu = 0.04;
CLalpha = 3.95;
CLalphadot = 0;
CLq = 0;
CMu = 0.012;
CMalpha = -0.45;
CMalphadot = -0.7;
CMq = -3.8;

% LONGITUDINAL CONTROL DERIVATIVES
CDdelta_e = 0;
CLdelta_e = 0.6;
CMdelta_e = -0.83;

% LATERAL STABILITY DERIVATIVES
Cy_beta = -0.88;
Cy_p = 0;
Cy_r = 0;
Cl_beta = -0.115;
Cl_p = -0.25;
Cl_r = 0.18;
Cn_beta = 0.105;
Cn_p = -0.01;
Cn_r = -0.34;

% LATERAL CONTROL DERIVATIVES
CYdelta_a = -0.015;
CLdelta_a = 0.055;
CNdelta_a = -0.002;
CYdelta_r = 0.23;
CLdelta_r = 0.007;
CNdelta_r = -0.105;

% ================================================================================================================================================
% I - Longitudinal Dynamic Stability Of Airplane
% ================================================================================================================================================

%-------------------------------- 1 ---------------------------------------

% A = [Xu           Xw           0            -gcos(teta0)
%      Zu           Zw           U0           -gsin(teta0)
%      Mu+Zu*Mw_dot Mw+Zw*Mw_dot Mq+U0*mw_dot 0
%      0            0            1            0           ]

q_bar = 0.5*rho*U0^2;
Xu = -((q_bar*S)/(m*U0))*(2*Cd0 + CDu);
Xw = ((q_bar*S)/(m*U0))*(Cl0*(1-(2*CLalpha/(pi*e*A)))); 
Zu = -((q_bar*S)/(m*U0))*(2*Cl0 + CLu); 
Zw = -((q_bar*S)/(m*U0))*(Cd0 + CLalpha);
Mu = ((q_bar*S*c_bar*CMu)/(Iyy* U0));
Mw = ((q_bar*S*c_bar*CMalpha)/(Iyy* U0));
Mw_dot = ((q_bar*S*c_bar^2*CMalphadot)/(2*Iyy* U0^2));
Mq = ((q_bar*S*c_bar^2*CMq)/(Iyy*U0*2));

Along = [Xu           Xw           0            -g*cos(TETA0);
     Zu           Zw           U0           -g*sin(TETA0);
     Mu+Zu*Mw_dot Mw+Zw*Mw_dot Mq+U0*Mw_dot 0;
     0            0            1            0           ];

%B = [Xdelta_e                 Xdelta_t
%     Zdelta_e                 Zdelta_t
%     Mdelta_e+Zdelta_e*Mw_dot Mdelta_t+Zdelta_t*Mw_dot
%     0                        0];

Xdelta_t = 0;
Zdelta_t = 0;
Mdelta_t = 0;
Xdelta_e = ((q_bar*S)/(m*U0))*(CDdelta_e);
Zdelta_e = ((q_bar*S)/(m*U0))*(CLdelta_e);
Mdelta_e = ((q_bar*S*c_bar)/(Iyy*U0))*(CMdelta_e);

B = [Xdelta_e                 Xdelta_t
     Zdelta_e                 Zdelta_t
     Mdelta_e+Zdelta_e*Mw_dot Mdelta_t+Zdelta_t*Mw_dot
     0                        0];

%-------------------------------- 2 ---------------------------------------

eq_char = charpoly(Along);

%-------------------------------- 3 ---------------------------------------

[eigenvectors, eigenvalues] = eig(Along);

%-------------------------------- 4 ---------------------------------------
% A

%I.	Short period mode (Natural Frequency, Damping Factor)
fprintf("\nSHORT PERIOD MODE");
WNsp = sqrt(eigenvalues(1,1)*eigenvalues(2,2));
Esp = (-eigenvalues(1,1)-eigenvalues(2,2))/(2*WNsp);

fprintf('\nWNsp = %f', WNsp);
fprintf('\nEsp = %f', Esp);

if Esp >0
    fprintf('\nDamping factor positive ==> MODE STABLE\n')
else
    fprintf('\nDamping factor positive ==> MODE UNSTABLE\n')
end

%Period
Tsp = (2*pi)/(WNsp*sqrt(1-Esp^2));

% b.	Phugoid mode (Natural Frequency, Damping Factor) 
fprintf("\nPHUGOID MODE");

WNp = sqrt(eigenvalues(3,3)*eigenvalues(4,4));
Ep = (-eigenvalues(3,3)-eigenvalues(4,4))/(2*WNp);

fprintf('\nWNp = %f', WNp);
fprintf('\nEp = %f', Ep);

if Ep >0
    fprintf('\nDamping factor positive ==> MODE STABLE\n')
else
    fprintf('\nDamping factor positive ==> MODE UNSTABLE\n')
end

%Period
Tp = (2*pi)/(WNp*sqrt(1-Ep^2));

%-------------------------------- 5 ---------------------------------------
% B
t = 0:1:2000;
t_sp = 0:0.01:10;

% Axial velocity
DeltaU = real(eigenvectors(1, 3))*exp(eigenvalues(3,3)*t) + real(eigenvectors(1, 4))*exp(eigenvalues(4,4)*t);

% Pitch rate
DeltaQ = real(eigenvectors(3, 1))*exp(eigenvalues(1, 1)*t_sp) +  real(eigenvectors(3, 2))*exp(eigenvalues(2, 2)*t_sp);

% Angle of attack
DeltaW = real(eigenvectors(2, 1))*exp(eigenvalues(1, 1)*t_sp) +  real(eigenvectors(2, 2))*exp(eigenvalues(2, 2)*t_sp);

% Pitch angle
DeltaTETA = real(eigenvectors(4, 3))*exp(eigenvalues(3, 3)*t) + real(eigenvectors(4, 4))*exp(eigenvalues(4, 4)*t);


% SHORT PERIOD MODE

plot(t_sp, DeltaW);
xlabel('Time (s)');
ylabel('DeltaW');
title('Short period mode - Angle of attack');
legend('DeltaW');
grid();

figure(1)
plot(t_sp, DeltaQ);
hold on;
plot(t_sp, DeltaW);
xlabel('Time (s)');
ylabel('DeltaQ, DeltaW');
title('Short period mode - Pitch rate & Angle of attack');
legend('DeltaQ', 'DeltaW');
grid();


% PHUGOID MODE

figure(2)
plot(t, DeltaTETA);
hold on;
plot(t, DeltaU);
xlabel('Time (s)');
ylabel('DeltaTETA, DeltaU');
title('Phugoid mode - Pitch angle & axial velocity');
legend('DeltaTETA', 'DeltaU');
grid();


%-------------------------------- 6 ---------------------------------------
syms delta;

s = tf('s');
TF = minreal(inv(s*eye(4)-Along)*B, 1e-4);

DeltaU_Dde = zpk(TF(1,1));
DeltaU_Ddt = zpk(TF(1,2));
DeltaW_Dde = zpk(TF(2,1));
DeltaW_Ddt = zpk(TF(2,2));
DeltaQ_Dde = zpk(TF(3,1));
DeltaQ_Ddt = zpk(TF(3,2));
DeltaTETA_Dde = zpk(TF(4,1));
DeltaTETA_Ddt = zpk(TF(4,2));


figure(3)
bode(DeltaU_Dde)
grid on;
title("bode DeltaU_Dde");

figure(4)
bode(DeltaW_Dde)
grid on;
title("bode DeltaW_Dde");

figure(5)
bode(DeltaQ_Dde)
grid on;
title("bode DeltaQ_Dde");

figure(6)
bode(DeltaTETA_Dde)
grid on;
title("bode DeltaTETA_Dde");

%--------------------------------------------------------------------------
syms k1 k2

% Given matrices and values
M1 = [-0.7349 85.00; -0.0381 -0.6405];
M2 = [0.1096; -0.0719];
K = [k1 k2]; % Assuming k1 and k2 are symbolic variables

% Compute the matrix operation
result = M1 - M2 * K;

% Définir les équations du système
eq1 = 0.1096*k1 - 0.0719*k2 + 1.3755 == 3.36;
eq2 = 3.7092675 - 6.0413012*k1 -0.05702226*k2 == 20;

% Résoudre le système
solution = solve([eq1, eq2], [k1, k2]);
k1_decimal = double(solution.k1);
k2_decimal = double(solution.k2);

disp(k1_decimal);
disp(k2_decimal);

K =[0 k1_decimal k2_decimal 0; 0 0 0 0];
A_de = Along -B*K;


[eigenvectors, eigenvalues] = eig(A_de);

[eigenvectors1, eigenvalues1] = eig(Along);

%-------------------------------- 4 ---------------------------------------
% A

%Period
Tp = (2*pi)/(WNp*sqrt(1-Ep^2));

%-------------------------------- 5 ---------------------------------------
% B
t = 0:1:2000;
t_sp = 0:0.01:10;

% Pitch rate rectified
DeltaQ_rectified = real(eigenvectors(3, 1))*exp(eigenvalues(1, 1)*t_sp) +  real(eigenvectors(3, 2))*exp(eigenvalues(2, 2)*t_sp);
% Pitch rate
DeltaQ = real(eigenvectors1(3, 1))*exp(eigenvalues1(1, 1)*t_sp) +  real(eigenvectors1(3, 2))*exp(eigenvalues1(2, 2)*t_sp);

% Angle of attack
DeltaW_rectified = real(eigenvectors(2, 1))*exp(eigenvalues(1, 1)*t_sp) +  real(eigenvectors(2, 2))*exp(eigenvalues(2, 2)*t_sp);
% Angle of attack
DeltaW = real(eigenvectors1(2, 1))*exp(eigenvalues1(1, 1)*t_sp) +  real(eigenvectors1(2, 2))*exp(eigenvalues1(2, 2)*t_sp);



% SHORT PERIOD MODE
figure(1)
plot(t_sp, DeltaW_rectified, 'g');
hold on;
plot(t_sp ,DeltaW, 'r')
xlabel('Time (s)');
ylabel('DeltaW');
title('Short period mode - Angle of attack');
legend('DeltaW rectified', 'DeltaW');
grid();

figure(2)
plot(t_sp, DeltaQ_rectified, 'g');
hold on;
plot(t_sp, DeltaQ, 'r');
xlabel('Time (s)');
ylabel('DeltaQ');
title('Short period mode - Pitch rate');
legend('DeltaQ rectified', 'DeltaQ');
grid();


% PHUGOID MODE
%{
plot(t, DeltaU);
xlabel('Time (s)');
ylabel('DeltaU');
title('Phugoid mode - Axial velocity');
legend('DeltaU');
grid();

figure(2)
plot(t, DeltaTETA);
xlabel('Time (s)');
ylabel('DeltaTETA');
title('Phugoid mode - Pitch angle');
legend('DeltaTETA');
grid();
%}




% ================================================================================================================================================
% II - Lateral Dynamic Stability Of Airplane
% ================================================================================================================================================

%-------------------------------- 1 ---------------------------------------
Yv = (q_bar*S)/(m*U0)*Cy_beta;
Yp = (q_bar*S*b)/(2*m*U0)*Cy_p;
Yr = (q_bar*S*b)/(2*m*U0)*Cy_r;
Lv = (q_bar*S*b)/(Ixx*U0)*Cl_beta;
Lp = (q_bar*S*b^2)/(2*Ixx*U0)*Cl_p;
Lr = (q_bar*S*b^2)/(2*Ixx*U0)*Cl_r;
Nv = (q_bar*S*b)/(Izz*U0)*Cn_beta;
Np = (q_bar*S*b^2)/(2*Izz*U0)*Cn_p;
Nr = (q_bar*S*b^2)/(2*Izz*U0)*Cn_r;

Alateral = [Yv           Yp           -(U0-Yr)      g*cos(TETA0);
             Lv           Lp           Lr           0;
             Nv           Np           Nr           0;
             0            1            0            0           ];

Ydelta_r = (q_bar*S)/(m*U0)*CYdelta_r;
Ydelta_a = (q_bar*S)/(m*U0)*CYdelta_a;
Ldelta_r = (q_bar*S*b)/(Ixx*U0)*CLdelta_r;
Ldelta_a = (q_bar*S*b)/(Ixx*U0)*CLdelta_a;
Ndelta_r = (q_bar*S*b)/(Izz*U0)*CNdelta_r;
Ndelta_a = (q_bar*S*b)/(Izz*U0)*CNdelta_a;

Blateral = [Ydelta_r                 Ydelta_a
     Ldelta_r                 Ldelta_a
     Ndelta_r                 Ndelta_a
     0                        0];

%-------------------------------- 2 ---------------------------------------

eq_char_lat = charpoly(Alateral);

%-------------------------------- 3 ---------------------------------------

[eigenvectors2, eigenvalues2] = eig(Alateral);

%-------------------------------- 4 ---------------------------------------

% Identification des valeurs propres pour chaque mode
rolling_eigenvalue = eigenvalues2(1, 1);
spiral_eigenvalue = eigenvalues2(2, 2);
dutch_roll_eigenvalue = eigenvalues2(3, 3);

% Calcul des fréquences naturelles et des facteurs d'amortissement
rolling_natural_frequency = abs(imag(rolling_eigenvalue));
rolling_damping_ratio = -real(rolling_eigenvalue) / abs(rolling_eigenvalue);

spiral_natural_frequency = abs(imag(spiral_eigenvalue));
spiral_damping_ratio = -real(spiral_eigenvalue) / abs(spiral_eigenvalue);

dutch_roll_natural_frequency = abs(imag(dutch_roll_eigenvalue));
dutch_roll_damping_ratio = -real(dutch_roll_eigenvalue) / abs(dutch_roll_eigenvalue);

% Affichage des résultats
disp('Rolling Mode:');
disp(['Natural Frequency: ', num2str(rolling_natural_frequency), ' rad/s']);
disp(['Damping Ratio: ', num2str(rolling_damping_ratio)]);

disp('Spiral Mode:');
disp(['Natural Frequency: ', num2str(spiral_natural_frequency), ' rad/s']);
disp(['Damping Ratio: ', num2str(spiral_damping_ratio)]);

disp('Dutch Roll Mode:');
disp(['Natural Frequency: ', num2str(dutch_roll_natural_frequency), ' rad/s']);
disp(['Damping Ratio: ', num2str(dutch_roll_damping_ratio)]);

%-------------------------------- 5 ---------------------------------------

% Time vector
t = 0:0.01:30;
r_lat = roots(eq_char_lat);

% Dutch roll mode
v_dutch = real(eigenvectors2(1,3))*exp(eigenvalues2(3,3)*t)+real(eigenvectors2(1,2))*exp(eigenvalues2(2,2)*t);
figure (7)
plot(t,v_dutch)
hold on
p_dutch = real(eigenvectors2(2,3))*exp(eigenvalues2(3,3)*t)+real(eigenvectors2(2,2))*exp(eigenvalues2(2,2)*t);
plot(t,p_dutch)
hold on
r_dutch = real(eigenvectors2(3,3))*exp(eigenvalues2(3,3)*t)+real(eigenvectors2(3,2))*exp(eigenvalues2(2,2)*t);
plot(t,r_dutch)
hold on
phi_dutch = real(eigenvectors2(4,3))*exp(eigenvalues2(3,3)*t)+real(eigenvectors2(4,2))*exp(eigenvalues2(2,2)*t);
plot(t, phi_dutch)
xlabel ('t')
title('Dutch roll mode')
grid()
legend('lat velocity', 'roll rate', 'yaw rate', 'roll angle')


%roll
figure(8)
t=0:0.1:4;
hold on;
v = real(eigenvectors2(1,1))*exp(t*r_lat(1));
plot(t,v)
hold on;
p = real(eigenvectors2(2,1))*exp(t*r_lat(1));
plot(t,p)
hold on;
r = real(eigenvectors2(3,1))*exp(t*r_lat(1));
plot(t,r)
hold on;
phi = real(eigenvectors2(4,1))*exp(t*r_lat(1));
plot(t,phi)
xlabel ('t')
title('Rolling mode')
grid()
legend('lat velocity', 'roll rate', 'yaw rate', 'roll angle')



%spiral
figure(9)
t=0:200;
x_sp = exp(r_lat(4)*t);
hold on;
v_sp = real(eigenvectors2(1,4))*exp(t*r_lat(4));
plot(t,v_sp)
hold on;
p_sp = real(eigenvectors2(2,4))*exp(t*r_lat(4));
plot(t,p_sp)
hold on;
r_sp = real(eigenvectors2(3,4))*exp(t*r_lat(4));
plot(t,r_sp)  
hold on;
phi_sp = real(eigenvectors2(4,4))*exp(t*r_lat(4));
plot(t,phi_sp)
xlabel ('t')
title('Spiral mode')
grid()
legend('lat velocity', 'roll rate', 'yaw rate', 'roll angle')



%-------------------------------- 6 ---------------------------------------

syms delta;

s = tf('s');
TF_lat = minreal(inv(s*eye(4)-Alateral)*Blateral, 1e-4);

DeltaV_Ddr = zpk(TF_lat(1,1));
DeltaV_Dda = zpk(TF_lat(1,2));
DeltaP_Ddr = zpk(TF_lat(2,1));
DeltaP_Dda = zpk(TF_lat(2,2));
DeltaR_Ddr = zpk(TF_lat(3,1));
DeltaR_Dda = zpk(TF_lat(3,2));
Deltaphi_Ddr = zpk(TF_lat(4,1));
Deltaphi_Dda = zpk(TF_lat(4,2));


figure(11)
bode(DeltaV_Ddr)
grid on;
title("bode DeltaV_Ddr - Lateral velocity rudder");

figure(12)
bode(DeltaP_Ddr)
grid on;
title("bode Deltap_Ddr_lat - Roll rate rudder");

figure(13)
bode(DeltaR_Ddr)
grid on;
title("bode DeltaR_Ddr - Yaw rate rudder");

figure(14)
bode(Deltaphi_Ddr)
grid on;
title("bode Deltaphi_Ddr - Roll angle rudder");

figure(15)
bode(DeltaV_Dda)
grid on;
title("bode DeltaV_Dda - Lateral velocity ailerons");

figure(16)
bode(DeltaP_Dda)
grid on;
title("bode DeltaP_Dda - Roll rate ailerons");

figure(17)
bode(DeltaR_Dda)
grid on;
title("bode DeltaR_Dda - Roll angle ailerons");

figure(18)
bode(Deltaphi_Dda)
grid on;
title("bode Deltaphi_Dda - Roll angle ailerons");

%}