%Parameters
A_O = 6;
alpha = 3.75*10^(-5);
C_min = 0.01;
C_struct = 0.20;
gamma = 0.5;
epsilon = 0.22;
I_sat = 200;
J_max = 1.4*10^(-4);
k_A = 0.6;
k_dw = 0.0785;
k_C = 2.1213;
k_N = 2.72;
m_1 = 0.1085;
m_2 = 0.03;
my_max = 0.18;
N_min = 0.01;
N_max = 0.022;
N_struct = 0.01;
P_1 = 1.22*10^(-3);
P_2 = 1.44*10^(-3);
a_1 = 0.85;
a_2 = 0.3;
R_1 = 2.785*10^(-4);
R_2 = 5.429*10^(-4);
T_R1 = 285;
T_R2 = 290;
T_AP = 1694.4;
T_APH = 25924;
T_APL = 27774;
T_AR = 11033;
U_065 = 0.03;
K_X = 4;

% state = [Area;
%          Nitrogen;
%          Carbon];

%% Forward Euler

NumberIterations = numColsTime;



% Statevariables
X = [];
for n = 1:NumberIterations
    T = Envdata.T(:,:,1,n) + 272.15;
    U = 1;
    X_NO3 = Envdata.NO3(:,:,1,n);
    I = Envdata.PAR(:,:,1,n);
    X(:, end+1) = [T U X_NO3 I];
end

% Standard variation and mean
[T_std, T_mu] = std(X(1, :));
[U_std, U_mu] = std(X(2, :));
[NO3_std, NO3_mu] = std(X(3, :));
[I_std, I_mu] = std(X(4, :));

Nsample = 10;
Y_A = zeros(Nsample, NumberIterations);
Y_N = zeros(Nsample, NumberIterations);
Y_C = zeros(Nsample, NumberIterations);
% for n=1:Nsample
%     Y = ForwardEuler (N_struct, C_struct, k_A, k_N, k_C, k_dw, N_min, C_min, m_1, m_2, A_O, epsilon, K_X, N_max, J_max, U_065, R_1, T_AR, T_R1, gamma, alpha, I_sat, P_1, T_AP, T_APL,T_APH, NumberIterations, X)
%     Y_A(n,:) = Y(1,:);
%     Y_N(n,:) = Y(2,:);
%     Y_C(n,:) = Y(3,:);
% end

h = 1;
Y = zeros(3, NumberIterations);
%Y_A = zeros(Nsample, NumberIterations);
% Y_N = zeros(Nsample, NumberIterations);
% Y_C = zeros(Nsample, NumberIterations);
% netCarbonFixed = zeros(1, NumberIterations-1);
% grossFrond = zeros(1, NumberIterations-1);
A_0 = 0.01;
N_0 = 0.008;
C_0 = 0.01;
Y(:,1) = [A_0;N_0;C_0];
% Y_A(1) = A_0;
% Y_N(1,1) = N_0;
% Y_C(1,1) = C_0;

for j = 1:Nsample
for n = 1:(NumberIterations-1)
    A = Y(1,n);
    N = Y(2,n);
    C = Y(3,n);
    %     T = X(1, n);
    %     U = X(2,n);
    %     NO3 = X(3,n);
    %     I = X(4,n);
    %     T = (X(1, n) + normrnd(T_mu, T_std))/2;
    %     U = (X(2,n) + normrnd(U_mu, U_std))/2;
    %     NO3 = (X(3,n) + normrnd(NO3_mu, NO3_std))/2;
    %     I = (X(4,n) + normrnd(I_mu, I_std))/2;
    T = X(1, n) + randn;
    U = X(2,n) + randn;
    NO3 = X(3,n) + randn/2;
    I = X(4,n) + randn;
    state_dot = kelp(N_struct, C_struct, k_A, k_N, k_C, k_dw, N_min, C_min, m_1, m_2, A_O, epsilon, K_X, N_max, J_max, U_065, R_1, T_AR, T_R1, gamma, alpha, I_sat, P_1, T_AP, T_APL,T_APH, U, NO3, T, I, N, C, A);
    %     netCarbonFixed(n+1) = netCarbon(W_s, state_dot(3), n, n+1) + netCarbonFixed(n);
    %     grossFrond(n+1) = totalGross(my, A, n, n+1) + grossFrond(n);
    Y(:, n+1) = Y(:,n) + state_dot*h;
    %     Y_A(n+1) = Y_A(n) + state_dot(1)*h;
    %     Y_N(n+1) = Y_N(n) + state_dot(2)*h;
    %     Y_C(n+1) = Y_C(n) + state_dot(3)*h;
end
Y_A(j,:) = Y(1,:);
Y_N(j,:) = Y(2,:);
Y_C(j,:) = Y(3,:);
end

plot(time, Y_A);
plot(time, Y_N);
plot(time, Y_C);

%% Plots
% figure(1);
% plot(time, Y(1,:));
% title('Area');
% 
% figure(2);
% plot(time, Y(2,:));
% title('Nitrogen');
% 
% figure(3);
% plot(time, Y(3,:));
% title('Carbon');

% figure (4)
% plot(time, netCarbonFixed);
%
% figure (5)
% plot(time, grossFrond);




%% Net Carbon Fixed
function C = netCarbon (W_s, dCdt, tmin, tmax)
syms t
fun = W_s * dCdt;
C = int(fun, t, tmin, tmax);
end

%% Total gross frond area
function A_tot = totalGross (my, A, tmin, tmax)
syms t
fun = my*A;
A_tot = int(fun, t, tmin, tmax);
end

%% Effect of temperature
function f_temp = EffectTemp (T_k)
T = T_k - 272.15;
if (T>=-1.8) && (T < 10)
    f_temp = 0.08*T + 0.2;
elseif (T>= 10) && (T<= 15)
    f_temp = 1;
elseif (T>15) && (T<=19)
    f_temp = 19/4 - T/4;
elseif (T>19)
    f_temp = 0;
end
end

%% Gross photosynthesis
function P = grossPhotosynthesis ( alpha, I_sat, P_1, T_AP, T_R1, T_APL,T_APH, T, I)
syms B;
T_PL = 271;
T_PH = 296;
P_max = (P_1*exp(T_AP/T_R1-T_AP/T)) / (1+exp(T_APL/T- T_APL/T_PL) + exp(T_APH/T_PH - T_PH/T));

fun = @(B) - (alpha*I_sat) / (log(1+alpha/B)) * (alpha/(alpha+B)) * ((B/(alpha+B))^(B/alpha)) - P_max;
B_0 = 10^(-9);
B = fminsearch(fun, B_0);

P_S = (alpha*I_sat)/(log(1+alpha/B));
P = P_S * (1-exp(-(alpha*I)/P_S)) * exp(- (B*I)/P_S);
end

%% Specific Growth Rate
function my = specificGrowthRate (f_area, f_photo, f_temp, N_min, C_min, N, C)
if ((1-N_min/N) <= (1-C_min/C))
    my = f_area * f_photo * f_temp * (1-N_min/N);
else
    my = f_area * f_photo * f_temp * (1-C_min/C);
end
end

%% Model
function state_dot = kelp(N_struct, C_struct, k_A, k_N, k_C, k_dw, N_min, C_min, m_1, m_2, A_O, epsilon, K_X, N_max, J_max, U_065, R_1, T_AR, T_R1, gamma, alpha, I_sat, P_1, T_AP, T_APL,T_APH, U, NO3, T, I, N, C, A)

% Effect of size on growth rate
f_area = m_1*exp(-(A/A_O)^2) + m_2;

%Effect of temperature on growth rate
f_temp = EffectTemp(T);

% Seasonal influence on growth rate
%f_photo = a_1*(1+sign(lambda)*abs(lambda)^0.5) + a_2;
f_photo = 1;

% Frond erosion
ny = (10^(-6)*exp(epsilon*A)) / (1+10^(-6)*(exp(epsilon*A)-1));

% Nitrate uptake rate
J = J_max* (NO3/(K_X + NO3)) * ((N_max-N)/(N_max-N_min)) * (1-exp(-U/U_065));

% Gross photosynthesis
P = grossPhotosynthesis(alpha, I_sat, P_1, T_AP, T_R1, T_APL,T_APH, T, I);

% Temperature dependent respiration
R = R_1 * exp((T_AR/T_R1)-(T_AR/T));

% Carbon exudation
E = 1 - exp(gamma*(C_min-C));

% Specific growth rate
my = specificGrowthRate(f_area, f_photo, f_temp, N_min, C_min, N, C);

% Amount of frond area lost
%     A_lost = A*(C_min - C) / C_struct;

% Structural weight
%     W_s = k_A*A;

% Total dry weight
%     W_d = k_A * (1+ k_N*(N-N_min) + N_min + k_C*(C-C_min) + C_min)*A;

% Total wet weight
%     W_w = k_A * (k_dw^(-1) + k_N*(N-N_min) + N_min + k_C*(C-C_min) + C_min)*A;

% Total carbon content
%     C_total = (C + C_struct) * W_s;

% Total nitrogen content
%     N_total = (N + N_struct) * W_s;

% Rate of change of frond area
% Rate of change in nitrogen reserves
% Rate of change in carbon reserves
state_dot = [(my-ny)*A;
    k_A^(-1)*J-my*(N+N_struct);
    k_A^(-1)*(P*(1-E-R)) - (C+C_struct)*my];

end

%% KelpModel with values from dataset withour offset
% for n = 1:(NumberIterations-1)
%     A = Y(1,n);
%     N = Y(2,n);
%     C = Y(3,n);
%     [state_dot, A_lost, W_s, W_d, W_w, C_total, N_total, my] = kelp(N_struct, C_struct, k_A, k_N, k_C, k_dw, N_min, C_min, m_1, m_2, A_O, epsilon, K_X, N_max, J_max, U_065, R_1, T_AR, T_R1, gamma, alpha, I_sat, P_1, T_AP, T_APL,T_APH, X(2, n), X(3, n), X(1, n), X(4, n), N, C, A);
% %     netCarbonFixed(n+1) = netCarbon(W_s, state_dot(3), n, n+1) + netCarbonFixed(n);
% %     grossFrond(n+1) = totalGross(my, A, n, n+1) + grossFrond(n);
%     Y(:, n+1) = Y(:,n) + state_dot*h;
% end





