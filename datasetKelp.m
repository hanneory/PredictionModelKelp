%Parameters
A_O = 6;
alpha = 3.75*10^(-5);
C_min = 0.01;
C_struct = 0.20;
gamma = 0.5;
epsilon = 0.22;
I_sat = 200;
J_max = 1.4*10^(-4);
k_A = 0.03;
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
R_1 = 2.2*10^(-4); 
%R_1 = 2.785*10^(-4);
R_2 = 5.429*10^(-4);
T_R1 = 285;
T_R2 = 290;
T_AP = 1694.4;
T_APH = 25924;
T_APL = 27774;
T_AR = 11033;
U_065 = 0.03;
K_X = 4;

% A_O = 8; 
% alpha = 1.204*10^(-5); 
% C_min = 0.01; 
% C_struct = 0.20; 
% gamma = 0.5; 
% epsilon = 0.22; 
% I_sat = 90; 
% J_max = 2.8*10^(-5); 
% k_A = 0.027;
% k_dw = 0.0685; 
% k_C = 2.1213; 
% k_N = 2.8929; 
% m_1 = 0.1085; 
% m_2 = 0.03; 
% my_max = 0.1823;
% N_min = 0.0088; 
% N_max = 0.0216; 
% N_struct = 0.0121; 
% P_1 = 1.44*10^(-3); 
% P_2 = 1.44*10^(-3); 
% a_1 = 0.85; %
% a_2 = 0.3; %
% R_1 = 2.2*10^(-4); 
% R_2 = 5.429*10^(-4); 
% T_R1 = 285; 
% T_R2 = 290; 
% T_AP = 1694.4; 
% T_APH = 25924; 
% T_APL = 27774; 
% T_AR = 6200; 
% U_065 = 0.03; 
% K_X = 56; 

Nsample = 200;
NumberIterations = size(Envdata.time, 2);

%% Temperature Dataset
t1 = datetime(Envdata.time(1,1),Envdata.time(2,1),Envdata.time(3,1),Envdata.time(4,1),Envdata.time(5,1),Envdata.time(6,1), 'Format','yyyy-MM-dd HH:mm:ss');
t2 = datetime(Envdata.time(1,end),Envdata.time(2,end),Envdata.time(3,end),Envdata.time(4,end),Envdata.time(5,end),Envdata.time(6,end), 'Format','yyyy-MM-dd HH:mm:ss');
time = t1:minutes(15):t2;

%% Initial value A_0
% A_0 = zeros(Nsample, 1) + 0.01;
%  for i=1:Nsample
%      A_0(i) = addRandomAddition(A_0(i), 1);
% end

A_0 = [startareas(:,1); startareas(:,2); startareas(:,3); startareas(:,4); startareas(:,5)];
[A0_std, A0_mu] = std(A_0(:, 1));
% for i = 1:Nsample
%     if A_0(i) > 50 
%         A_0(i, 1) = 20;
%     end
% end


%% Statevariable
% X = zeros(4, NumberIterations);
% for n = 1:NumberIterations
%     T = Envdata.T(:,:,2,n) + 273.15;
%     U = 0.06;
%     X_NO3 = Envdata.NO3(:,:,2,n);
%     I = Envdata.PAR(:,:,2,n);
%     X(:, n) = [T U X_NO3 I];
% end

X_T = zeros(Nsample, NumberIterations);
X_U = zeros(Nsample, NumberIterations);
X_XNO3 = zeros(Nsample, NumberIterations);
X_I = zeros(Nsample, NumberIterations);
pert = zeros(Nsample, NumberIterations);
for n = 1:NumberIterations
    X_T(:, n) = Envdata.T(:,:,2,n) ;
    X_U(:, n) = 0.06;
    X_XNO3(:, n) = Envdata.NO3(:,:,2,n);
    X_I(:, n) = Envdata.PAR(:,:,2,n);
end

%% Standard variation and mean
[T_std, T_mu] = std(X_T(1, :));
[U_std, U_mu] = std(X_U(1, :));
[NO3_std, NO3_mu] = std(X_XNO3(1, :));
[I_std, I_mu] = std(X_I(1, :));

[A_std, A_mu] = std(A_0);

%% Randomize state variable

for i = 2:Nsample
     for j = 1:NumberIterations
          pert(i,j+1) = gaussMarkov(pert(i,j), 0.4, 15/1440, T_mu*0.1, 0);
          X_T(i, j) = X_T(i, j) + pert(i,j);
      end
 end

% [pert, X_T] = randomizeStateVariable(X_T, 0.2, 15/1440, T_mu*0.1, 0, Nsample, NumberIterations);
% function [pert, pertMatrix] = randomizeStateVariable(matrix, beta, dt, sigma, mean, Nsample, NumberIterations)
% pert = zeros(Nsample, NumberIterations);
% for i = 2:Nsample
%      for j = 1:NumberIterations
%           pert(i,j+1) = gaussMarkov(pert(i,j), beta, dt, sigma, mean);
%           matrix(i, j) = matrix(i, j) + pert(i,j);
%           if (matrix(i, j) < 0)
%               matrix(i, j) = 0;
%           end
%       end
% end
% pertMatrix = matrix;
% end

%for i = 2:Nsample
%   for j = 1:NumberIterations
%        X_I(i,j) = addRandomAddition(X_I(i,j), I_mu*0.8);
%   end
%end

%% Gauss Markov noise
function x_next = gaussMarkov(x, beta, dt, sigma, mean)
f = exp(-beta*dt);
r = randn*sigma + mean;
x_next = f*x + sqrt(1-f^2)*r;
end

%% White noise
function y = addRandomAddition(x, sigma)
%y = x + randn*sigma;
y = x + 2;
if (y)<0 
    y = 0;
end
end

% T = randi([278 290], 1, 17520, 'double');
% N = 1:1:17520;
% TemperatureDataset = table(N', time', T', 'VariableNames', {'Number', 'Time', 'Temperature'});


% %% Currents Dataset
% U = randi([0 5], 1, 17520)/100;
% CurrentsDataset = table(N', time', U', 'VariableNames', {'Number', 'Time', 'Currentspeed'});
% 
% %% Nutrient Concentration Dataset
% X = randi([0 10], 1, 17520);
% NutrientDataset = table(N', time', X', 'VariableNames', {'Number', 'Time', 'Nutrient concentration'});
% 
% %% Irradiance Dataset
% I = randi([150 250], 1, 17520);
% IrradianceDataset = table(N', time', I', 'VariableNames', {'Number', 'Time', 'Irradiance'});