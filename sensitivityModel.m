% time_A_NO3 = time(3361:6144);
% A_NO3 = Y_A(:, 3361:6144);
%  
% time_N_NO3 = time(5377:7296);
% N_NO3 = Y_N(:, 5377:7296);
%  
% time_C_NO3 = time(577:1824);
% C_NO3 = Y_C(:, 577:1824);
%  
% randomNO3 = struct('time_A_NO3', time_A_NO3, 'A_NO3', A_NO3, 'time_N_NO3', time_N_NO3, 'N_NO3', N_NO3, 'time_C_NO3', time_C_NO3, 'C_NO3', C_NO3);
% save randomNO3 randomNO3

%% Area
% A_NO3 = randomNO3.A_NO3(:,:);
% A_I = randomI.A_I(:,:);
% A_U = randomU.A_U(:,:);
% A_T = randomT.A_T(:,:);
% time_A = randomNO3.time_A_NO3(:,:);
% Y_area = [A_NO3; A_I; A_U; A_T];
% 
% figure(10)
% plot(time_A, Y_area)
% grid on
% title('Comparison of area with rand');


%% Nitrogen
% N_NO3 = randomNO3.N_NO3(:,:);
% N_I = randomI.N_I(:,:);
% N_U = randomU.N_U(:,:);
% N_T = randomT.N_T(:,:);
% time_N = randomNO3.time_N_NO3(:,:);
% Y_nitrogen = [N_NO3; N_I; N_U; N_T];
% 
% figure(11)
% plot(time_N, Y_nitrogen)
% grid on
% title('Comparison of nitrogen with randn');

%% Carbon
% C_NO3 = randomNO3.C_NO3(:,:);
% C_I = randomI.C_I(:,:);
% C_U = randomU.C_U(:,:);
% C_T = randomT.C_T(:,:);
% time_C = randomNO3.time_C_NO3(:,:);
% Y_carbon = [C_NO3; C_I; C_U; C_T];
% 
% figure(12)
% plot(time_C, Y_carbon)
% grid on
% title('Comparison of carbon with randn');

%% 30.april 12:00 - 5329
% A_N03_30april = randomNO3.Y_A(:, 5328);
% A_I_30april = randomI.Y_A(:, 5328);
% A_U_30april = randomU.Y_A(:, 5328);
% A_T_30april = randomT.Y_A(:, 5328);
% Y_A_30april = [A_N03_30april; A_I_30april; A_U_30april; A_T_30april];
% o = zeros(1, size(Y_A_30april, 1)) + 5328;
% figure(14)
% scatter(o', Y_A_30april);
% title('Area scatter with random statevariables');
% 
% 
% C_N03_30april = randomNO3.Y_C(:, 5328);
% C_I_30april = randomI.Y_C(:, 5328);
% C_U_30april = randomU.Y_C(:, 5328);
% C_T_30april = randomT.Y_C(:, 5328);
% Y_C_30april = [C_N03_30april; C_I_30april; C_U_30april; C_T_30april];
% figure(15)
% scatter(o', Y_C_30april);
% title('Carbon scatter with random state variables');
% 
% 
% N_N03_30april = Y_N(:, 5328);
% N_I_30april = Y_N(:, 5328);
% 
% C_N03_30april = Y_C(:, 5328);
% C_I_30april = Y_C(:, 5328);

%% Simulation

%testI0711 = struct('Y_A', Y_A, 'Y_N', Y_N, 'Y_C', Y_C, 'X_T', X_T, 'X_I', X_I, 'X_NO3', X_XNO3, 'X_U', X_U);
%save testI0711 testI0711
q = zeros((Nsample-1), NumberIterations);
for i = 2:1:Nsample
    q(i, :) = Y_A(i,:) - Y_A(1,:);
end

ca = zeros((Nsample-1), NumberIterations);
for i = 2:1:Nsample
    ca(i, :) = C_content(i,:) - C_content(1,:);
end

ni = zeros((Nsample-1), NumberIterations);
for i = 2:1:Nsample
    ni(i, :) = N_content(i,:) - N_content(1,:);
end

% figure(31)
% t = tiledlayout(1,3);
% title(t,'Temperature')
% xlabel(t,'time')
% 
% nexttile
% plot(time(8275:end), pert(8:10,8275:end-1))
% grid on;
% title('Perturbasjon')
% nexttile
% plot(time(8275:end), X_T(8:10,8275:end))
% grid on;
% title('Temperatur med pertubasjon')
% nexttile
% plot(time(8275:end), q(8:10,8275:end))
% grid on;
% title('Avvik areal')

figure(31)
t = tiledlayout(1,2);
title(t,'Temperature')

nexttile
plot(time, pert(8:10,1:end-1))
grid on;
title('Perturbation')
ylabel('Degree Celsius')
nexttile
plot(time, q(8:10,:))
grid on;
title('Deviation area')
ylabel('dm^2')


figure(32)
t = tiledlayout(1,2);
title(t,'Carbon content (fraction of dry weight)')

nexttile
plot(time, pert(8:10,1:end-1))
grid on;
title('Perturbation')
ylabel('Degree Celsius')
nexttile
plot(time, ca(8:10,:))
grid on;
title('Deviation carbon')
ylabel('C(gsw)^{-1}');

figure(33)
t = tiledlayout(1,2);
title(t,'Nitrogen content (fraction of dry weight)')

nexttile
plot(time, pert(8:10,1:end-1))
grid on;
title('Perturbation')
ylabel('Degree Celsius')
nexttile
plot(time, ni(8:10,:))
grid on;
title('Deviation nitrogen')
ylabel('N(gsw)^{-1}');







