%% Temperature Dataset
t1 = datetime(Envdata.time(1,1),Envdata.time(2,1),Envdata.time(3,1),Envdata.time(4,1),Envdata.time(5,1),Envdata.time(6,1), 'Format','yyyy-MM-dd HH:mm:ss');
t2 = datetime(Envdata.time(1,end),Envdata.time(2,end),Envdata.time(3,end),Envdata.time(4,end),Envdata.time(5,end),Envdata.time(6,end), 'Format','yyyy-MM-dd HH:mm:ss');
time = t1:minutes(15):t2;
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