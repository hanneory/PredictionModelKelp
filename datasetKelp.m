%% Temperature Dataset
% t1 = datetime(2021,1,1,0,0,0, 'Format','yyyy-MM-dd HH:mm:ss');
% t2 = datetime(2021,12,31,23,59,0, 'Format','yyyy-MM-dd HH:mm:ss');
% time = t1:minutes(30):t2;
% T = randi([278 290], 1, 17520, 'double');
% N = 1:1:17520;
% TemperatureDataset = table(N', time', T', 'VariableNames', {'Number', 'Time', 'Temperature'});

size(Envdata.time)
[numRowsTime, numColsTime] = size(Envdata.time);
time = [];
for t = 1:numColsTime
    time(end+1) = datetime(Envdata.time(1:t), Envdata.time(2:t), Envdata.time(2:t), Envdata.time(3:t))
end


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