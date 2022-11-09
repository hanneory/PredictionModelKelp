% t = tiledlayout(1,3);
% title(t,'Perturbert Temp')
% xlabel(t,'time')
% 
% 
% nexttile
% plot(time, Y_A)
% title('Area')
% nexttile
% plot(time, Y_N)
% title('Nitrogen')
% nexttile
% plot(time, Y_C)
% title('Carbon')


%% Calculate difference

q = zeros((Nsample-1), NumberIterations);
for i = 2:1:Nsample
    q(i, :) = Y_A(i,:) - Y_A(1,:);
end

o = Y_A(2:end, 8275:end);
o_vec_std = zeros(1,size(o, 2));
o_vec_mu = zeros(1,size(o, 2));
o (:,2);
for i = 1:size(o, 2)
    [o_std, o_mu] = std(o(:, i));
     o_vec_std(i) = o_std;
     o_vec_mu(i) = o_mu;
end
figure(22)
size(time(8275:end))
size(o_vec_mu)
plot(time(8275:end), o_vec_mu)
hold on;
plot(time(8275:end), Y_A(1,8275:end))

size(o_vec_mu)

q = zeros((Nsample-1), NumberIterations-8274);

figure(23)
plot(time(8275:end), o_vec_mu-Y_A(1, 8275:end))

% plot(time(8275:end), Y_A(1, 8275:end))
% errorbar(o_vec_mu, o_vec_std);

% time = 8275