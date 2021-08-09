clc, clear all, close all
% Read the data
data = readCSAT('data/dronecagestest-113-30min.txt');

% Calculate wind direction from measurements
direction = atan2d(data.U_y, data.U_x);
% Extract times from measurements
t = datetime(data.Time, 'InputFormat', 'yyyy-MM-dd HH:mm:ss.S', 'Format', 'yyyy-MM-dd HH:mm:ss.S');

% Calculate the elapsed time
elapsedTime = t - t(1);
elapsedTime.Format = 'mm:ss.S';
% Convert elapsed time to a float
for k = numel(elapsedTime):-1:1
    tmp = char(elapsedTime(k));
    scan = textscan(tmp, '%f %f', 'Delimiter', ':');
    elapsedTimeFloat(k) = scan{1} + scan{2}/60;
end


%%
[prediction, t] = narnetWindPrediction(elapsedTimeFloat, direction, 10);

% oss - small block
% gdx - large block
% gdx - very large block


% Everything else - constant either at 0 or at the top

t = [elapsedTimeFloat; t'];

%
plot(t, [direction; cell2mat(prediction)'])
hold on
plot(elapsedTimeFloat, direction, 'k')
xline(elapsedTimeFloat(end), '-', 'Start NARNET');
hold off
xlabel("Elapsed Time $t$ [minutes]", 'interpreter', 'latex', 'FontSize', 20);
ylabel("Wind Direction $\theta$ [deg.]", 'interpreter', 'latex', 'FontSize', 20);
title("Using NARNET to Predict Wind Direction From CSAT Data", 'Interpreter', 'latex', 'FontSize', 20)
