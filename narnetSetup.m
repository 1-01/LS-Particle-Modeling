clc, clear, close all

% Read the data
data = readCSAT('data/dronecagestest-116-all.txt');

% Read in measurements
% measurement = atan2d(data.U_y, data.U_x);
% vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
% TEMP TEMP TEMP TEMP TEMP TEMP TEMP TEMP TEMP TEMP TEMP TEMP TEMP TEMP
measurement = data.U_x;
% TEMP TEMP TEMP TEMP TEMP TEMP TEMP TEMP TEMP TEMP TEMP TEMP TEMP TEMP 
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

% Extract times from measurements
timeStamps = datetime(data.Time, 'InputFormat', 'yyyy-MM-dd HH:mm:ss.S', 'Format', 'yyyy-MM-dd HH:mm:ss.S');

% Calculate the elapsed time
elapsedTimeStamps = timeStamps - timeStamps(1);
elapsedTimeStamps.Format = 'mm:ss.S';
% Convert elapsed time from 'duration' to 'double'
for k = numel(elapsedTimeStamps):-1:1
    scan = textscan(char(elapsedTimeStamps(k)), '%f %f', 'Delimiter', ':');
    elapsedTime(k,1) = scan{1} + scan{2}/60;
end

% Average data across a chosen interval of time (minutes) to significantly
% reduce the size of the training data.
% Note: If not wanting to create any averaging blocks for combining data,
%       specify the minutes to average together as 0.
Hz = 1/(mean(diff(elapsedTime))*60);
[elapsedTimeBlock, measurementBlock] = partitionAverage(elapsedTime, 1/Hz, 0.5, measurement);

%% Training
% Select a subset of the data to use for training (and use the remaining
% data to compare against to the prediction to see how well it turned out).
elapsedTimeTraining = elapsedTimeBlock(1:60);
measurementTraining = measurementBlock(1:60);
[measurementPrediction, elapsedTimePrediction] = narnetWindPrediction(elapsedTimeTraining, measurementTraining, 30);

% % Combine the training time and predicted time
% tNARNET = [elapsedTimeTraining; elapsedTimePrediction];
% xNARNET = [measurementTraining; measurementPrediction];

% Plot the complete data, training data, and predicted values
hold on
plot(elapsedTimeBlock, measurementBlock, 'color', '#0072BD')
plot(elapsedTimeTraining, measurementTraining, 'k')
xline(elapsedTimeTraining(end), '-', 'Start NARNET');
plot(elapsedTimePrediction, measurementPrediction, 'color', '#D95319')
hold off
xlabel("Elapsed Time $t$ [minutes]", 'interpreter', 'latex', 'FontSize', 20);
ylabel("Wind Direction $\theta$ [deg.]", 'interpreter', 'latex', 'FontSize', 20);
title("Using NARNET to Predict Wind Direction From CSAT Data", 'Interpreter', 'latex', 'FontSize', 20)
h = legend("True", "Training", "Predicted", 'interpreter', 'latex', 'fontsize', 20);
title(h, "Data Types", 'FontSize', 20)