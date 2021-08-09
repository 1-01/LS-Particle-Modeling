function [predictions, t] = narnetWindPrediction(time, measurements, futureTime)
% 
% Matt Werner (m.werner@vt.edu) - August 7, 2021
% 
% Use (N)onlinear (A)uto(r)egression (N)eural N(et)work (NARNET), a
% machine-learning procedure, to attempt prediction of the wind speed and
% (most importantly) wind direction given a set of prior measurements. 
% 
% These measurements must be in an array format with each point
% corresponding to a unique point in time.
% 
% Note: This function requires the deep learning toolbox.
% 
% More information about NARNET can be found at...
% https://www.mathworks.com/help/deeplearning/ref/narnet.html;jsessionid=c2cd272795f05788f97c8fdc71dd
% 
%    Inputs:
% 
%              time - Time elapsed (in minutes) since the recording of the
%                     first measurement. That is, the first measurement
%                     corresponds with a time of 0 and subsequent
%                     measurements correspond with the amount of time, in
%                     minutes, that have passed since that epoch.
%                     Size: n-by-1 (vector)
%                     Units: min (minutes)
% 
%      measurements - Time-history data of local wind properties (speed
%                     and/or direction, etc.) whose future values are to be
%                     predicted.
%                     Size: n-by-m (matrix)
%                     Units: ? (SI)
% 
%        futureTime - The amount of time by which the measurements should
%                     be extended through machine learning. This amount of
%                     time should be given in minutes relative to the time
%                     of the last recorded measurement.
%                     Size: 1-by-1 (sclar)
%                     Units: min (minutes)
% 
%    Outputs:
% 
%       predictions - Future values of the wind measurements predicted by
%                     NARNET deep learning.
%                     Size: k-by-m (vector)
%                     Units: ? (SI)
% 

%% Checks
% Ensure that the first element of time is 0.
assert(numel(time) > 0 & time(1) == 0, "Invalid start time.")
% Ensure that time is flowing forward without repeated values
assert(all(diff(time) > 0), "Invalid time stamps.")
% Ensure that there are many measurements as time logs
assert(size(time, 1) == size(measurements, 1), "Unequal time stamps to measurements.")
% Ensure that the future time is positive
assert(futureTime > 0, "Future time must be positive.")

%% Data Preparation
% Take all data in time and arrange it into a cell that is accepted by
% NARNET algorithms. Data points corresponding to matching instants in time
% must be grouped together into one cell, stacked vertically, such that
% time "flows" horizontally across the cells.
for n = 1:numel(time)
    T{1,n} = measurements(n, :)';
end

% T = simplenar_dataset;

% Determine how much time that a single increase in the index corresponds
% to. Note that this value cannot be zero since the initial check has
% already been passed.
dt = mean(diff(time));
% Calculate how many indices should be estimated to meet the requested
% future time. As this value is an index (integer-valued), interpret the
% result within a ceiling.
maxIndex = ceil(futureTime/dt);

%% Computations (Machine Learning)
% Create a NAR network. Define the feedback delays and size of the hidden
% layers. Also specify the default behaviors of open feedbackMode and the
% Leven-Marquardt training algorithm (trainFcn).
net = narnet(1:2,10, 'open', 'trainlm');
% Prepare the time series data using preparets. This function automatically
% shifts input and target time series by the number of steps needed to fill
% the initial input and layer delay states.
[Xs,Xi,Ai,Ts] = preparets(net,{},{},T);

% A recommended practice is to fully create the network in an open loop,
% and then transform the network to a closed loop for multistep-ahead 
% prediction. Then, the closed-loop network can predict as many future 
% values as you want. If you simulate the neural network in closed-loop
% mode only, the network can perform as many predictions as the number of
% time steps in the input series.
% 
% Train the NAR network. The train function trains the network in an open
% loop (series-parallel architecture), including the validation and testing
% steps.
net = train(net,Xs,Ts,Xi,Ai);
% % Display the trained network.
% view(net)
% Calculate the network output Y, final input states Xf, and final layer
% states Af of the open-loop network from the network input Xs, initial
% input states Xi, and initial layer states Ai.
[Y,Xf,Af] = net(Xs,Xi,Ai);
% Calculate the network performance.
perf = perform(net,Ts,Y)

% To predict the output for the next 20 time steps, first simulate the
% network in closed-loop mode. The final input states Xf and layer states
% Af of the open-loop network net become the initial input states Xic and
% layer states Aic of the closed-loop network netc.
[netc,Xic,Aic] = closeloop(net,Xf,Af);
% % Display the closed-loop network. The network has only one input. In 
% % closed-loop mode, this input connects to the output. A direct delayed
% % output connection replaces the delayed target input.
% view(netc)
% To simulate the network 20 time steps ahead, input an empty cell array of
% length 20. The network requires only the initial conditions given in Xic
% and Aic.
predictions = netc(cell(0,maxIndex),Xic,Aic);
t = time(end)+dt:dt:time(end)+futureTime;