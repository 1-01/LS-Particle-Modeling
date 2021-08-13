function [timePar, dataAvg] = partitionAverage(time, dt, minutesToAverage, data)
% 
% 
% 
% minutesToAverage - Blocks of time averages
% 

%% Checks
% Return the input time and data if the specified block size is 0
if (minutesToAverage == 0)
    timePar = time;
    dataAvg = data;
    return
end

%% Computation
% Determine the amount of indices to skip per block
indexSkip = minutesToAverage/(dt/60);

% Average the partitions
index = 1;
for i = 1:indexSkip:numel(data)
    range = i:i+indexSkip-1;
    if (max(range) > numel(data))
        range = i:numel(data);
    end
    dataAvg(index,1) = mean(data(range));
    timePar(index,1) = mean(time(range));
    index = index + 1;
end