function [timeAvg, dataAvg] = partitionAverage(time, dt, minutes, data)
% 
% Matt Werner (m.werner@vt.edu) - August 9, 2021
% 
% Reduce the size of a time-history quantity by segmenting its elements
% into blocks and taking an average. Each block corresponds to a duration
% in time. 
% 
% For very large vectors representing data measured at high frequency, a
% significant reduction in size can be realized if averaging over 30
% seconds (or so) worth of data at a time. This is especially useful if the
% data are evolving very slowly over each measurement.
% 
% Example: A data set is given from an hour's worth of measurements taken
%          at 10 Hz (36,000 data points). To reduce the amount of data, we
%          can partition and average over 30 second intervals such that
%          each data point now corresponds to 30 seconds worth of
%          measurements (averaged over one another). This reduces the data
%          set from 36,000 elements to just 120.
%          (120 = 36000/(30*10))
% 
% Example: Same scenario as before but now averaged over 1 second. This
%          reduces the data set from 36,000 elements to 3,600.
% 
% Inputs:
% 
%              time - The times at which all of the measurements are taken.
%                     Size: n-by-1 (vector)
%                     Units: min (minutes)
% 
%                dt - The timestep separating measurements from one
%                     another. For any given instrument that measures at a
%                     constant frequency, this value is equivalent to 1/Hz,
%                     i.e.
%                                         dt = 1/Hz,
%                     where Hz is the number of samples taken per SECOND.
%                     ~
%                     Note: The units of time and dt are DIFFERENT.
%                     Size: 1-by-1 (scalar)
%                     Units: s (seconds)
% 
%           minutes - Specifies the amount of time (in minutes) for which
%                     each averaging block represents.
%                     Size: 1-by-1 (scalar)
%                     Units: min (minutes)
% 
%              data - The time-series measurements/data points to be
%                     partioned and averaged into equally* distributed
%                     blocks of time.
%                     Size: n-by-1 (vector)
%                     Units: ? (SI)
% 
%                   * All block sizes are of equal size EXCEPT for possibly
%                     (probably) the last one. The last block's size is
%                     adjusted so that it doesn't extend past the available
%                     data. As such, the last block is always either of
%                     equal or less size with/than the rest.
% 
%    Outputs:
% 
%           timeAvg - Times (in minutes) to which all of the block-averaged
%                     data points correspond. 
%                     Note: The times in each block are all averaged
%                           together to produce this quantity. For data
%                           corresponding to equally spaced timesteps, this
%                           quantity is equivalent to the mean and the
%                           median of the times in the block.
%                     Size: m-by-1 (vector)
%                     Units: min (minutes)
% 
%           dataAvg - Data points resulting from partioning the original
%                     data into blocks and averaging over each one. Each
%                     averaged data point corresponds to a single time,
%                     explained above.
%                     Size: m-by-1 (vector)
%                     Units: ? (SI)
% 

%% Checks
% Ensure that time and the data have the same number of elements
assert(numel(time) == numel(data), "Unequal time and data sizes.")

% Return the input time and data if the specified block size is 0
if (minutes == 0)
    timeAvg = time;
    dataAvg = data;
    return
end

%% Computation
% Get the number of data points
N = numel(time);
% Determine the amount of indices to skip per block
blockSize = ceil(minutes/(dt/60));

% Make sure that the block size isn't bigger than the provided data
if (blockSize > N)
    error("Requested block size exceeds data size.")
end

% Calculate how many elements the result will have and allocate memory for
% them
M = ceil(N/blockSize);
dataAvg = NaN(M, 1);
timeAvg = NaN(M, 1);

% Partition the data and average them together
% Note: Potentially confusing with the sizes of variables in the
%       description (n and m) versus how they are about to be defined. 
%       Both n and m end up as N and M, which are the sizes of the input
%       data and output data. There should now be no confusion.
m = 1;
for n = 1:blockSize:N
    % Create the indices for this block
    block = n:n+blockSize;
    % Check if the block extends past the end of the data
    if (block(end) > N)
        % Reduce the block size to match the end of the data. 
        % Note: This block is the only one whose size is unequal from the
        %       rest (it's smaller).
        block = n:N;
    end
    % Assign a single data point corresponding to the average of the data
    % points within this block
    dataAvg(m,1) = mean(data(block));
    timeAvg(m,1) = mean(time(block));
    m = m + 1;
end