function data = readCSAT(file)
% 
% Matt Werner (m.werner@vt.edu) - August 8, 2021
% 
% Read a data file provided by the CSAT measuring device.
% 
% Note: This file is expected to be provided in a standard format given by
%       the CSAT device. Any variables not found within the CSAT sample
%       data file (described below) will not be read and exported.
% 
%       Modifications made to the input file may not result in expected
%       results.
% 
%       The sample data file that was used to create this routine was
%       provided by Manu Nimmala of Virginia Tech taken from an experiment
%       in Virginia Tech's drone cage.
% 
% The CSAT is capable of taking high-frequency measurements of the local
% wind speed (m/s) and temperature (C). The header for this file
% contains the measured quantities and their units as...
% 
% "TIMESTAMP","RECORD","Ux","Uy","Uz","Ts","diag_csat","t_hmp","e_hmp"
% "TS","RN","m/s","m/s","m/s","C","unitless","C","kPa"
% 
% where 
%   TS - TIMESTAMP
%   RN - RECORD NUMBER.
% 
% The record number increases linearly from 0 to however many measurements
% were taken; the timestamp marks these measurements with the format
% "YYYY-MM-DD HH:MM:SS.s", where HH is a 24-hour time. Invalid measurements
% are marked by "NAN". Each measurement corresponds to a single line in the
% data file.
% 
%    Inputs:
% 
%              file - The data file provided by the CSAT that is to be
%                     read. This file is expected to provide a standard
%                     (fixed) amount of variables. Any variables not found
%                     within the CSAT sample data file will not be read and
%                     exported.
% 
%    Outputs:
% 
%              data - CSAT measurement data.
% 

%% Checks
% No checks

%% Read the file
% Open the file
fid = fopen(file);

% Count how many lines are in the file
fgetl(fid);fgetl(fid);fgetl(fid);fgetl(fid);
blocks = 0;
while (~feof(fid))
    fgetl(fid);
    blocks = blocks + 1;
end

% Reset the file to the beginning
frewind(fid);

% Create counter and containers in anticipation of data being read
k = 1;
TIME = strings(blocks, 1);
DATA = NaN(blocks, 8);

% Expect the first 4 lines to be the header. Since the variables are fixed,
% also expect the order in which they appear to be fixed. Thus, the file
% can immediately be read, assuming it is not empty.
fgetl(fid);fgetl(fid);fgetl(fid);fgetl(fid);
while (~feof(fid))
    % Read the next line
    line = fgetl(fid);
    % Remove all quotes (") from the line
    line = erase(line, '"');
    % Read each line with values delimited by the comma (,)
    lineData = textscan(line, '%s %f %f %f %f %f %f %f %f', 'Delimiter', ',');
    % Insert time into its own array
    TIME(k, 1) = convertCharsToStrings(lineData{1}{1});
    % Insert each element (not time) into the 
    DATA(k, :) = cell2mat(lineData(2:end));
    % Increment the counter
    k = k + 1;
end

%% Computation
% Organize the data into a table with labelled columns
data = [table(TIME, 'VariableNames', "Time"), array2table(DATA, 'VariableNames', ["Record", "U_x", "U_y", "U_z", "T_s", "diag_csat", "t_hmp", "e_hmp"])];

%% Closing
% Close the file
fclose(fid);