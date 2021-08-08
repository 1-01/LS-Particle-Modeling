function averagedData = readSODAR(file)
% 
% Matt Werner (m.werner@vt.edu) - July 26, 2021
% 
% Read a data file (.dat extension) provided by the SODAR measuring system.
% 
% Note: This file is expected to be provided in a standard format given by
%       the SODAR device. Any variables not found within the SODAR user
%       manual (link provided below) will not be read and exported.
% 
%       Modifications made to the input file may not result in expected
%       results. 
% 
%       The sample data file that was used to create this routine was
%       provided by Dr. Gonzalez-Rocha of the University of California,
%       Riverside (UCR) taken from an experiment in 2020 at KEAS.
% 
% The SODAR is capable of taking periodic measurements of atmospheric
% conditions (wind speed and direction, standard deviations, covariance,
% etc.)
% 
% For more information about the SODAR device, refer to its manual at...
%  https://www.remtechinc.com/sites/default/files/2017-05/USER_MANUAL_0.pdf
% 
%    Inputs:
% 
%              file - The data file provided by the SODAR that is to be
%                     read. This file is expected to provide a standard
%                     (fixed) amount of variables. Any variables not found
%                     within the SODAR user manual will not be read and
%                     exported.
% 
%    Outputs:
% 
%      averagedData - SODAR measurement data
% 

%% Checks
% No checks

%% Read the file
% Open the file
fid = fopen(file);

% Determine how many measurements (in time) were made. The end of a
% measurement period is marked by the special character $ in the first
% and only column of a new line.
blocks = 0;
while (~feof(fid))
    % Read the current line
    line = fgetl(fid);
    % Search for $
    if (strcmp(line, '$'))
        % Increment the count
        blocks = blocks + 1;
    end
end

% Restart at the beginning of the file
frewind(fid);

% Read the file for each block.
keyphr = "ALT    CT SPEED   DIR     W    SW    SU    SV  U'V'  U'W'  ETAM";
for block = 1:blocks
    while true
        % Get this line
        line = fgetl(fid);
        % Read the block looking for the line containing the key phrase
        % "ALT    CT SPEED   DIR     W    SW    SU    SV  U'V'  U'W'  ETAM"
        % with exactly this particular spacing. In this case, the following
        % line will be blank and the line following that will contain the
        % measured data in the same appearing order across, decreasing in
        % altitude as the line number continues to advance until hitting
        % the block escape character $.
        if (contains(line, keyphr))
            % Set the row counter to 1
            row = 1;
            % Skip 2 lines to get to the data
            fgetl(fid); line = fgetl(fid);
            % Read the data
            while (~strcmp(line, '$'))
                % Values are delimited by at least one space character ' '.
                % Do NOT add the name/value pair ('Delimiter', ' ') to the
                % options of textscan, as this option will literally
                % interpret ' ' in the data file as NaN.
                %
                % Convert the textscan output (a cell) to an array since
                % all entries are expected to be numeric. Also transpose
                % the result to place it as a column vector
                data = cell2mat(textscan(line, '%f'))';
                % Change all -9999 entries (default value of no
                % measurement) to NaN
                data(data == -9999) = nan;
                % Add data in appearing order of the key phrase
                DATA(row, :, block) = data;
                % Read the next line and increment the row counter
                line = fgetl(fid);
                row = row + 1;
            end
            % Break out of the block since the escape character $ contained
            % within every block has been reached
            % reached
            break
        end
    end
end

%% Computations
% Average the data in time (blocks), ignoring NaN (-9999) entries in the
% averaging
averagedData = mean(DATA, 3, 'omitnan');
% Place into a table with the original keyphrase heading
headings = textscan(keyphr, '%s');
averagedData = array2table(averagedData, 'VariableNames', headings{1});

%% Closing
% Close the file
fclose(fid);