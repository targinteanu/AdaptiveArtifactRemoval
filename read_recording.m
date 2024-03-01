%% read_recording
function data = read_recording(fileName)
%function data = read_recording(fileName, t)
if isfolder(fileName)
    data = [];
else

% Open the file
fileID = fopen(fileName,'r');

% Read the recording parameters
single_line = fgets(fileID);
formatSpec = '%s %s %s %s %s %s';
C = textscan(single_line, formatSpec);

% Read the channel names (lt erb-rterb should be a single channel but got separated)
single_line = fgets(fileID);
single_line = strrep(single_line, 'lt erb', 'lterb');
formatSpec = '%s %s %s %s %s %s';
C2 = textscan(single_line,formatSpec,1);
C2 = cellfun(@char, C2, 'UniformOutput', false);
% Find the index for C4-C3 or C3-C4
%Lt Median
channel_index = find(strcmp(C2, 'C4-C3'));
%Rt Median
if isempty(channel_index)
    channel_index = find(strcmp(C2, 'C3-C4'));
end

% Read the six channels of data
formatSpec = '%f %f %f %f %f %f';
C3 = cell2mat(textscan(fileID,formatSpec));
%C3 = cell2mat(textscan(fileID,formatSpec,t));
fclose(fileID);

data = C3(:, channel_index);

end

end