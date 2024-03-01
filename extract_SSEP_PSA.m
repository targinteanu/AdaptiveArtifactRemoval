%{
This code extracts the PSA of SSEP for each subject
%}
%% Initialization
clear; clc;
data_folder = 'C:\Users\targi\Documents\Human_SSEP_Data';
subject_folders = dir(data_folder);
subject_folders = subject_folders(~ismember({subject_folders(:).name},{'.','..'}));
subject_ID = {subject_folders(:).name}';
%subject_ID = {'0082', '0092', '0094', '0105', '0102', '0112', '0123', '0126'};
% subject_ID = {'0078', '0082', '0083', '0084', '0085', '0090', '0092', '0094', '0099', '0103', '0105', '0106', '0109', ...
%     '0120', '0121', '0133', '0134', ...
%     '0101', '0102', '0112', '0116', '0119', '0123', '0126'};
n_subject = length(subject_ID);

output_folder = 'SSEP PSA';
fs = 10000; % Hz

%%
PSA = nan(n_subject, 2);
for subject_pointer = 1:n_subject
    subject_name = subject_ID{subject_pointer};

    % Identify how many SSEP recordings there are
    recordings = dir(fullfile(data_folder, subject_name));
    recordings = recordings(~ismember({recordings(:).name},{'.','..'}));
    
    % calculate the average recording for the left and right median nerve
    data_L = [];
    data_R = [];
    for recording_pointer = 1:length(recordings)
        if contains(recordings(recording_pointer).name, 'L')
            data_L = [data_L, read_recording(fullfile(data_folder, subject_name, recordings(recording_pointer).name))];
        else
            data_R = [data_R, read_recording(fullfile(data_folder, subject_name, recordings(recording_pointer).name))];
        end
    end
    data_L = mean(data_L, 2);
    data_R = mean(data_R, 2);   

    slope_L = gradient(data_L);
    slope_R = gradient(data_R);

    if ~isempty(data_L)
        [~, PSA_L] = convhull(slope_L, data_L);
    end

    if ~isempty(data_R)    
        [~, PSA_R] = convhull(slope_R, data_R);
    end
    
    %if ~isempty(data_L) | ~isempty(data_R)
    PSA(subject_pointer, :) = [PSA_L, PSA_R];
    %end
end

PSA = table(subject_ID, PSA);
save(fullfile(output_folder, 'PSA'), 'PSA');

