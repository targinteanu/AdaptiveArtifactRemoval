%% Initialization
clear; clc;
data_folder = 'C:\Users\targi\Documents\Human_SSEP_Data';
subject_folders = dir(data_folder);
subject_folders = subject_folders(~ismember({subject_folders(:).name},{'.','..'}));
subject_ID = '0117'; %4 digit identifier
fs = 10000; % Hz

%%
ti=4*10; %starting time (ms) [edit number in front of '*10']
tf=90*10; %ending time (ms)
    
% calculate the average recording for the left and right median nerve
data_L = [];
data_R = [];

    % Identify how many SSEP recordings there are
    recordings = dir(fullfile(data_folder, subject_ID));
    recordings = recordings(~ismember({recordings(:).name},{'.','..'}));

    ti=4*10; %starting time (ms) [edit number in front of '*10']
    tf=90*10; %ending time (ms)
    
    % calculate the average recording for the left and right median nerve
    data_L = [];
    data_R = [];
    for recording_pointer = 1:length(recordings)
        if contains(recordings(recording_pointer).name, 'L')
            data_L = [data_L, read_recording(fullfile(data_folder, subject_ID, recordings(recording_pointer).name), tf)];
        else
            data_R = [data_R, read_recording(fullfile(data_folder, subject_ID, recordings(recording_pointer).name), tf)];
        end
    end

    data_L = mean(data_L, 2);
    data_R = mean(data_R, 2);

    if ~isempty(data_L)
    data_L(1:ti,:)=[];
    end
    if ~isempty(data_R)
    data_R(1:ti,:)=[];
    end

    slope_L = gradient(data_L);
    slope_R = gradient(data_R);


%% plot phase space and convex hull
window = 1; %(*0.1ms)
data_L=movmean(data_L,window);
phaseSpace = [data_L,slope_L];
[k,av] = convhull(phaseSpace, 'Simplify', true);

figure
subplot(1,2,1)

plot(phaseSpace(:,1),phaseSpace(:,2),".")
title('Left', ['ID:',subject_ID, ' Area:', num2str(av)]); xlabel('x'); ylabel('dx/dt');
hold on
plot(phaseSpace(k,1),phaseSpace(k,2))
hold off

data_R=movmean(data_R,window);
phaseSpace = [data_R,slope_R];
[k,av] = convhull(phaseSpace, 'Simplify', true);

subplot(1,2,2)

plot(phaseSpace(:,1),phaseSpace(:,2),".")
title('Right', ['ID:',subject_ID, ' Area:', num2str(av)]); xlabel('x'); ylabel('dx/dt');
hold on
plot(phaseSpace(k,1),phaseSpace(k,2))
hold off


av