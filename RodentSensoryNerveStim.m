%% load file
load('rat1/amplifier_data/Naive/rat1_N_CMAP_1_230731_113111'); % load the "filt" 2D array 
load('rat1/trigger_data/Naive/rat1_N_CMAP_1_230731_113111_trigger'); % load the "board_adc_data" 1D array 

%% define timing--no response yet from Kiara's team?
Fs = 1; % samples per second 
dt = 1/Fs; % time step (seconds) 
t = 0:dt:length(d_unfilt); % fill this in 

%% processing loaded data into signal object
% "filt" is the neural recording data, which we want to change into "d_unfilt"
% "board_adc_data" is the noise reference data, which we want to change into "g"
d_unfilt = cell(2, length(filt));
for
    i=1:length(filt);
    d_unfilt(1,i)=filt(1,i);
    d_unfilt(2,i)=filt(2,i);
end

g = cell(2, length(board_adc_data));
for
    i=1:length(board_adc_data);
    g(1,i)=board_adc_data(i);
    g(2,i)=board_adc_data(i);
end

chA = buildChannelObj('insert_name_here', 0,0,0, 'Cartesian');
chB = buildChannelObj('insert_name_here', 0,0,0, 'Cartesian');

% (t and g should have the same number of rows as d_unfilt see repmat)

sig = buildSignalObj([], d_unfilt, t, g, Fs, [chA; chB], ...
                     tTrainBnd, tTrainBnd, 2);