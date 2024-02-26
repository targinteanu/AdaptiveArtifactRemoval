%% Goals: 
% We want to visualize the phase space of several examples and use these to
% figure out the appropriate time points and other settings for further
% analysis. 
% To Do: Complete this code and run several times for different examples. 

%% start by loading in the file 
%
% We have folders of .mul files that contain patient SSEP data. There are
% multiple SSEP recordings for each patient pertaining to different test
% dates/times and the right/left median nerve. 
%
% Folders are named with the patient's study ID (e.g. 0055).
% Files are named 'YYYYMMDD-HHMM-<side>-<nerve>.mul' (e.g. '20181119-0845-LT-median.mul'). 
% 
% The data is in a text file. 
% Top row is the channel: C3-Fpz  Cz-Fpz  C3-C4  C4-lt erb  C5s-Fpz  rterb-lt erb
% We care about C3 or C4 depending on whether this is right or left side. 
% The Sampling Interval is also at the top. 

% use either 
folder = uigetdir; 
filelabel = folder.name; 
% or
file = uigetfile; 
filelabel = file.name; 

% Use Ze/Mingfeng's code and/or chatGPT to parse the text document in
% folder or file. 

%% obtain and convert the data 
% We want the data of interest as a time-series x (volts) 
% with associated times t (seconds). 
x = []; % be sure to convert to volts 
dt = []; % Sampling Interval; be sure to convert to seconds 
t = 0:length(x) * dt;

%% calculate and visualize the phase space 

% Define times of interest.
% To Do: adjust these as needed
ti = 0; % starting time; ignore signal before this time
ts = 0.009; % defines sub-PSA as up to 9 ms
tf = 0.4; % ending time; ignore signal after this time
winLen = 10; % moving-average window length 

figure; 

% start by viewing x over t and finding the SSEP peaks (n20, p23)
subplot(221);
measureERP(t,x, .02, .023, [ti, tf], winLen, true); 
title('used for full PSA');
subplot(222); 
measureERP(t,x, .02, .023, [ti, ts], winLen, true); 
title('used for sub-PSA')

% calculate PSA and view
subplot(223);
fullPSA = getPhaseSpace(t,x, [ti, tf], winLen, true,filelabel);
subplot(224);
subPSA  = getPhaseSpace(t,x, [ti, ts], winLen, true,filelabel);