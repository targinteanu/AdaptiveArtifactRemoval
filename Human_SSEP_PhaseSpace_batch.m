%% setup before loop 
codedir = cd;
folder = uigetdir; 
files = dir(folder);

PStable = table; % table is initiated here.  
% Currently, the table has columns of Phase Space Area (fullPSA) and Phase
% Space sub-Area (subPSA), and rows corresponding to each recording of each
% subject. 
% To Do: consider alternative ways to organize the table. For example, it
% could be cleaner to have each subject as a row, and different sheets or
% columns corresponding to different recordings such as right vs left side.

%% main loop
for f = 1:length(files)
    subjname = files(f).name; 
    if isfolder([folder,filesep,subjname])
        %
        % [folder,filesep,subjname] should be a path to a folder full of .mul
        % files that contain patient SSEP data.
        %
        % subjname is the subject's study ID (e.g. 0055). 
        %
        % Files are named 'YYYYMMDD-HHMM-<side>-<nerve>.mul' (e.g.
        % '20181119-0845-LT-median.mul'). 
        %
        subjfiles = dir([folder,filesep,subjname]); % all files for this subject. 

        for ff = 1:length(subjfiles)
            subjfilename = subjfiles(ff).name;

            %% load in the data 
            % The data is in a text file. 
            % Top row is the channel: C3-Fpz  Cz-Fpz  C3-C4  C4-lt erb  C5s-Fpz  rterb-lt erb
            % We care about C3 or C4 depending on whether this is right or
            % left side. 
            % The Sampling Interval is also at the top. 

            % To Do: pull in the data of interest as a time-series x (volts) 
            % with associated times t (seconds). Use Ze/Mingfeng's code
            % and/or chatGPT to parse the text document. 
            x = []; % be sure to convert to volts 
            dt = []; % Sampling Interval; be sure to convert to seconds 
            t = 0:length(x) * dt;

            %% measure and record phase space areas

            % Define times of interest. 
            % To Do: adjust these as needed 
            ti = 0; % starting time; ignore signal before this time
            ts = 0.009; % defines sub-PSA as up to 9 ms
            tf = 0.4; % ending time; ignore signal after this time

            % calculate PSA 
            fullPSA = getPhaseSpace(t,x,[ti, tf],10);
            subPSA  = getPhaseSpace(t,x,[ti, ts],10);

            % populate table 
            % Feel free to change the organization as desired. 
            % Here is one way to do it: 
            tblRow = table; % start a new row 
            tblRow.StudyID = subjname; % record subject study ID 
            tblRow.Recording = subjfilename; % record details about what recording this is
            tblRow.fullPSA = PSA; 
            tblRow.subPSA = subPSA; 
            PSAtable = [PSAtable; tblRow]; % stack row onto main table
        end
    end
end

%% save results to file 
% To Do: save PSAtable(s) to excel (or mat?) file using writetable() 
doc writetable