% TDTPATH = 'TDTMatlabSDK';
% addpath(genpath(TDTPATH));
% data = TDTbin2mat('Rodent SSEP Data/AC5-230830-130841');

folder = uigetdir;
files = dir(folder);
%check if mkdir already exists 
if isfolder([folder,'/','saved_signal_obj']) | isfolder([folder,'\','saved_signal_obj'])
    error('already saved')
else
    mkdir(folder,'saved_signal_obj');
end

codedir = cd;
savedir = [folder,'/','saved_signal_obj/'];

if contains(folder,'rat') % for Kiara files
    for idx = 3:size(files,1)
        cd(folder); %brings loop back to the outer folder
    if ~ isfolder(files(idx,1).name)
        load(files(idx,1).name);
        cd(codedir);
            %% define timing
            Fs = 20000; % samples per second 
            dt = 1/Fs; % time step (seconds) 
            samples=1:length(filt);
            samples = samples-1;
            t = samples/Fs;
            t = t';
            %% processing loaded data into signal object
            % "filt" is the neural recording data, which we want to change into "d_unfilt"
            % "board_adc_data" is the noise reference data, which we want to change into "g"
            d_unfilt = filt;
            d_unfilt = (d_unfilt-32768)*0.0003125; %filt in multiple of volt units
            d_unfilt = d_unfilt';
            
            g = [board_adc_data;board_adc_data];
            g = g';
            
            chA = buildChannelObj('insert_name_here', 0,0,0, 'Cartesian');
            chB = buildChannelObj('insert_name_here', 0,0,0, 'Cartesian');
            
            % (t and g should have the same number of rows as d_unfilt see repmat)
            %% define parameters for filter and training 
            trainfrac = .1;
            nUpdates = 100;
            splIdx = floor(trainfrac*size(t,1));
            tTrainBnd = [t(1), t(splIdx)];
            
            sig = buildSignalObj([], d_unfilt, t, g, Fs, [chA; chB], ...
                                 tTrainBnd, tTrainBnd, 2);
            saveSignalObj([savedir,files(idx,1).name,'_saved'],sig);
    end
    end
cd(codedir);

else %for Mingfeng files 
TDTPATH = 'TDTMatlabSDK';
addpath(genpath(TDTPATH));

PeriodA = 0.8;  % ms - between pulses?
CountA  = 1;    % pulse count? 
AmpA    = 1500; % uA
DurA    = 0.2;  % ms
DelayA  = 0;    % ms
ChanA   = 1;

for index = 3:size(files,1)
    cd(folder); %brings loop back to the outer folder
    fname = files(index,1).name;
    if isfolder([fname,'\',fname])
        data = TDTbin2mat([fname,'\',fname]);
        cd(codedir);


        dta = data.snips.SSEP.data;
        t_stim = data.snips.SSEP.ts;
        fs = data.snips.SSEP.fs;
        chan = data.snips.SSEP.chan; uchan = unique(chan);

        dta_t_chan = cell(2, length(uchan));
        for idx = 1:length(uchan)
            ch = uchan(idx);
            chIdx = chan == ch;
            dta_t_chan{1, idx} =  dta(chIdx,:);
            dta_t_chan{2, idx} = t_stim(chIdx);
        end

        dta = zeros([size(dta_t_chan{1,1}), length(uchan)]);
        t_stim = zeros([length(dta_t_chan{2,1}), length(uchan)]);
        for idx = 1:length(uchan)
            dta(:,:,idx)  = dta_t_chan{1, idx};
            t_stim(:,idx) = dta_t_chan{2, idx};
        end

        t_trl = (1:size(dta, 2))/fs - .3; % ~ -.3 to +1 s
        g_trl = (AmpA/1000)*((t_trl >= 0)&(t_trl < DurA/1000)); % noise reference, mA
        G = repmat(g_trl, size(dta,1), 1, size(dta,3));

        T = zeros(size(dta));
        for idx = 1:length(uchan)
            T(:,:,idx) = t_stim(:,idx) + t_trl;
        end

        chA = buildChannelObj('A', 1, 1,0,'Cartesian'); % right somatosensory
        chB = buildChannelObj('B',-1, 1,0,'Cartesian'); % left somatosensory
        chC = buildChannelObj('C', 1,-1,0,'Cartesian'); % ?
        chD = buildChannelObj('D',-1,-1,0,'Cartesian'); % ?

        sig = buildSignalObj([], d_unfilt, t, g, Fs, [chA; chB; chC; chD], ...
            tTrainBnd, tTrainBnd, 2);

        saveSignalObj([savedir,files(index,1).name,'_saved'],sig);
    end
end

end