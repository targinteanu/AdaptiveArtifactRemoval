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
savedir = [folder,filesep,'saved_signal_obj',filesep];

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

for index = 3:size(files,1)
%    cd(folder); %brings loop back to the outer folder
    fname = files(index,1).name;
    if isfolder([folder,filesep,fname])
%        cd(codedir);

        sig = RodentSSEPtoSig([folder,filesep,fname]);

        saveSignalObj([savedir,files(index,1).name,'_saved'],sig);
    end
end

end