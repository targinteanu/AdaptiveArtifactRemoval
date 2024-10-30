%% load data from file(s) 

% find folder 
codedir = cd;
folder = uigetdir; 
cd(folder);
files = dir('*filtered_sigObj_ERP_SNR_Spreadsheet.xlsx');
cd(codedir); 

% choose file(s)
filesel = listdlg("PromptString","Select Files", "SelectionMode","multiple", ...
    "ListString",{files.name});
files = files(filesel);

%% setup 
TTf = []; TTu = [];

%% loop
for f = files'
    f.name

    % parse details from file name 
    [nameParsed, nameParsedLen] = sscanf(f.name, '%c %d %[-] %d %[-] %d %[-] %d');
    subjName = [char(nameParsed(1)),num2str(nameParsed(2))]
    Tdate = nameParsed(4);
    Tstr = nameParsed(6);
    if nameParsedLen >= 8
        Tend = nameParsed(8);
    else
        Tend = nan;
    end

    % timing 
    Tstr = datetime([num2str(Tdate),' | ',num2str(Tstr)], ...
        'InputFormat', 'yyMMdd | HHmmss');
    if ~isnan(Tend)
    Tend = datetime([num2str(Tdate),' | ',num2str(Tend)], ...
        'InputFormat', 'yyMMdd | HHmmss');
    else
        Tend = NaT;
    end

    % read data - filt
    T = readtable(fullfile(f.folder, f.name), "Sheet","Channel B Filt");
    t = T.StimTime;
    t = Tstr + seconds(t);
    if ~isnat(Tend)
        if sum(t > Tend)
            warning(['Stim times may be after end time: ',f.name])
        end
    end
    T.Time  = t;
    T = table2timetable(T);
    TTf = [TTf; T];

    % read data - unfilt
    T = readtable(fullfile(f.folder, f.name), "Sheet","Channel B Unfilt");
    t = T.StimTime;
    t = Tstr + seconds(t);
    if ~isnat(Tend)
        if sum(t > Tend)
            warning(['Stim times may be after end time: ',f.name])
        end
    end
    T.Time  = t;
    T = table2timetable(T);
    TTu = [TTu; T];
end

%% plotting 
figure; 

subplot(2,1,1); 
plot(TTu.Time, TTu.n07amp, '.');
grid on; hold on;
plot(TTu.Time, TTu.n15amp, '.');
ylabel('Amplitude (V)'); legend('N7', 'N10');
title('Unfilt');

subplot(2,1,2); 
plot(TTf.Time, TTf.n07amp, '.');
grid on; hold on;
plot(TTf.Time, TTf.n15amp, '.');
ylabel('Amplitude (V)'); legend('N7', 'N10');
title('Filt');