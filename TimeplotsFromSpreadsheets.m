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
TTf = {}; TTu = {}; TTtime = [];

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
    TTtime = [TTtime; Tstr];

    % read data
    T = getTimeTable(f, "Channel B Filt", Tstr, Tend);
    TTf = [TTf; {T}];
    T = getTimeTable(f, "Channel B Unfilt", Tstr, Tend);
    TTu = [TTu; {T}];
end

%% plotting 
% consider first table = baseline, second = injury 
TTsel = TTf; % pick which one to plot 
winsize = 100; % samples 
varplt = {'n15amp',  'n07amp'}; 
varnam = {'N10',     'N7'};
varmkr = {'.',       '.'};
varclr = {'#D95319', '#0072BD'}; % red, blue
CAtime = TTtime(2); REtime = TTtime(3);

% summary stats helpers 
avgTbl = @(T) mean(T, 'omitnan');
stdTbl = @(T) varfun(@(x)std(x,'omitnan'),T); 
numTbl = @(T) varfun(@(x)sum(~isnan(x)),T);
ci95rng = @(T) tinv(.975, numTbl(T).Variables-1) .* stdTbl(T)./sqrt(numTbl(T).Variables);

% summary stats 
BLavg = avgTbl(TTsel{1}); BLstd = stdTbl(TTsel{1}); BLci = ci95rng(TTsel{1});
CAavg = avgTbl(TTsel{2}); CAstd = stdTbl(TTsel{2}); CAci = ci95rng(TTsel{2});

% concatenate table 
TTcat = [];
for TTi = TTsel'
    TTcat = [TTcat; TTi{:}];
end
TTcat.Time = TTcat.Time - CAtime; REtime = REtime - CAtime;
REtime_lbl = round(minutes(REtime));

figure; 

% plot summary stats 
for v = 1:length(varplt)
    x = TTcat.Time; x = [min(x), max(x)]; xP = [x(1), x(1), x(2), x(2)];
    y = BLavg.(varplt{v}); 
    yP = BLci.(['Fun_',varplt{v}]); 
    yP = [y-yP, y+yP, y+yP, y-yP];
    y = [y, y];
    plot(x,y, '--', 'Color',varclr{v}, 'LineWidth',2); hold on;
    patch('XData',xP, 'YData',yP, 'FaceColor',varclr{v}, ...
        'FaceAlpha',.1, 'EdgeColor','none');
    %{
    y = CAavg.(varplt{v}); 
    yP = CAci.(['Fun_',varplt{v}]); 
    yP = [y-yP, y+yP, y+yP, y-yP];
    y = [y, y];
    plot(x,y, '-.', 'Color',varclr{v}, 'LineWidth',2); hold on;
    patch('XData',xP, 'YData',yP, 'FaceColor',varclr{v}, 'FaceAlpha',.1);
    %}
end

% plot data points and moving mean 
for v = 1:length(varplt)
    x = TTcat.Time; y = TTcat.(varplt{v});
    %plot(x, y, varmkr{v}, 'Color', varclr{v}); hold on;
    ysel = ~isnan(y); 
    x = x(ysel); y = y(ysel); 
    yavg = movmean(y,winsize); ystd = movstd(y,winsize);
    yci = tinv([.025,.975], winsize-1) .* ystd/sqrt(winsize);
    %errorbar(x, yavg, yci(:,1), yci(:,2), '-', 'Color',varclr{v}, 'LineWidth',1.5);
    %%{
    yci = yci + yavg;
    plot(x, yavg, '-', 'Color',varclr{v}, 'LineWidth',1.5); hold on;
    %plot(x, yci, '--', 'Color',varclr{v}, 'LineWidth',1);
    xP = [x; flipud(x)]; yP = [yci(:,1); flipud(yci(:,2))];
    patch('XData', xP, 'YData', yP, ...
        'FaceColor',varclr{v}, 'FaceAlpha',.2, 'EdgeColor','none');
    %}
end

% plot start of resusc
yrng = ylim; 
plot(seconds(0)*ones(size(yrng)), yrng, '--k', 'LineWidth',1.5);
plot(    REtime*ones(size(yrng)), yrng,  '-k', 'LineWidth',1.5);

% labels 
grid on;
set(gca, 'FontSize',16);
xlabel('time from CA onset', 'FontSize',20); 
ylabel('Amplitude (V)', 'FontSize',20);
title([num2str(REtime_lbl),'-Minute CA Recovery'], 'FontSize',28);
legend('N10 Baseline Avg', 'N10 Baseline 95% C.I.', 'N17 Baseline Avg', 'N17 Baseline 95% C.I.', ...
    'N10 Moving Average', 'N10 Moving 95% C.I.', 'N7 Moving Average', 'N7 Moving 95% C.I.', ...
    'CA onset', 'Resuscitation onset', ...
    'Location','westoutside', 'FontSize',18);

%% helpers 

function T = getTimeTable(f, sheetname, Tstr, Tend)
    T = readtable(fullfile(f.folder, f.name), "Sheet",sheetname);
    t = T.StimTime;
    t = Tstr + seconds(t);
    if ~isnat(Tend)
        if sum(t > Tend)
            warning(['Stim times may be after end time: ',f.name])
        end
    end
    T.Time  = t;
    T = table2timetable(T);
end