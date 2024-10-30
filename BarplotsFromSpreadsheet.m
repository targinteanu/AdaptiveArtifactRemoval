keys = cell(2,1); % [AC; CA] 
keys{1} = readtable("Key.xlsx", "Sheet","AC", "ReadRowNames",true);
keys{2} = readtable("Key.xlsx", "Sheet","CA", "ReadRowNames",true);

sheets = cell(2,2,2);
    % (1,:,:) = AC; (2,:,:) = CA
    % (:,1,:) = unfilt; (:,2,:) = filt 
    % (:,:,1) = ERP; (:,:,2) = PSA

sheets{1,1,1} = readtable("ERP_SNR_Spreadsheet (2).xlsx", "Sheet","Channel B Unfilt", "ReadRowNames",true);
sheets{1,2,1} = readtable("ERP_SNR_Spreadsheet (2).xlsx", "Sheet","Channel B Filt", "ReadRowNames",true);
sheets{1,1,2} = readtable("PSA_Spreadsheet (4).xlsx", "Sheet","Channel B Unfilt", "ReadRowNames",true);
sheets{1,2,2} = readtable("PSA_Spreadsheet (4).xlsx", "Sheet","Channel B Filt", "ReadRowNames",true);

sheets{2,1,1} = readtable("ERP_SNR_Spreadsheet - CA (2).xlsx", "Sheet","Channel B Unfilt", "ReadRowNames",true);
sheets{2,2,1} = readtable("ERP_SNR_Spreadsheet - CA (2).xlsx", "Sheet","Channel B Filt", "ReadRowNames",true);
sheets{2,1,2} = readtable("PSA_Spreadsheet - CA (2).xlsx", "Sheet","Channel B Unfilt", "ReadRowNames",true);
sheets{2,2,2} = readtable("PSA_Spreadsheet - CA (2).xlsx", "Sheet","Channel B Filt", "ReadRowNames",true);

% append sheets with identifier from keys 
for l = 1:size(sheets,1)
    for m = 1:size(sheets,2)
        for n = 1:size(sheets,3)
            thisSheet = sheets{l,m,n};
            thisKey = keys{l};
            thisSheet = [thisKey, thisSheet];
            sheets{l,m,n} = thisSheet;
        end
    end
end

summarySheetsAC = cell(2,2); 
summarySheetsCA = cell(2,2);
    % (1,:) = unfilt; (2,:) = filt 
    % (:,1) = ERP; (:,2) = PSA
for m = 1:size(sheets,2)
    for n = 1:size(sheets,3)

        % AC: get summary stats 
        ACsheet = sheets{1,m,n};
        ACsummarySheet = table;
        LVL = ACsheet(:,1); %ACsheet = ACsheet(:,2:end);
        LVLu = unique(LVL{:,1})';
        for thisLVL = LVLu
            thisIdx = LVL{:,1} == thisLVL;
            thisTbl = ACsheet(thisIdx,:);
            thisNewRow = thisTbl(1,:); 
            for v = 2:width(thisTbl)
                vname = thisTbl.Properties.VariableNames{v};
                vvals = thisTbl{:,v};
                if contains(vname, 'Mean', 'IgnoreCase',true)
                    vval = mean(vvals);
                elseif (contains(vname, 'SD', 'IgnoreCase',true) | contains(vname, 'std', 'IgnoreCase',true))
                    vval = rms(vvals);
                elseif contains(vname, 'num', 'IgnoreCase',true)
                    vval = sum(vvals);
                else
                    vval = nan;
                end
                % fix units 
                if contains(vname, 'amp')
                    vval = vval*1e6; % Volts -> uV
                end
                if contains(vname, 'lat')
                    vval = vval*1e3; % seconds -> ms
                end
                thisNewRow{1,v} = vval;
            end
            ACsummarySheet = [ACsummarySheet; thisNewRow];
        end
        summarySheetsAC{m,n} = ACsummarySheet;

        % CA: ? 
        CAsheet = sheets{2,m,n};
        % fix units 
        for v = 1:width(CAsheet)
            vname = CAsheet.Properties.VariableNames{v};
            if contains(vname, 'amp')
                CAsheet{:,v} = CAsheet{:,v}*1e6; % Volts -> uV
            end
            if contains(vname, 'lat')
                CAsheet{:,v} = CAsheet{:,v}*1e3; % seconds -> ms
            end
        end
        summarySheetsCA{m,n} = CAsheet;

    end
end

%% AC PSA Plots 
figure('Units','normalized', 'Position',[.1,.1,.8,.8]);
H = size(summarySheetsAC, 1);  
rNames = {'Unfiltered'; 'Filtered'}; 
vDesired = {{'PSAsub', 'PSA'}, {'PSA09', 'PSA30'}}; 
vNewLbl = {{'5-9ms', '5-30ms'}, {'0-9ms', '0-30ms'}};
P = length(vDesired);
for r = 1:H % [unfilt; filt]
    psaTbl = summarySheetsAC{r,2};
    for p = 1:P
        subplot(P, H, (p-1)*H + r);
        Yval = []; Yerr = [];
        for vname = vDesired{p}
            vname = vname{:};
                yval = eval(['psaTbl.',vname,'mean']); 
                yerr = eval(['psaTbl.',vname,'std']); 
                ynum = eval(['psaTbl.',vname,'num']);
                ysem = yerr./sqrt(ynum);
                yerr = ysem.*tinv(.975,ynum);
            Yval = [Yval, yval]; Yerr = [Yerr, yerr];
        end
        plotBarWithErr(Yval, Yerr, 1.5);
        % label plot 
        X = psaTbl(:,1); 
        xlabel(X.Properties.VariableNames{1}); 
        X = X{:,1}; X = arrayfun(@(x) [num2str(x),'%'], X, 'UniformOutput',false);
        xticklabels(X);
        ylabel('PSA'); title(rNames{r});
        %lgd = [vNewLbl{p}, '±1 S.D.'];
        lgd = [vNewLbl{p}, '95% C.I.'];
        legend(lgd, 'Location','eastoutside');
    end
end

%% AC ERP Plots 
figure('Units','normalized', 'Position',[.1,.1,.8,.8]);
H = size(summarySheetsAC, 1);  
rNames = {'Unfiltered'; 'Filtered'}; 
vDesired = {'n07', 'n15'}; vNewLbl = {'N7', 'N10'};
pDesired = {'amp', 'lat', 'num'}; pNewLbl = {'Amplitude (\muV)', 'Latency (ms)', 'Amount'};
P = length(pDesired);
for r = 1:H % [unfilt; filt]
    erpTbl = summarySheetsAC{r,1};
    for p = 1:P
        subplot(P, H, (p-1)*H + r);
        Yval = []; Yerr = [];
        for v = vDesired
            vname = [v{:},pDesired{p}];
            if p == P
                yval = eval(['erpTbl.',vname]); 
                yerr = 0*yval;
            else
                yval = eval(['erpTbl.',vname,'Mean']); 
                yerr = eval(['erpTbl.',vname,'SD']);   
                ynum = eval(['erpTbl.',v{:},'num']);
                ysem = yerr./sqrt(ynum);
                yerr = ysem.*tinv(.975,ynum);
            end
            Yval = [Yval, yval]; Yerr = [Yerr, yerr];
        end
        plotBarWithErr(Yval, Yerr, 1.5);
        % label plot 
        X = erpTbl(:,1); 
        xlabel(X.Properties.VariableNames{1}); 
        X = X{:,1}; X = arrayfun(@(x) [num2str(x),'%'], X, 'UniformOutput',false);
        xticklabels(X);
        ylabel(pNewLbl{p}); title(rNames{r});
        %lgd = [vNewLbl, '±1 S.D.'];
        lgd = [vNewLbl, '95% C.I.'];
        legend(lgd, 'Location','eastoutside');
    end
end

%% CA PSA Plots 
figure('Units','normalized', 'Position',[.1,.1,.8,.8]);
H = size(summarySheetsAC, 1);  
rNames = {'Unfiltered'; 'Filtered'}; 
vDesired = {'PSAsub', 'PSA', 'PSA09', 'PSA30'}; 
vNewLbl = {'5-9ms', '5-30ms', '0-9ms', '0-30ms'};
P = length(vDesired);
for r = 1:H % [unfilt; filt]
    psaTbl = summarySheetsCA{r,2};
    grp = psaTbl{:,1}; grpu = unique(grp)';
    for p = 1:P
        subplot(P, H, (p-1)*H + r);
        pname = vDesired{p}; 
        pIdx = contains(psaTbl.Properties.VariableNames, pname);
        pTbl = psaTbl(:,pIdx);
        Yval = []; Yerr = [];
        for gname = grpu
            gname = gname{:};
            gIdx = strcmp(grp, gname);
            gTbl = pTbl(gIdx,:);
                yval = eval(['gTbl.',pname,'mean']); 
                yerr = eval(['gTbl.',pname,'std']);   
            Yval = [Yval, yval]; Yerr = [Yerr, yerr];
        end
        plotBarWithErr(Yval, Yerr, 1);
        % label plot 
        X = psaTbl(:,2); 
        xlabel(X.Properties.VariableNames{1}); 
        X = X{:,1}; 
        xticklabels(X);
        ylabel(['PSA ',vNewLbl{p}]); title(rNames{r});
        lgd = [grpu, '±1 S.D.'];
        legend(lgd, 'Location','eastoutside');
    end
end

%% CA ERP Plots 
figure('Units','normalized', 'Position',[.05,.05,.9,.9]);
H = size(summarySheetsAC, 1);  
rNames = {'Unfiltered'; 'Filtered'}; 
vDesired = {'n07', 'n15'}; vNewLbl = {'N7', 'N10'};
pDesired = {'amp', 'lat', 'num'}; pNewLbl = {'Amplitude (\muV)', 'Latency (ms)', 'Amount'};
P = length(pDesired); V = length(vDesired);
for r = 1:H % [unfilt; filt]
    erpTbl = summarySheetsCA{r,1};
    grp = psaTbl{:,1}; grpu = unique(grp)';
    for p = 1:P
    for v = 1:V
        subplot(P, H+V, (p-1)*(H+V) + (v-1)*H + r); 
        vname = [vDesired{v},pDesired{p}];
        Yval = []; Yerr = [];
        for gname = grpu
            gname = gname{:};
            gIdx = strcmp(grp, gname);
            gTbl = erpTbl(gIdx,:);
            if p == P
                yval = eval(['gTbl.',vname]); 
                yerr = 0*yval;
            else
                yval = eval(['gTbl.',vname,'Mean']); 
                yerr = eval(['gTbl.',vname,'SD']);   
            end
            Yval = [Yval, yval]; Yerr = [Yerr, yerr];
        end
        plotBarWithErr(Yval, Yerr, 1);
        % label plot 
        X = erpTbl(:,2); 
        xlabel(X.Properties.VariableNames{1}); 
        X = X{:,1}; 
        xticklabels(X);
        ylabel([vNewLbl{v},' ',pNewLbl{p}]); title(rNames{r});
        lgd = [grpu, '±1 S.D.'];
        legend(lgd, 'Location','eastoutside');
    end
    end
end

%% helpers 
function plotBarWithErr(Yval, Yerr, lw)
    % bar plot values
    b = bar(Yval, 'EdgeColor','k', 'LineWidth',lw); grid on;
    if sum(abs(Yerr(:)))
        hold on;
        for c = 1:length(b)
            x = b(c).XEndPoints; 
            yv = Yval(:,c);
            ye = Yerr(:,c); 
            errorbar(x,yv,ye,ye, '.k', 'LineWidth',lw); 
        end
    end
end