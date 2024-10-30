FontSizeAxisTitle = 16; 
FontSizeTickLbl = 12;
FontSizeTitle = 16;

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

%% number of N7s and SNR 

% summary sheets: 
    % (1,:) = unfilt; (2,:) = filt 
    % (:,1) = ERP; (:,2) = PSA

% number 
num_n07_unfilt = sum(summarySheetsAC{1,1}.n07num) + sum(summarySheetsCA{1,1}.n07num)
num_n15_unfilt = sum(summarySheetsAC{1,1}.n15num) + sum(summarySheetsCA{1,1}.n15num)
num_n07_filt = sum(summarySheetsAC{2,1}.n07num) + sum(summarySheetsCA{2,1}.n07num)
num_n15_filt = sum(summarySheetsAC{2,1}.n15num) + sum(summarySheetsCA{2,1}.n15num)

% mean lat 
lat_n07_unfilt = sum([summarySheetsAC{1,1}.n07latMean.*summarySheetsAC{1,1}.n07num; ...
                  summarySheetsCA{1,1}.n07latMean.*summarySheetsCA{1,1}.n07num]) ...
                 / num_n07_unfilt
lat_n15_unfilt = sum([summarySheetsAC{1,1}.n15latMean.*summarySheetsAC{1,1}.n15num; ...
                  summarySheetsCA{1,1}.n15latMean.*summarySheetsCA{1,1}.n15num]) ...
                 / num_n15_unfilt
lat_n07_filt = sum([summarySheetsAC{2,1}.n07latMean.*summarySheetsAC{2,1}.n07num; ...
                summarySheetsCA{2,1}.n07latMean.*summarySheetsCA{2,1}.n07num]) ...
               / num_n07_filt
lat_n15_filt = sum([summarySheetsAC{2,1}.n15latMean.*summarySheetsAC{2,1}.n15num; ...
                summarySheetsCA{2,1}.n15latMean.*summarySheetsCA{2,1}.n15num]) ...
               / num_n15_filt
% SD 
latSD_n07_unfilt = sum([summarySheetsAC{1,1}.n07latSD.*summarySheetsAC{1,1}.n07num; ...
                    summarySheetsCA{1,1}.n07latSD.*summarySheetsCA{1,1}.n07num].^2) ...
                   / num_n07_unfilt;
latSD_n15_unfilt = sum([summarySheetsAC{1,1}.n15latSD.*summarySheetsAC{1,1}.n15num; ...
                    summarySheetsCA{1,1}.n15latSD.*summarySheetsCA{1,1}.n15num].^2) ...
                   / num_n15_unfilt;
latSD_n07_filt = sum([summarySheetsAC{2,1}.n07latSD.*summarySheetsAC{2,1}.n07num; ...
                  summarySheetsCA{2,1}.n07latSD.*summarySheetsCA{2,1}.n07num].^2) ...
                 / num_n07_filt;
latSD_n15_filt = sum([summarySheetsAC{2,1}.n15latSD.*summarySheetsAC{2,1}.n15num; ...
                  summarySheetsCA{2,1}.n15latSD.*summarySheetsCA{2,1}.n15num].^2) ...
                 / num_n15_filt;
for v = {'latSD_n07_unfilt', 'latSD_n15_unfilt', 'latSD_n07_filt', 'latSD_n15_filt'}
    vname = v{:};
    eval([vname,' = sqrt(',vname,');']);
end
% SEM and CI 
for v = {'_n07_unfilt', '_n15_unfilt', '_n07_filt', '_n15_filt'}
    vname = v{:};
    eval(['latSEM',vname,' = latSD',vname,'/sqrt(num',vname,');']);
    eval(['latCI',vname,' = latSEM',vname,'*tinv([.025,.975],num',vname,')']);
    %eval(['latCI',vname,' = latCI',vname,' + lat',vname])
end

% SNR 
snr_unfilt = sum([summarySheetsAC{1,1}.snrMean.*summarySheetsAC{1,1}.n15num; ...
                  summarySheetsCA{1,1}.snrMean.*summarySheetsCA{1,1}.n15num]) ...
                 / num_n15_unfilt
snr_filt = sum([summarySheetsAC{2,1}.snrMean.*summarySheetsAC{2,1}.n15num; ...
                summarySheetsCA{2,1}.snrMean.*summarySheetsCA{2,1}.n15num]) ...
               / num_n15_filt
snrSD_unfilt = sum([summarySheetsAC{1,1}.snrSD.*summarySheetsAC{1,1}.n15num; ...
                    summarySheetsCA{1,1}.snrSD.*summarySheetsCA{1,1}.n15num].^2) ...
                   / num_n15_unfilt;
snrSD_filt = sum([summarySheetsAC{2,1}.snrSD.*summarySheetsAC{2,1}.n15num; ...
                  summarySheetsCA{2,1}.snrSD.*summarySheetsCA{2,1}.n15num].^2) ...
                 / num_n15_filt;
snrSD_unfilt = sqrt(snrSD_unfilt); snrSD_filt = sqrt(snrSD_filt);
snrSEM_unfilt = snrSD_unfilt/sqrt(num_n15_unfilt);
snrSEM_filt = snrSD_filt/sqrt(num_n15_filt);
snrCI_unfilt = snrSEM_unfilt*tinv([.025,.975],num_n15_unfilt)
snrCI_filt = snrSEM_filt*tinv([.025,.975],num_n15_filt)
S = (snrSD_unfilt^2)/num_n15_unfilt + (snrSD_filt^2)/num_n15_filt; % combined var
DoFd = ((snrSD_unfilt/num_n15_unfilt)^2)/(num_n15_unfilt-1) + ...
       ((snrSD_filt/num_n15_filt)^2)/(num_n15_filt-1) ; % denom
DoF = (S^2)/DoFd; % deg of freedom
T = (snrSD_filt-snrSD_unfilt)/sqrt(S);
p = tcdf(T,DoF,'upper')

%% AC PSA Plots 
figure('Units','normalized', 'Position',[.1,.1,.8,.4]);

psaTbl = summarySheetsAC{1,2}; % unfilt AC PSA
ax(1) = subplot(1,2,1);
Yval = []; Ystd = []; Ynum = [];
for vname = {'PSAsub', 'PSA'}
    vname = vname{:};
        yval = eval(['psaTbl.',vname,'mean']); 
        ystd = eval(['psaTbl.',vname,'std']); 
        ynum = eval(['psaTbl.',vname,'num']);
        ysem = ystd./sqrt(ynum);
        %yerr = ysem.*tinv(.975,ynum);
    Yval = [Yval, yval]; Ystd = [Ystd, ystd]; Ynum = [Ynum, ynum];
end
plotBarWithErr(Yval, Ystd, Ynum, 1.5, 5e-7);
% label plot 
ax(1).FontSize = FontSizeTickLbl;
X = psaTbl(:,1); 
xlabel(X.Properties.VariableNames{1}, 'FontSize',FontSizeAxisTitle); 
X = X{:,1}; X = arrayfun(@(x) [num2str(x),'%'], X, 'UniformOutput',false);
xticklabels(X);
ylabel('PSA', 'FontSize',FontSizeAxisTitle); 
title('PSA, Post-Artifact', 'FontSize',FontSizeTitle);
lgd = {'5-9ms', '5-30ms', '±1 S.E.'};
%lgd = {'5-9ms', '5-30ms', '95% C.I.'};
legend(lgd, 'Location','eastoutside', 'FontSize',FontSizeTickLbl);

psaTbl = summarySheetsAC{2,2}; % filt AC PSA
ax(2) = subplot(1,2,2);
Yval = []; Ystd = []; Ynum = [];
for vname = {'PSA09', 'PSA30'}
    vname = vname{:};
        yval = eval(['psaTbl.',vname,'mean']); 
        ystd = eval(['psaTbl.',vname,'std']); 
        ynum = eval(['psaTbl.',vname,'num']);
        ysem = ystd./sqrt(ynum);
        %yerr = ysem.*tinv(.975,ynum);
    Yval = [Yval, yval]; Ystd = [Ystd, ystd]; Ynum = [Ynum, ynum];
end
plotBarWithErr(Yval, Ystd, Ynum, 1.5, 5e-7);
% label plot 
ax(2).FontSize = FontSizeTickLbl;
X = psaTbl(:,1); 
xlabel(X.Properties.VariableNames{1}, 'FontSize',FontSizeAxisTitle); 
X = X{:,1}; X = arrayfun(@(x) [num2str(x),'%'], X, 'UniformOutput',false);
xticklabels(X);
ylabel('PSA', 'FontSize',FontSizeAxisTitle); 
title('PSA, Artifact-Removed', 'FontSize',FontSizeTitle);
lgd = {'0-9ms', '0-30ms', '±1 S.E.'};
%lgd = {'0-9ms', '0-30ms', '95% C.I.'};
legend(lgd, 'Location','eastoutside', 'FontSize',FontSizeTickLbl);

linkaxes(ax, 'y'); clear ax ax_
print('AC_PSA.png','-dpng')

%% AC ERP Plots 
figure('Units','normalized', 'Position',[.1,.1,.8,.4]);
vDesired = {'n07', 'n15'}; vNewLbl = {'N7', 'N10'};

erpTbl = summarySheetsAC{1,1}; % unfilt AC ERP
ax(1) = subplot(1,2,1);
Yval = []; Ystd = []; Ynum = [];
for v = vDesired
    vname = [v{:},'amp'];
        yval = eval(['erpTbl.',vname,'Mean']); 
        ystd = eval(['erpTbl.',vname,'SD']);   
        ynum = eval(['erpTbl.',v{:},'num']);
        ysem = ystd./sqrt(ynum);
        %yerr = ysem.*tinv(.975,ynum);
    Yval = [Yval, yval]; Ystd = [Ystd, ystd]; Ynum = [Ynum, ynum];
end
plotBarWithErr(Yval, Ystd, Ynum, 1.5, 3);
% label plot 
ax(1).FontSize = FontSizeTickLbl;
X = erpTbl(:,1); 
xlabel(X.Properties.VariableNames{1}, 'FontSize',FontSizeAxisTitle); 
X = X{:,1}; X = arrayfun(@(x) [num2str(x),'%'], X, 'UniformOutput',false);
xticklabels(X);
ylabel('Amplitude (\muV)', 'FontSize',FontSizeAxisTitle); 
title('Unfiltered', 'FontSize',FontSizeTitle);
lgd = [vNewLbl, '±1 S.E.'];
%lgd = [vNewLbl, '95% C.I.'];
legend(lgd, 'Location','eastoutside', 'FontSize',FontSizeTickLbl);

erpTbl = summarySheetsAC{2,1}; % filt AC ERP
ax(2) = subplot(1,2,2);
Yval = []; Ystd = []; Ynum = [];
for v = vDesired
    vname = [v{:},'amp'];
        yval = eval(['erpTbl.',vname,'Mean']); 
        ystd = eval(['erpTbl.',vname,'SD']);   
        ynum = eval(['erpTbl.',v{:},'num']);
        ysem = ystd./sqrt(ynum);
        %yerr = ysem.*tinv(.975,ynum);
    Yval = [Yval, yval]; Ystd = [Ystd, ystd]; Ynum = [Ynum, ynum];
end
plotBarWithErr(Yval, Ystd, Ynum, 1.5, 3);
% label plot 
ax(2).FontSize = FontSizeTickLbl;
X = erpTbl(:,1); 
xlabel(X.Properties.VariableNames{1}, 'FontSize',FontSizeAxisTitle); 
X = X{:,1}; X = arrayfun(@(x) [num2str(x),'%'], X, 'UniformOutput',false);
xticklabels(X);
ylabel('Amplitude (\muV)', 'FontSize',FontSizeAxisTitle); 
title('Artifact-Removed', 'FontSize',FontSizeTitle);
lgd = [vNewLbl, '±1 S.E.'];
%lgd = [vNewLbl, '95% C.I.'];
legend(lgd, 'Location','eastoutside', 'FontSize',FontSizeTickLbl);

linkaxes(ax, 'y'); clear ax
print('AC_ERP.png','-dpng')

%% CA PSA Plots 
figure('Units','normalized', 'Position',[.05,.1,.9,.8]);

psaTbl = summarySheetsCA{1,2}; % unfilt CA PSA
grp = psaTbl{:,1}; grpu = unique(grp)';
for g = 1:length(grpu)
    gname = grpu{g};
    gIdx = strcmp(grp, gname);
    gTbl = psaTbl(gIdx,:);
    ax(g) = subplot(2,2, g);
    Yval = []; Ystd = []; Ynum = [];
    for vname = {'PSAsub', 'PSA'}
        vname = vname{:};
            yval = eval(['gTbl.',vname,'mean']); 
            ystd = eval(['gTbl.',vname,'std']);   
            ynum = eval(['gTbl.',vname,'num']);
        Yval = [Yval, yval]; Ystd = [Ystd, ystd]; Ynum = [Ynum, ynum];
    end
    plotBarWithErr(Yval, Ystd, Ynum, 1.5, 1e-6);
    % label plot 
    ax(g).FontSize = FontSizeTickLbl;
    X = gTbl(:,2); 
    xlabel(X.Properties.VariableNames{1}, 'FontSize',FontSizeAxisTitle); 
    X = X{:,1}; 
    xticklabels(X);
    title(['PSA - ',gname,' Injury'], 'FontSize',FontSizeTitle); 
    ylabel('PSA, Post-Artifact', 'FontSize',FontSizeAxisTitle);
    lgd = {'5-9ms', '5-30ms', '±1 S.E.'};
    legend(lgd, 'Location','eastoutside', 'FontSize',FontSizeTickLbl);
end

linkaxes(ax,'y'); clear ax

psaTbl = summarySheetsCA{2,2}; % filt CA PSA
grp = psaTbl{:,1}; grpu = unique(grp)';
for g = 1:length(grpu)
    gname = grpu{g};
    gIdx = strcmp(grp, gname);
    gTbl = psaTbl(gIdx,:);
    ax(g) = subplot(2,2, g+2);
    Yval = []; Ystd = []; Ynum = [];
    for vname = {'PSA09', 'PSA30'}
        vname = vname{:};
            yval = eval(['gTbl.',vname,'mean']); 
            ystd = eval(['gTbl.',vname,'std']);   
            ynum = eval(['gTbl.',vname,'num']);
        Yval = [Yval, yval]; Ystd = [Ystd, ystd]; Ynum = [Ynum, ynum];
    end
    plotBarWithErr(Yval, Ystd, Ynum, 1.5, 2e-6);
    % label plot 
    ax(g).FontSize = FontSizeTickLbl;
    X = gTbl(:,2); 
    xlabel(X.Properties.VariableNames{1}, 'FontSize',FontSizeAxisTitle); 
    X = X{:,1}; 
    xticklabels(X);
    title(['PSA - ',gname,' Injury'], 'FontSize',FontSizeTitle); 
    ylabel('PSA, Artifact-Removed', 'FontSize',FontSizeAxisTitle);
    lgd = {'0-9ms', '0-30ms', '±1 S.E.'};
    legend(lgd, 'Location','eastoutside', 'FontSize',FontSizeTickLbl);
end

linkaxes(ax,'y'); clear ax
print('CA_PSA.png','-dpng')

%% CA ERP Plots 
figure('Units','normalized', 'Position',[.05,.1,.9,.8]);
vDesired = {'n07', 'n15'}; vNewLbl = {'N7', 'N10'};

erpTbl = summarySheetsCA{1,1}; % unfilt CA ERP
grp = erpTbl{:,1}; grpu = unique(grp)';
for g = 1:length(grpu)
    gname = grpu{g};
    gIdx = strcmp(grp, gname);
    gTbl = erpTbl(gIdx,:);
    ax(g) = subplot(2,2, g);
    Yval = []; Ystd = []; Ynum = [];
    for v = vDesired
        vname = [v{:},'amp'];
            yval = eval(['gTbl.',vname,'Mean']); 
            ystd = eval(['gTbl.',vname,'SD']);   
            ynum = eval(['gTbl.',v{:},'num']);
        Yval = [Yval, yval]; Ystd = [Ystd, ystd]; Ynum = [Ynum, ynum];
    end
    plotBarWithErr(Yval, Ystd, Ynum, 1.5, 10);
    % label plot 
    ax(g).FontSize = FontSizeTickLbl;
    X = gTbl(:,2); 
    xlabel(X.Properties.VariableNames{1}, 'FontSize',FontSizeAxisTitle); 
    X = X{:,1}; 
    xticklabels(X);
    title(['Amplitude (\muV) - ',gname,' Injury'], 'FontSize',FontSizeTitle); 
    ylabel('Amplitude (\muV), Unfiltered', 'FontSize',FontSizeAxisTitle);
    lgd = [vNewLbl, '±1 S.E.'];
    legend(lgd, 'Location','eastoutside', 'FontSize',FontSizeTickLbl);
end

linkaxes(ax,'y'); clear ax

erpTbl = summarySheetsCA{2,1}; % filt CA ERP
grp = erpTbl{:,1}; grpu = unique(grp)';
for g = 1:length(grpu)
    gname = grpu{g};
    gIdx = strcmp(grp, gname);
    gTbl = erpTbl(gIdx,:);
    ax(g) = subplot(2,2, g+2);
    Yval = []; Ystd = []; Ynum = [];
    for v = vDesired
        vname = [v{:},'amp'];
            yval = eval(['gTbl.',vname,'Mean']); 
            ystd = eval(['gTbl.',vname,'SD']);   
            ynum = eval(['gTbl.',v{:},'num']);
        Yval = [Yval, yval]; Ystd = [Ystd, ystd]; Ynum = [Ynum, ynum];
    end
    plotBarWithErr(Yval, Ystd, Ynum, 1.5, 10);
    % label plot 
    ax(g).FontSize = FontSizeTickLbl;
    X = gTbl(:,2); 
    xlabel(X.Properties.VariableNames{1}, 'FontSize',FontSizeAxisTitle); 
    X = X{:,1}; 
    xticklabels(X);
    title(['Amplitude (\muV) - ',gname,' Injury'], 'FontSize',FontSizeTitle); 
    ylabel('Amplitude (\muV), Artifact-Removed', 'FontSize',FontSizeAxisTitle);
    lgd = [vNewLbl, '±1 S.E.'];
    legend(lgd, 'Location','eastoutside', 'FontSize',FontSizeTickLbl);
end

linkaxes(ax,'y'); clear ax
print('CA_ERP.png','-dpng')

%% helpers 
function plotBarWithErr(Yval, Ystd, Ynum, lw, ysigspc)
    Yerr = Ystd./sqrt(Ynum); 
    pAlpha = .05; 
    %ysigspc = .1;

    ymax = Yval+Yerr; ymax = max(ymax(:)); 
    ymin = Yval-Yerr; ymin = min(ymin(:)); ymin = min(0,ymin);
    %ysigspc = ysigspc*(ymax-ymin);

    % bar plot values
    b = bar(Yval, 'EdgeColor','k', 'LineWidth',lw); grid on;
    if sum(abs(Ystd(:)))
        hold on;
        for c = 1:length(b)
            x = b(c).XEndPoints; 

            % error bars 
            yv = Yval(:,c);
            ye = Yerr(:,c); 
            errorbar(x,yv,ye,ye, '.k', 'LineWidth',lw); 

            % calc significance  
            P = ones(size(Yval,1));
            for r1 = 1:size(Yval,1)
                s1 = Ystd(r1,c)^2; n1 = Ynum(r1,c); m1 = Yval(r1,c);
                for r2 = (r1+1):size(Yval,1)
                    s2 = Ystd(r2,c)^2; n2 = Ynum(r2,c); m2 = Yval(r2,c);
                    S = s1/n1 + s2/n2; % combined var
                    DoFd = ((s1/n1)^2)/(n1-1) + ((s2/n2)^2)/(n2-1); % denom
                    DoF = (S^2)/DoFd; % deg of freedom
                    T = (m1 - m2)/sqrt(S);
                    if T > 0
                        P(r1,r2) = tcdf(T,DoF,'upper');
                    else
                        P(r1,r2) = tcdf(T,DoF);
                    end
                end
            end

            % significance bars 
            r1=1; r2=2;
            Q = true(size(P)); % inds of P left to plot 
            incrementY = true;
            for d = 1:size(Yval,1)
                Q(d,1:d) = false;
            end
            while sum(Q(:))
                if (r1 > size(Yval,1)) | (r2 > size(Yval,1)) | ~Q(r1,r2)
                    [r1,r2] = find(Q); 
                    r1 = r1(1); r2 = r2(1); 
                    incrementY = true;
                end
                p = P(r1,r2); Q(r1,r2) = false;
                if p <= pAlpha
                    x1 = x(r1); x2 = x(r2); 
                    x1 = x1 + .05*b(c).BarWidth; x2 = x2 - .05*b(c).BarWidth;
                    xc = mean([x1,x2]); xb = abs(x2-xc);
                    if incrementY
                        ymax = ymax + ysigspc;
                        incrementY = false;
                    end
                    errorbar(xc,ymax,0,0,xb,xb, 'k','LineWidth',lw);
                    text(xc,ymax, ...
                        ['p = ',num2str(p,1)], ...
                        'FontSize',12, 'FontWeight','bold', ...
                        'HorizontalAlignment','center', ...
                        'VerticalAlignment','bottom');
                    r1 = r2; r2 = r1+1;
                else
                    r2 = r2+1;
                end
                if r2 > size(Yval,1)
                    r1 = r1+1; r2 = r1+1;
                    incrementY = true;
                end
            end
        end
    end
end