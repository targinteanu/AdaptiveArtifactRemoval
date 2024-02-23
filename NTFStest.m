startdir = cd; 
searchdir = 'C:\Users\targi\Documents\Rodent_SSEP_Data\Data_AnesthesiaControl\saved_signal_obj\filtered_signal_obj';
cd(searchdir); 
[fn,fp] = uigetfile('*_filtered_sigObj.mat');
cd(startdir);
sig = loadSignalObj([fp,'\',fn])
Fs = sig.SampleRate;
sig = extractChannel(sig, 2);
%%
N = size(sig.Data_Unfiltered, 1) - size(sig.Data_LMS, 1) + 1;
tBeforeTrig = .04;
[t_PrePost, d_PrePost, e_PrePost] = getPrePostStim(tBeforeTrig, ...
    sig.Noise_Reference, sig.Data_HPF, sig.Data_LMS_LPF, Fs, sig.Channels, N);
%PrePostAvgAll_v2(tBeforeTrig,t_PrePost,d_PrePost,e_PrePost,Fs,sig.Channels,10);

%% NTFS 
% first, try getting it from averaged spectrogram 
[SUB, SUA, SFB, SFA] = PrePostAvgAll(tBeforeTrig,t_PrePost,d_PrePost,e_PrePost,Fs,sig.Channels,10);

%%
figure;
SS = SUA{4};
St = SUA{3};
Sf = SUA{2};
trl = randi(size(SS,3));
Sch = SS(:,:,trl);
Sti = (St>=.03); St = St(Sti);
Sfi = (Sf<=1500)&(Sf>=150); Sf = Sf(Sfi);
Sch = Sch(Sfi, Sti); 
imagesc(St, Sf, log(Sch)); colorbar;
xlabel('t (s)'); ylabel('f (Hz)'); title(['SUA trial ',num2str(trl)]);

%%
getAllTFS = @(S, fRng, tRng) S{1,4}( ...
                                     (S{1,2}>=min(fRng)-10) & (S{1,2}<=max(fRng))+10, ...
                                     (S{1,3}>=min(tRng)-10) & (S{1,3}<=max(tRng))+10);
NAllTFS = @(SM) squeeze(sum(sum(SM,1),2)/(size(SM,1)*size(SM,2)));
getAllNTFS = @(S, fRng, tRng) NAllTFS(getAllTFS(S, fRng, tRng));
getAvgNTFS = @(S, fRng, tRng) mean(getAllTFS(S, fRng, tRng), 'all', 'omitnan');
figure('Units','normalized', 'Position',[.1,.1,.8,.8])
Snames = {'SUB', 'SUA', 'SFB', 'SFA'};
for s = 1:length(Snames)
    Sname = Snames{s};
    Ss = eval(Sname);
    for ch = 1:size(Ss,1)
        subplot(size(Ss,1), length(Snames), length(Snames)*(ch-1) + s);
        Sch = Ss{ch,4}; St = Ss{ch,3}; Sf = Ss{ch,2};
        
        Sch = mean(Sch,3);
        Sti = (St>=.01); St = St(Sti);
        Sfi = (Sf<=950)&(Sf>=150); Sf = Sf(Sfi);
        Sch = Sch(Sfi, Sti); 
        imagesc(St, Sf, 10*log10(Sch)); colorbar;
        xlabel('t (s)'); ylabel('f (Hz)'); title([Sname,' ch',num2str(ch)]);
    end
end

figure('Units','normalized', 'Position',[.05,.1,.9,.8])
for ch = 1:size(SFA,1)
    NTFS_ch = cell(4,2); % cols = [200_400, 400_800] | rows = [UB; UA; FB; FA]

    NTFS_ch{1,1} = getAllNTFS(SUB(ch,:), [200,400], [.015,.025]);
    NTFS_ch{2,1} = getAllNTFS(SFB(ch,:), [200,400], [.015,.025]);
    NTFS_ch{3,1} = getAllNTFS(SUA(ch,:), [200,400], [.015,.025]);
    NTFS_ch{4,1} = getAllNTFS(SFA(ch,:), [200,400], [.015,.025]);

    NTFS_ch{1,2} = getAllNTFS(SUB(ch,:), [400,800], [.015,.025]);
    NTFS_ch{2,2} = getAllNTFS(SFB(ch,:), [400,800], [.015,.025]);
    NTFS_ch{3,2} = getAllNTFS(SUA(ch,:), [400,800], [.015,.025]);
    NTFS_ch{4,2} = getAllNTFS(SFA(ch,:), [400,800], [.015,.025]);

    subplot(size(SFA,1),2,2*(ch-1)+1);
    boxplot(cell2mat(NTFS_ch(:,1)'));
    xticklabels({'UB','FB','UA','FA'});
    ylabel(num2str(ch));
    title('200-400 Hz');

    subplot(size(SFA,1),2,2*ch);
    boxplot(cell2mat(NTFS_ch(:,2))+100);
    xticklabels({'UB','FB','UA','FA'});
    ylabel(num2str(ch));
    title('400-800 Hz');
end