%   Extracts raw data and removes artifact from neurophysiology data
%
%   Written by Alex Fanning, 6/21/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all

% Specify recording and artifact channels
prompt = {'Ephys ch#: ','Force plate ch #: ','EMG ch#: ','EMG ch#:','Artifact ch#: ','Artifact ch#: '};
dlgtitle = 'Recording channel details';
default = {'1','144','5','6','2','3'};
tempParams = inputdlg(prompt,dlgtitle,1,default);
chSig = str2double(tempParams);

% Extract data
params = {}; data = {}; count = 0; sig = cell(1,length(default));
params.chRef = 0;
for i = 1:3
    [data,params] = getData(params,data,i);
    for ii = 1:length(params(i).chIDs)
        count = count + 1;
        chIdx = find(chSig(count) == params(i).chIDs);
        sig{count} = data{i}(:,chIdx);
    end
end

% Choose correct EMG channel
figure()
for i = 1:length(params(2).chIDs)-1
    subplot(length(params(2).chIDs)-1,1,i)
    plot(sig{i+2})
    ylabel('\muVs')
    title(sprintf('%s %d','Ch',i))
    ylim([(prctile(sig{i+2},5)*10) (prctile(sig{i+2},95)*10)])
    if i == 2
        xlabel('Timepts')
    end
    set(gca,'FontSize',16)
end

prompt = {'EMG ch: '};
dlgtitle = 'EMG channel to use';
default = {'1'};
temp = inputdlg(prompt,dlgtitle,1,default);
chSig(end+1) = str2double(temp)+2;

% Upsample data
for i = 1:length(sig)
    if i == 1
        sig{2,i} = sig{1,1};
    else
        sig{2,i} = upsamp(sig{1,i},sig{1,1})';
    end
end

%% Filter out rhythmic artifact from EMG signal

tempSig = sig{1,chSig(end)};
signalAnalyzer
waitfor(gcf)
sig{3,3} = tempSig;
sig{3,4} = upsamp(sig{3,3},sig{1,1});

clear i ii temp tempParams prompt dlgtitle a1 a2 a3 a4 ax ax2 default a b order chIdx count tempSig

%% Identify best artifact channel

figure('Name','Artifact channels')

    ax = subplot(2,1,1); hold on
    plot(sig{1,5})
    title('Artifact channel 1')
    ylabel('\muVs')
    ylim([min(sig{1,5}) max(sig{1,5})])
    set(gca,'FontSize',16)
    
    ax2 = subplot(2,1,2); hold on
    plot(sig{1,6})
    ylabel('\muVs')
    title('Artifact channel 2')
    ylim([min(sig{1,6}) max(sig{1,6})])
    xlabel('Timepoints')
    set(gca,'FontSize',16)
    
    linkaxes([ax ax2],'x')

prompt = {'Artifact ch# to use: '};
dlgtitle = 'Artifact channel';
default = {'1'};
chSig(end+1) = str2double(inputdlg(prompt,dlgtitle,1,default))+4;

figure('Name','Raw data')

a1 = subplot(4,1,1); hold on
plot(sig{2,1},'k')
title('Ephys')
ylabel('\muVs')
ylim([-2000 2000])
set(gca,'FontSize',16)

a2 = subplot(4,1,2); hold on
plot(sig{2,2},'r')
title('Force plate')
ylabel('\muVs')
ylim([-500 1000])
set(gca,'FontSize',16)

a3 = subplot(4,1,3); hold on
plot(sig{3,4},'r')
title('EMG')
ylabel('\muVs')
ylim([-2000 2000])
set(gca,'FontSize',16)

a4 = subplot(4,1,4); hold on
plot(sig{2,chSig(end)})
title('Artifact')
xlabel('Timepoints')
ylabel('\muVs')
ylim([-2000 2000])
set(gca,'FontSize',16)
linkaxes([a1 a2 a3 a4],'x')

%% Remove artifact from neurophysiology channel

[newSig,segEndPts]  = removeArtifact(sig{1,1},sig{2,chSig(end)},params);

%% Process EMG signal

sig5rect = abs(sig{3,3});
sig5smooth = smooth(sig5rect,15);

data{5} = smooth(sig5smooth,100);

sig{4,4} = upsamp(sig5smooth,sig{1})';

clear a b order sig5rect artCh dlgtitle default a1 a2 a3 a4 ax ax2 prompt

save([params(1).recName])
