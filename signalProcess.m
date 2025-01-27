%   Signal processing for ephys
%
%   Written by Alex Fanning, 1/27/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all;

% Load relevant workspace
load(uigetfile("*.mat"))

% Import spike sorted data
nexFileData = readNexFile([params(1).recName '.nex']);

%% Grab waveform and ephys data

% Find waveform channels
params(1).wfCh = cell(1);
params(1).wfCh{1} = 'chan1a_wf';
params(1).wfCh{2} = 'chan1b_wf';
for tt = 1:length(params(1).wfCh)
    for i = 1:length(nexFileData.waves)
        if matches(nexFileData.waves{i}.name,params(1).wfCh{tt})
            params(1).wfCh{tt+2} = i;
        end
    end
end

% Grab waveform timepoints and waveforms
waveformData = struct();
waveformData(1).waveformTimepts = cell(1); waveformData(1).waveformTimeptsMS = cell(1); waveformData(1).waveforms = cell(1);
for i = 1:length(params(1).wfCh) - 2
    waveformData(1).waveformTimepts{i} = fix(nexFileData.waves{params(1).wfCh{i+2},1}.timestamps * params(1).sf); % Timepts in seconds
    waveformData(1).waveformTimeptsMS{i} = nexFileData.waves{params(1).wfCh{i+2},1}.timestamps * 1000; % Timepts in milliseconds
    waveformData(1).waveforms{i} = nexFileData.waves{params(1).wfCh{i+2},1}.waveforms;
end

iter = length(waveformData(1).waveformTimepts);

% Find ephys channel
for i = 1:length(nexFileData.contvars)
    if contains('chan1',nexFileData.contvars{i}.name)
        params(1).ephysChNum = i;
    end
end

sig{3,1} = nexFileData.contvars{params(1).ephysChNum,1}.data * 5000;

%% Remove artifact and discontinuous data

% Remove artifact from filtered signal
for i = 1:length(sig{2,1})
    if isnan(newSig(i))
        sig{2,1}(i) = NaN;
        sig{3,1}(i) = NaN;
        sig{4,4}(i) = NaN;
    end
end

for i = 1:length(segEndPts)
    sig{2,1}(segEndPts(i,1):segEndPts(i,2)) = NaN;
    sig{3,1}(segEndPts(i,1):segEndPts(i,2)) = NaN;
    sig{4,4}(segEndPts(i,1):segEndPts(i,2)) = NaN;
end

% Remove discontinuous data
params(1).disconThresh = params(1).sf/2;
nanTimepts = find(isnan(sig{2,1}));
chunkLngth = nanTimepts(2:end) - nanTimepts(1:end-1);
a = 1;
for i = 1:length(chunkLngth)
    if chunkLngth(i) ~= 1
        nanEnds(a) = i;
        a = a + 1;
    end
end

nanEnds = nanTimepts(nanEnds);
chunkLngth(chunkLngth==1) = [];
goodChnkStart = nanEnds + 1;
goodChnkStart(chunkLngth >= params(1).disconThresh) = [];
goodChnkLngths = chunkLngth;
goodChnkLngths(goodChnkLngths >= params(1).disconThresh) = [];

for i = 1:length(goodChnkStart)
    sig{2,1}(goodChnkStart(i):goodChnkStart(i) + goodChnkLngths) = NaN;
    sig{3,1}(goodChnkStart(i):goodChnkStart(i) + goodChnkLngths) = NaN;
    sig{4,4}(goodChnkStart(i):goodChnkStart(i) + goodChnkLngths) = NaN;
end

% Remove waveforms that have been removed from the raw ephys trace
for ii = 1:iter
    for i = 1:length(waveformData(1).waveformTimepts{ii})
        if isnan(sig{2,1}(waveformData(1).waveformTimepts{ii}(i)))
            waveformData(1).waveforms{ii}(:,i) = NaN;
        end
    end
end

%% Computer firing frequency after calculating total sampling time

% Calculate total sampling time after artifact has been removed
params(1).totalSampTime = length(sig{1,1}) - sum(isnan(sig{2,1}));

% Compute firing frequency
data2export = struct();
for i = 1:iter
    tempWfs = waveformData(1).waveforms{i}(~isnan(waveformData(1).waveforms{i}));
    tempWfs2 = reshape(tempWfs,size(waveformData(1).waveforms{i},1),length(tempWfs)/ size(waveformData(1).waveforms{i},1));
    newWaveforms{i} = tempWfs2;
    data2export(1).firingFreq(i) = length(newWaveforms{i}) / (params(1).totalSampTime / params(1).sf);
end

clear a nanEnds chunkLngth goodChnkStart goodChnkLngths a tempWfs tempWfs2 totalSampTime nanTimepts newSig

%% Inter-spike interval calculations

for ii = 1:iter
    data2export(1).isi{ii} = NaN(1,length(waveformData(1).waveformTimeptsMS{ii}) - 1);
    for i = 2:length(waveformData(1).waveformTimeptsMS{ii}) 
        data2export.isi{ii}(i-1) = waveformData(1).waveformTimeptsMS{ii}(i) - waveformData(1).waveformTimeptsMS{ii}(i - 1);
    end

    data2export.isiPkIdx = [];
    figure('Name','Fitted curve ISI'); hold on
    histfit(data2export.isi{ii},1000,'kernel');
    pd = fitdist(data2export.isi{ii}','kernel');
    x = linspace(1,150,500);
    y = pdf(pd,x);
%     plot(x,y,'r','LineWidth',2)
    xlim([-5 150])
    xlabel('ISI (ms)')
    ylabel('Number of occurrences')
    set(gca,'FontSize',16)
    [~,data2export.isiPkIdx(ii,:)] = findpeaks(y,x,'MinPeakProminence',0.0005);
    txt = text(75,50,sprintf('%d %s', round(data2export.firingFreq(ii)),'Hz'),"FontSize",16);

    % CV = Irregularity of spiking across record
    data2export.cv(ii) = std(data2export.isi{ii}) / mean(data2export.isi{ii});

    % Cv2 = Irregularity of adjacent spike times
    for i = 2:length(data2export.isi{ii})
        tempCv2(i-1) = (std([data2export.isi{ii}(i-1) data2export.isi{ii}(i)]) / mean([data2export.isi{ii}(i-1) data2export.isi{ii}(i)])) * sqrt(2);
    end
    data2export.cv2(ii) = mean(tempCv2);
end

%% Simple spike bursts

consecTimepts = cell(1); consecTimeptsIdxs = cell(1);
ssBurstIdx = cell(1); burstSrndTimepts = cell(1); ssBurstPause = cell(1);
burstSrnd = cell(1);
for ii = 1:length(data2export.isi)
    endpts = [];
    start = 1;
    ending = 1;
    count = 1;
    
    if data2export.isiPkIdx(ii,1) < 5
        params(1).burstThr = 4;
    else
        params(1).burstThr = 20;
    end
    
    for i = 1:length(data2export.isi{ii})
     
        if data2export.isi{ii}(i) <= params(1).burstThr
            ending = i;
        else
            endpts = [endpts; start, ending];
            start = i;
            ending = i;
        end
    
    end
    consecTimepts{ii} = endpts(:,2) - endpts(:,1);
    
    a = 1;
    for i = 1:length(consecTimepts{ii})
        if consecTimepts{ii}(i) > 0
            ssBurstIdx{ii}(a,1) = endpts(i,1) + 1;
            ssBurstIdx{ii}(a,2) = endpts(i,2) + 1;
            a = a + 1;
        end
    end

    ssBurstCount(ii,1) = 0;
    a = 1;
    for i = 1:7
    
        if i <= 6
            tempConsec = find(consecTimepts{ii} == i);
            tempIdx = endpts(tempConsec,1) + 1;
            temp2idx = endpts(tempConsec,2) + 1;
        else
            tempConsec = find(consecTimepts{ii} >= i);
            tempIdx = endpts(tempConsec,1) + 1;
            temp2idx = endpts(tempConsec,2) + 1;
        end

        % consecTimeptsIdxs refers to index for waveform variables
        consecTimeptsIdxs{ii}{1,i} = tempIdx;
        consecTimeptsIdxs{ii}{2,i} = temp2idx;
        ssBurstCount(ii,1) = ssBurstCount(1) + numel(tempConsec);
        ssBurstCount(ii,1+i) = length(tempConsec);

    end

    clear tempConsec tempIdx temp2idx tempCv2 numBins x y pd txt

    data2export.burstFreq(ii) = sum(ssBurstCount(ii,:)) / (length(sig{2,1}) / params(1).sf);
    
    figure('Name','Burst size freq')
    subplot(1,ii,ii); hold on
    bar(ssBurstCount(ii,:),'FaceColor','k')
    hMax(ii) = max(ssBurstCount(ii,:));
    xticks(1:8)
    xticklabels({2:8 '9+'})
    xlabel('Number of spikes in burst')
    ylabel('Number of bursts')
    title('Burst size occurrence')
    txt = text(6,hMax(ii)-50,sprintf('%d %s', round(data2export.burstFreq(ii)),'Hz'),"FontSize",16);
    set(gca,'FontSize',16)

    % Determine pauses after bursts
    [burstSrndTimepts{ii},ssBurstPause{ii}] = peristimData(waveformData(1).waveformTimepts{ii},ssBurstIdx{ii},sig{2,1}, params,1);
    data2export.burstPauseAvg(ii) = mean(ssBurstPause{ii});

    burstSrnd{ii} = cat(2,burstSrndTimepts{ii}{:});
    burstSrnd{ii}(burstSrnd{ii}==0) = [];
    burstSrnd{ii} = (burstSrnd{ii} / params(1).sf) * 1000;

    figure('Name','Burst pause to next AP'); hold on
    histogram(ssBurstPause{ii},1000)
    xlim([0 100])
    xlabel('Time to next action potential (ms)')
    ylabel('Number of events')
    title('Time to next action potential after burst')
    set(gca,'FontSize',16)

    figure('Name','Burst pause'); hold on
    hp = histogram(burstSrnd{ii},250);
    maxHp = median(hp.Values) + 50;
    xlim([-100 100])
    ylim([0 maxHp])
    xlabel('Time (ms)')
    ylabel('Number of events')
    set(gca,'FontSize',16)
end

clear count a tt i ii hMax hp h ending endpts maxHist maxHp numBins start tempCv2 txt prompt

%% Extract complex spike events

if params(1).burstThr == 5
    ch_data = sig{1,1}';
    ch_time = 1:length(ch_data);
    ch_time = ch_time / params(1).sf;
    sample_rate = params(1).sf;

    save([params(1).recName 'cs.mat'],'ch_data','ch_time','sample_rate')

    clear ch_data ch_time sample_rate
end

%% AUTOCORRELATION
% test = autocorr(isi)

%% Compute PSD

[data2export.mPSD(:,1),data2export.maxPSD(1),data2export.maxPSDidx(1),data2export.nPSD(:,1),data2export.nPSDmax(1),data2export.nPSDmaxIdx(1),params] = computePSD(data,params,20,2,2,1,1);
[data2export.mPSD(:,2),data2export.maxPSD(2),data2export.maxPSDidx(2),data2export.nPSD(:,2),data2export.nPSDmax(2),data2export.nPSDmaxIdx(2),params] = computePSD(sig5smooth,params,,data2export,20,2,5,1,1);

%% Extract tremor events
data2export.tremEpochs = {}; epochIdxs = {}; moveEpochs = {};
[data2export.tremEpochs,moveEpochs,epochIdxs,data] = extractFPtremor(sig{1,2},data,params,data2export,1);
[data2export.tremEpochs,moveEpochs,epochIdxs,data] = extractFPtremor(data{5},data,params,data2export,2);

save([params(1).recName])

%% Further filter tremor events
type = 1; tremIdxs = [];
waitfor(highlightSegments(type,sig{2,2},tremIdxs))
type = 2;
waitfor(highlightSegments(type,sig{2,3},tremIdxs.goodTremIdxs))

%% 
tempFiltIdx{1} = find(tremIdxs.goodTremIdxs(1,:)==1);
tempFiltIdx{2} = find(tremIdxs.goodTremIdxs(2,:)==1);
for ii = 1:2
    a = 1;
    for i = tempFiltIdx{ii}
        filtTremIdxs{ii}(a,1) = epochIdxs{ii,1}(i,1);
        filtTremIdxs{ii}(a,2) = epochIdxs{ii,1}(i,2);
        filtEphysEpochs{ii,a} = sig{2,1}(filtTremIdxs{ii}(a,1):filtTremIdxs{ii}(a,2));
        filtTremPeaks{ii,a} = data2export.tremEpochs{ii,2}{i};
        filtTremLocs{ii,a} = data2export.tremEpochs{ii,3}{i};
        filtTremRel{ii,a} = data2export.tremEpochs{ii,8}(i);
        a = a + 1;
    end
end
 
% Find spike times for filtered tremor events
[filtTremSpikes{1},~] = peristimData(waveformData(1).waveformTimepts{1},filtTremIdxs{1},sig{2,1},params,4);
[filtTremSpikes{2},~] = peristimData(waveformData(1).waveformTimepts{1},filtTremIdxs{2},sig{2,1},params,4);
for ii = 1:2
    for i = 1:length(filtTremSpikes{ii})
        filtTremSpikeTimes{ii,i} = waveformData(1).waveformTimepts{1}(filtTremSpikes{ii}{i});
        tempIdx = find(filtTremIdxs{ii}(i,1)==epochIdxs{ii,1}(:,1));
        filtTremEpochs{ii,i} = data2export.tremEpochs{ii,1}{tempIdx};
        filtEphysTimes{ii,i} = filtTremSpikeTimes{ii,i} - filtTremIdxs{ii}(i,1);
        invFiltTremEpochs{ii,i} = filtTremEpochs{ii,i}*-1;
    end
end

% Find peaks of each tremor event
for ii = 1:2
    % Parameters for finding peaks
    prompt = {'minPeakProm: ','minPeakInt: ','maxPeakWidth: ','minPeakWidth: '};
    dlgtitle = 'FindPeaks params';
    if ii == 1
        default = {'40','750','3000','50'};
    elseif ii == 2
        default = {'5','250','3000','50'};
    end
    tempParams = inputdlg(prompt,dlgtitle,1,default);
    parameters(1).minPeakProm = str2double(tempParams{1});
    parameters(1).minPeakDist =  str2double(tempParams{2});
    parameters(1).maxPeakWidth =  str2double(tempParams{3});
    parameters(1).minPeakWidth =  str2double(tempParams{4});
    for i = 1:length(filtTremIdxs{ii})
        [pks,idxs,w,p] = findpeaks(invFiltTremEpochs{ii,i},1,'MinPeakProminence',parameters(1).minPeakProm,'MinPeakDistance',parameters(1).minPeakDist,'Annotate','extents','MaxPeakWidth',parameters(1).maxPeakWidth,'MinPeakWidth',parameters(1).minPeakWidth);
        invPeaks{ii,i} = pks; invLocs{ii,i} = idxs; invWidths{ii,i} = w; invOtherPeaks{ii,i} = p;
        % For trough detection: Flip the sign of each tremorEpoch and rerun
        % peak detection to get locs and peaks
    end
    clear pks idxs w p
end

%% 
figure()
counter = [1,1;1,1];
for tt = 1:2
    for i = 1:length(filtTremIdxs{tt})
        for ii = 1:length(filtEphysTimes{tt,i})
            filtTremEphysDiff{tt,i}(ii,:) = filtEphysTimes{tt,i}(ii) - filtTremLocs{tt,i}(:);
            filtTremEphysDiff2{tt,i}(ii,:) = filtEphysTimes{tt,i}(ii) - invLocs{tt,i}(:);
        end
        if isempty(filtTremEphysDiff{tt,i}) || isempty(filtTremEphysDiff2{tt,i})
            continue
        else
            for k = 1:size(filtTremEphysDiff{tt,i},1)
                [filtTremEphysTimeTemp{tt,i}(k),filtTremEphysTimeTempIdx{tt,i}(k)] = min(abs(filtTremEphysDiff{tt,i}(k,:)));
            end
        
            for jj = 1:size(filtTremEphysDiff2{tt,i},1)
                [filtTremEphysTimeTemp2{tt,i}(jj),filtTremEphysTimeTempIdx2{tt,i}(jj)] = min(abs(filtTremEphysDiff2{tt,i}(jj,:)));
            end
        
            for mm = 1:length(filtTremEphysTimeTemp{tt,i})
                if filtTremEphysTimeTemp{tt,i}(mm) < filtTremEphysTimeTemp2{tt,i}(mm)
                    filtTremEphysTime{tt,i}(mm) = filtTremEphysDiff{tt,i}(mm,filtTremEphysTimeTempIdx{tt,i}(mm));
                    counter(tt,1) = counter(tt,1) + 1;
                else
                    filtTremEphysTime{tt,i}(mm) = filtTremEphysDiff2{tt,i}(mm,filtTremEphysTimeTempIdx2{tt,i}(mm));
                    counter(tt,2) = counter(tt,2) + 1;
                end
            end
        end
    end
    filtTremEphysTiming{tt} = cat(2,filtTremEphysTime{tt,:});
    filtTremEphysTimingMS{tt} = (filtTremEphysTiming{tt} / params(1).sf) * 1000;

    
    subplot(1,2,tt); hold on
    histogram(filtTremEphysTimingMS{tt},50)
end

    
%% Produce peristimulus histograms centered on all tremor epochs

for tt = 1:2
    if exist('epochIdxs')
        [tremSrndSStimepts{tt},~] = peristimData(waveformData(1).waveformTimepts{1},epochIdxs{tt,1},sig{2,1},params,2);
        [movSrndSStimepts{tt},~] = peristimData(waveformData(1).waveformTimepts{1},epochIdxs{tt,2},sig{2,1},params,2);
    end
    
    if exist('data2export.tremEpochs')
        if ~isempty(data2export.tremEpochs{tt,8})
            for i = 1:length(tremSrndSStimepts{tt})
                if epochIdxs{tt,1}(i,1) - 2999 < 1 || epochIdxs{tt,1}(i,2) + 3000 > length(sig{2,1})
                    continue
                else
                    tempEphys = sig{2,1}(epochIdxs{tt,1}(i,1) - 2999:epochIdxs{tt,1}(i,2) + 3000);
                    tremTrigEphys{tt,i} = circshift(tempEphys,data2export.tremEpochs{tt,8})';
                    tremSrndSStimeptsTest{tt,i} = tremSrndSStimepts{tt}{i} - data2export.tremEpochs{tt,8}(i);
                end
            end
        end
    end
    
    % Create variable argutment in and out function
    % Shift data to align tremor peaks
    % Add tremor timepts to histogram
    
    if exist('epochIdxs')
        tremSrndSS{tt} = cat(2,tremSrndSStimeptsTest{tt,:});
        tremSrndSS{tt} = (tremSrndSS{tt} / params(1).sf) * 1000;
        movSrndSS{tt} = cat(2,movSrndSStimepts{tt}{:});
        movSrndSS{tt} = (movSrndSS{tt} / params(1).sf) * 1000;
        
        figure('Name','Tremor triggered neurophysiology')
        subplot(1,2,1); hold on
        histogram(tremSrndSS{tt},150)
        xlabel('Time (ms)')
        ylabel('Number of events')
        title('Tremor PSTH')
        set(gca,'FontSize',16)
        
        subplot(1,2,2); hold on
        histogram(movSrndSS{tt},150)
        xlabel('Time (ms)')
        title('Movement PSTH')
        set(gca,'FontSize',16)
    end
end

%% Analyze PSTH centered on filtered tremor epochs

for tt = 1:2
    if exist('filtTremIdxs')
        [tremSrndEphys{tt},~] = peristimData(waveformData(1).waveformTimepts{1},filtTremIdxs{tt},sig{2,1},params,2);
    end
    
    if exist('filtTremRel')
        if ~isempty(filtTremRel{tt})
            for i = 1:length(tremSrndEphys{tt})
                if filtTremIdxs{tt}(i,1) - 2999 < 1 || filtTremIdxs{tt}(i,2) + 3000 > length(sig{2,1})
                    continue
                else
                    tempFiltEphys = sig{2,1}(filtTremIdxs{tt}(i,1) - 2999:filtTremIdxs{tt}(i,2) + 3000);
                    tremTrigFiltEphys{tt,i} = circshift(tempFiltEphys,filtTremRel{tt})';
                    tremSrndEphysTimepts{tt,i} = tremSrndEphys{tt}{i} - filtTremRel{tt,i};
                end
            end
        end
    end
    
    % Create variable argutment in and out function
    % Shift data to align tremor peaks
    % Add tremor timepts to histogram
    
    if exist('filtTremIdxs')
        tremSrnd{tt} = cat(2,tremSrndEphys{tt}{:});
        tremSrnd{tt} = (tremSrnd{tt} / params(1).sf) * 1000;
        
        figure('Name','Tremor triggered neurophysiology')
        hold on
        histogram(tremSrnd{tt},50)
        xlabel('Time (ms)')
        ylabel('Number of events')
        title('Tremor PSTH')
        set(gca,'FontSize',16)
    end
end

%% Burst triggered tremor epochs

for tt = 1:length(ssBurstIdx)
    ssBurstTimepts{tt}(:,1) = waveformData(1).waveformTimepts{tt}(ssBurstIdx{tt}(:,1));
    ssBurstTimepts{tt}(:,2) = waveformData(1).waveformTimepts{tt}(ssBurstIdx{tt}(:,2));
    
    [tremSrndBurstSStimepts{tt},~] = peristimData(ssBurstTimepts{tt},sig{3,4},sig{3,4},params,3);
    
    figure('Name','Neurophys-triggered tremor'); hold on
        for i = 1:length(tremSrndBurstSStimepts{tt})
            plot(tremSrndBurstSStimepts{tt}{i})
        end
        xlabel('Time (ms)')
        ylabel('\muVs')
        title('Neurophys-triggered tremor epochs')
        set(gca,'FontSize',16)
    
    % Find rhythmic bursts
    ssBurstTimeptsMs{tt} = (ssBurstTimepts{tt} / params(1).sf) * 1000;
    a = 1;
    for i = 1:length(ssBurstTimeptsMs{tt})-1
        if ssBurstTimeptsMs{tt}(i,2) > (ssBurstTimeptsMs{tt}(i+1,1) - 75)
            tempRhythBurst{tt}(a) = i;
            a = a + 1;
        end
    end
    
    endIdxs = [];
    start = 1;
    ending = 1;
    count = 1;
    for i = 2:length(tempRhythBurst{tt})
     
        if tempRhythBurst{tt}(i) == tempRhythBurst{tt}(i-1) + 1
            ending = i;
        else
            endIdxs = [endIdxs; start, ending];
            start = i;
            ending = i;
        end
    
    end
    tempBurstPts{tt} = endIdxs(:,2) - endIdxs(:,1);
    burstPtsTemp{tt} = find(tempBurstPts{tt}>0);
    tempBurstIdxs{tt}(:,1) = endIdxs(burstPtsTemp{tt},1);
    tempBurstIdxs{tt}(:,2) = endIdxs(burstPtsTemp{tt},2);
    for i = 1:size(tempBurstIdxs{tt},1)
        burstIdxs{tt}(i,1) = tempRhythBurst{tt}(tempBurstIdxs{tt}(i,1));
        burstIdxs{tt}(i,2) = tempRhythBurst{tt}(tempBurstIdxs{tt}(i,2));
        rhythmicBurstTimepts{tt}(i,1) = ssBurstTimepts{tt}(burstIdxs{tt}(i,1),1);
        rhythmicBurstTimepts{tt}(i,2) = ssBurstTimepts{tt}(burstIdxs{tt}(i,2),2);
    end
    
    figure(); hold on
    for i = 1:size(rhythmicBurstTimepts{tt},1)
        plot(sig{2,2}(rhythmicBurstTimepts{tt}(i,1):rhythmicBurstTimepts{tt}(i,2)))
    end
end

%% Selective burst-triggered tremor epochs

temp2ndSpike = {}; temp2ndSpikePost = {};
for tt = 1:length(waveformData(1).waveformTimepts)
    a = 1; b = 1; c = 1; d = 1; e = 1;
    for i = 1:length(ssBurstIdx{tt})
        if waveformData(1).waveformTimepts{tt}(ssBurstIdx{tt}(i,1) - 1) >= 1 && waveformData(1).waveformTimepts{tt}(ssBurstIdx{tt}(i,2) + 1) <= length(sig{2,1})
            tempSpike(1) = (waveformData(1).waveformTimepts{tt}(ssBurstIdx{tt}(i,1) - 1)) / params(1).sf * 1000;
            tempSpike(2) = (waveformData(1).waveformTimepts{tt}(ssBurstIdx{tt}(i,2) + 1)) / params(1).sf * 1000;
            tempSpike(3) = (waveformData(1).waveformTimepts{tt}(ssBurstIdx{tt}(i,1))) / params(1).sf * 1000;
            tempSpike(4) = (waveformData(1).waveformTimepts{tt}(ssBurstIdx{tt}(i,2))) / params(1).sf * 1000;

            if tempSpike(3) - tempSpike(1) < 70 && tempSpike(3) - tempSpike(1) > 30
                rhythmSpike{tt,1}(a) = i;
                a = a + 1;
            end
            if tempSpike(2) - tempSpike(4) < 70 && tempSpike(2) - tempSpike(4) > 30
                rhythmSpike{tt,2}(b) = i;
                b = b + 1;
            end

            if tempSpike(3) - tempSpike(1) < 70 && tempSpike(3) - tempSpike(1) > 30 && tempSpike(2) - tempSpike(4) < 70 && tempSpike(2) - tempSpike(4) > 30
                rhythmSpike{tt,3}(c) = i;
                c = c + 1;

                temp2ndSpike{tt,i} = find(waveformData(1).waveformTimeptsMS{tt} < (tempSpike(1) - 30) & waveformData(1).waveformTimeptsMS{tt} > (tempSpike(1) - 80));
                temp2ndSpikePost{tt,i} = find(waveformData(1).waveformTimeptsMS{tt} > (tempSpike(2) + 30) & waveformData(1).waveformTimeptsMS{tt} < (tempSpike(2) + 80));
                
                if ~isempty(temp2ndSpike{tt,i})
                    secSpikePre(tt,d) = i;
                    d = d + 1;
                end

                if ~isempty(temp2ndSpikePost{tt,i})
                    secSpikePost(tt,e) = i;
                    e = e + 1;
                end

            end
        end
    end
    secSpikePre(secSpikePre==0) = NaN;
    secSpikePost(secSpikePost==0) = NaN;

    rhythmSpike{tt,4} = intersect(secSpikePre(tt,:),secSpikePost(tt,:));

    for i = 1:length(rhythmSpike{tt,4})
        rhythmSpike{tt,5}(i,1) = ssBurstTimepts{tt}(rhythmSpike{tt,4}(i),1);
        rhythmSpike{tt,5}(i,2) = ssBurstTimepts{tt}(rhythmSpike{tt,4}(i),2);

        tempWfmIdx = find(rhythmSpike{tt,5}(i,1)==waveformData(1).waveformTimepts{tt});
        rhythmSpike{tt,6}(i,1) = waveformData(1).waveformTimepts{tt}(tempWfmIdx - 2) - 750;
        rhythmSpike{tt,6}(i,2) = waveformData(1).waveformTimepts{tt}(tempWfmIdx + 2) + 750;

        data2export.tremEpochs{tt,10}{i} = data{6}(rhythmSpike{tt,6}(i,1):rhythmSpike{tt,6}(i,2));
    end


    endIdxs = [];
    start = 1;
    ending = 1;
    count = 1;
    for i = 2:length(rhythmSpike{tt,4}) - 1
     
        if rhythmSpike{tt,4}(i) == rhythmSpike{tt,4}(i-1) + 1
            ending = i;
        else
            endIdxs = [endIdxs; start, ending];
            start = i;
            ending = i;
        end
    
    end

    tempBursts{tt} = endIdxs(:,2) - endIdxs(:,1);
    burstPts{tt} = find(tempBursts{tt}>0);
    tempBurstTimes{tt}(:,1) = endIdxs(burstPts{tt},1);
    tempBurstTimes{tt}(:,2) = endIdxs(burstPts{tt},2);
    for i = 1:length(tempBurstTimes{tt})
        filtBurstIdxs{tt}(i,1) = rhythmSpike{tt,4}(tempBurstTimes{tt}(i,1));
        filtBurstIdxs{tt}(i,2) = rhythmSpike{tt,4}(tempBurstTimes{tt}(i,2));
        rhythmicBurstIdxs{tt}(i,1) = ssBurstTimepts{tt}(filtBurstIdxs{tt}(i,1),1);
        rhythmicBurstIdxs{tt}(i,2) = ssBurstTimepts{tt}(filtBurstIdxs{tt}(i,2),2);
    end
    figure()
    t = tiledlayout('flow','TileSpacing','tight');
    for i = 1:length(rhythmSpike{tt,6})
        nexttile; hold on
        plot(sig{2,2}(rhythmSpike{tt,6}(i,1):rhythmSpike{tt,6}(i,2)))
        plot(data{6}(rhythmSpike{tt,6}(i,1):rhythmSpike{tt,6}(i,2)))
    end
    title(t,'Neurophys-triggered tremor','FontSize',16)
    xlabel(t,'Timepoints','FontSize',16)
    ylabel(t,'\muVs','FontSize',16)

    if length(rhythmicBurstIdxs) == 1
        tt = 1;
    else

        figure()
        t = tiledlayout('flow','TileSpacing','tight');
        for i = 1:length(rhythmicBurstIdxs{tt})
            nexttile; hold on
            plot(sig{2,2}(rhythmicBurstIdxs{tt}(i,1):rhythmicBurstIdxs{tt}(i,2)))
            plot(data{6}(rhythmicBurstIdxs{tt}(i,1):rhythmicBurstIdxs{tt}(i,2)))
        end
        title(t,'Neurophys-triggered tremor','FontSize',16)
        xlabel(t,'Timepts','FontSize',16)
        ylabel(t,'\muVs','FontSize',16)
    end
end

%% Find spike times for burst-triggered tremor epochs

for ii = 1:length(waveformData(1).waveformTimepts)
    [filtBurstSpikes{ii},~] = peristimData(waveformData(1).waveformTimepts{ii},rhythmSpike{1,6},sig{2,1},params,4);
    for i = 1:length(filtBurstSpikes{1})
        filtBurstSpikeTimes{ii,i} = waveformData(1).waveformTimepts{1}(filtBurstSpikes{ii}{i});
        burstSpikeTimes{ii,i} = filtBurstSpikeTimes{ii,i} - rhythmSpike{1,6}(i,1);
        data2export.tremEpochs{1,10}{2,i} = data2export.tremEpochs{1,10}{ii,i}*-1;
    end
end

%% Find peaks of each tremor event

proceed = 'n';
while proceed == 'n'
    for ii = 1:2
        % Parameters for finding peaks
        prompt = {'minPeakProm: ','minPeakInt: ','maxPeakWidth: ','minPeakWidth: '};
        dlgtitle = 'FindPeaks params';
        if ii == 1
            default = {'4','500','3000','50'};
        elseif ii == 2
            default = {'4','500','3000','50'};
        end
        tempParams = inputdlg(prompt,dlgtitle,1,default);
        parameters(1).minPeakProm = str2double(tempParams{1});
        parameters(1).minPeakDist =  str2double(tempParams{2});
        parameters(1).maxPeakWidth =  str2double(tempParams{3});
        parameters(1).minPeakWidth =  str2double(tempParams{4});
        for i = 1:length(filtBurstSpikeTimes)
            [pks,idxs,w,p] = findpeaks(data2export.tremEpochs{1,10}{ii,i},1,'MinPeakProminence',parameters(1).minPeakProm,'MinPeakDistance',parameters(1).minPeakDist,'Annotate','extents','MaxPeakWidth',parameters(1).maxPeakWidth,'MinPeakWidth',parameters(1).minPeakWidth);
            invPeaks{ii,i} = pks; invLocs{ii,i} = idxs; invWidths{ii,i} = w; invOtherPeaks{ii,i} = p;
            % For trough detection: Flip the sign of each tremorEpoch and rerun
            % peak detection to get locs and peaks
        end
        clear pks idxs w p
    
        fig = figure();
        for i = 1:20
            nexttile; hold on
            plot(data2export.tremEpochs{1,10}{ii,i},'LineWidth',2,'Color','k')
            peakIdx = invLocs{ii,i}; peakAmp = invPeaks{ii,i};
            for m = 1:length(invLocs{ii,i})
                scatter(peakIdx(m),peakAmp(m),100,'LineWidth',2)
            end
        end
    
        prompt = {'Keep peaks? (y/n): '};
        dlgtitle = 'Peak analysis';
        proceed = inputdlg(prompt,dlgtitle);
        proceed = proceed{1};
    end
end

%% 
counter = [1,1;1,1];
for tt = 1
    for i = 1:length(burstSpikeTimes)
        for ii = 1:length(burstSpikeTimes{tt,i})
            if isempty(invPeaks{tt,i}(:))
                burstTremDiff{tt,i} = NaN;
            else
                burstTremDiff{tt,i}(ii,:) = burstSpikeTimes{tt,i}(ii) - invLocs{1,i}(:);
            end
            if isempty(invPeaks{2,i})
                burstTremDiff2{tt,i}(ii,:) = NaN;
            else
                burstTremDiff2{tt,i}(ii,:) = burstSpikeTimes{tt,i}(ii) - invLocs{2,i}(:);
            end
        end
        if isempty(burstTremDiff{tt,i}) || isempty(burstTremDiff2{tt,i})
            continue
        else
            for k = 1:size(burstTremDiff{tt,i},1)
                [burstTremTemp{tt,i}(k),burstTremTempIdx{tt,i}(k)] = min(abs(burstTremDiff{tt,i}(k,:)));
            end
        
            for jj = 1:size(burstTremDiff2{tt,i},1)
                [burstTremTemp2{tt,i}(jj),burstTremTempIdx2{tt,i}(jj)] = min(abs(burstTremDiff2{tt,i}(jj,:)));
            end
        
            for mm = 1:length(burstTremTemp{tt,i})
                if burstTremTemp{tt,i}(mm) < burstTremTemp2{tt,i}(mm)
                    burstTremTime{tt,i}(mm) = burstTremDiff{tt,i}(mm,burstTremTempIdx{tt,i}(mm));
                    counter(tt,1) = counter(tt,1) + 1;
                else

                    burstTremTime{tt,i}(mm) = burstTremDiff2{tt,i}(mm,burstTremTempIdx2{tt,i}(mm));
                    counter(tt,2) = counter(tt,2) + 1;
                end
            end
        end
    end
    burstTremTiming{tt} = cat(2,burstTremTime{tt,:});
    burstTremTimingMS{tt} = (burstTremTiming{tt} / params(1).sf) * 1000;

    
    figure(); hold on
    histogram(burstTremTimingMS{tt},200)
    xlabel('Time of AP re: peak or trough of oscillation(ms)')
    ylabel('Num. events')
    set(gca,'FontSize',16)

end

minLoc = 1000;
for i = 1:length(invLocs)
    if isempty(invLocs{1,i})
        temp1stPk(i) = 3000;
    else
        temp1stPk(i) = min(invLocs{1,i});
        if temp1stPk(i) < minLoc
            minLoc = temp1stPk(i);
            minLocIdx = i;
        end
    end
end

for i = 1:length(invLocs)
    if exist('minLocIdx')
        if i ~= minLocIdx
            shiftTime = minLoc - temp1stPk(i);
            tempShift = circshift(data2export.tremEpochs{1,10}{1,i},shiftTime);
            data2export.tremEpochs{1,11}{1,i} = tempShift;
            data2export.tremEpochs{1,11}{2,1}(i) = shiftTime;
        end
    end
end

figure(); hold on
for i = 1:length(invLocs)
    plot(data2export.tremEpochs{1,11}{1,i})
end

if exist('invLocs')
    [burstTremorBurstSrnd{1},~] = peristimData(waveformData(1).waveformTimepts{1},rhythmSpike{1,6},sig{2,1},params,2);
end
    
if exist('data2export.tremEpochs')
    if ~isempty(data2export.tremEpochs{1,10})
        for i = 1:length(burstTremorBurstSrnd{1})
            if rhythmSpike{1,6}(i,1) - 2999 < 1 || rhythmSpike{1,6}(i,2) + 3000 > length(sig{2,1})
                continue
            else
                burstTremorBurstShift{1,i} = burstTremorBurstSrnd{1}{i} - data2export.tremEpochs{1,11}{2,1}(i);
            end
        end
    end
end

% Create variable argutment in and out function
% Shift data to align tremor peaks
% Add tremor timepts to histogram

if exist('rhythmSpike')
    brstTremBrst{1} = cat(2,burstTremorBurstShift{1,:});
    brstTremBrst{1} = (brstTremBrst{1} / params(1).sf) * 1000;
    figure('Name','Tremor triggered neurophysiology')
    hold on
    histogram(brstTremBrst{1},1000)
    xlabel('Time (ms)')
    ylabel('Number of events')
    title('Tremor PSTH')
    xlim([-150 150])
    set(gca,'FontSize',16)

end

%% Vector with non-consecutive burst event indices
% filtBurstIdxs
for tt = 1:length(waveformData(1).waveformTimepts)
    a = 1;
    for i = rhythmSpike{tt,4}
        for ii = 1:length(filtBurstIdxs{tt})
            if i > filtBurstIdxs{tt}(ii,1) && i <= filtBurstIdxs{tt}(ii,2)
                tempRemove{tt}(a) = i;
                tempI{tt}(a) = find(tempRemove{tt}(a) == rhythmSpike{tt,4});
                a = a + 1;
            end
        end
    end

    rhythmSpike{tt,4}(tempI{tt}) = [];
    for i = 1:length(rhythmSpike{4})
        rhythmSpike{tt,7}(i,1) = ssBurstTimepts{tt}(rhythmSpike{tt,4}(i),1);
        rhythmSpike{tt,7}(i,2) = ssBurstTimepts{tt}(rhythmSpike{tt,4}(i),2);
    end
    
    [burstSrndTremor{tt},~] = peristimData(rhythmSpike{tt,6},data{6},data{6},params,3);
    [filtBurstSpikes{tt},~] = peristimData(waveformData(1).waveformTimepts{tt},rhythmSpike{tt,7},sig{2,1},params,4);

end

%% 

summary = cell(13,1);
summary(:,1) = {'ffreq','peakIsi','cv','cv2','burstPauseAvg','burstFreq','PeakTroughCounts','meanPSD','maxPSD','maxPSDidx','normPSD','normPSDmax','nPSDmaxIdx'};
for i = 1
    summary(:,i+1)={data2export.firingFreq,data2export.isiPkIdx,data2export.cv,data2export.cv2,data2export.burstPauseAvg,data2export.burstFreq,counter,data2export.mPSD,data2export.maxPSD,data2export.maxPSDidx,data2export.nPSD,data2export.nPSDmax,data2export.nPSDmaxIdx};
end
writecell(summary,[params(1).recName '.xlsx'],'Sheet',1,'Range','A1');

save([params(1).recName])

%% Raster plot

% dimXmax = 0;
% for i = 1:length(tremSrndSStimepts)
%     dimXtemp = length(tremSrndSStimepts{i});
%     if dimXtemp > dimXmax
%         dimXmax = dimXtemp;
%     end
% end
% rasterArr = ones(length(tremSrndSStimepts),dimXmax);
% for i = 1:length(tremSrndSStimepts)
%     rasterArr(i,1:length(tremSrndSStimepts{i})) = tremSrndSStimepts{i}';
% end
% 
% 
% rasterplot(rasterArr,49,10000)

%%  Tremor epoch summary statistics

% figure('Name','Peak Amplitudes'); hold on
% histogram(tremor{6},30)
% xlabel('Peak amplitude (\muVs)')
% ylabel('Number of events')
% title('Peak amplitudes')
% set(gca,'FontSize',14)
% 
% figure('Name','Num peaks per epoch'); hold on
% histogram(tremorDur)
% xlabel('Number of peaks')
% ylabel('Number of epochs')
% title('Number of peaks per epoch')
% set(gca,'FontSize',14)
% 
% figure('Name','Tremor frequencies'); hold on
% histogram(tremFreqActual,11)
% xlabel('Tremor frequency')
% ylabel('Number of epochs')
% title('Tremor frequencies')
% set(gca,'FontSize',14)

%% Compute coherence

% coh = NaN(nCoh/2+1,endNumber);
% [Cy1y2,freqCoh] = mscohere(y1,y2,hanning(nCoh),1/2*nCoh,nCoh,parameters(num).sf); % (y1,y2 coherece, with hamming window, take NCoh points, half points overlap, number of fft:NCoh, sampling frequency:fs
% coh(:,s) = Cy1y2;
