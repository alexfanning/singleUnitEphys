function [newSig,segEndPts] = removeArtifact(signal,artCh,parameters) 

%   Removes artifact timepoints from neurophys data
%
%
%   Alex Fanning, July 2023
% *************************************************************************

prompt = {'Window size (ms): ', 'diff lower limit: ', 'diff upper limit: '};
dlgtitle = 'Sliding window';
default = {'200','5','95'};
paramsTemp = inputdlg(prompt,dlgtitle,1,default);
parameters(1).wndwSize = str2double(paramsTemp{1});
parameters(1).diffRange = [str2double(paramsTemp{2}) str2double(paramsTemp{3})];

% Loop through each sample to calculate f0
%%%%%% Use multiplication to enhance separation of artifact from baseline %%%%%%%%%%%%%%
[f0,g0] = slideWndw(parameters,artCh);

%% Calculate difference between bounds and thresholds
% Bounds difference and tremor start/stop threshold
bndDiff = NaN(1,length(f0));
bndDiff = g0 - f0;

% baseline threshold
bndDiffThr = 1.25 * nanmedian(bndDiff);

% Find threshold for movement event detection based on cumulative sum
movThr = 3 * nanmedian(bndDiff);

%% Plot raw data, bandpass vector, and difference of bounds

figure('Name','Bandpass data and bounds')
    ax = subplot(3,1,1); hold on
    plot(artCh)
    ylabel('\muVs')
    ylim([-2000 2000])
    set(gca,'FontSize',14)
    
    ax1 = subplot(3,1,2); hold on
    plot(artCh)
    plot(f0)
    plot(g0)
    ylim([-2000 2000])
    ylabel('\muVs')
    set(gca,'FontSize',14)
    
    ax2 = subplot(3,1,3); hold on
    plot(bndDiff)
    yline(bndDiffThr,'LineWidth',2,'Color','r')
    yline(movThr,'LineWidth',2)
    ylim([0 2000])
    linkaxes([ax ax1 ax2], 'x')
    xlabel('Time (ms)')
    ylabel('Difference')
    legend('','Tremor endpts threshold','Movement threshold')
    set(gca,'FontSize',14)

%% Plot histogram of difference of bounds with threshold indicated

figure('Name','Movement threshold'); hold on
    histogram(bndDiff,1000,'FaceColor','k')
    xline(movThr,'r','LineWidth',1)
    xlabel('Difference value')
    ylabel('Number of datapoints')
    title('Difference between bounds')
    xlim([0 2000])
    set(gca,'FontSize',14)

%% Extract potential tremor epoch start and end times

%%%%%% bndDiff < 5 or is NaN

% Find points above movement threshold
tempVec = zeros(1,length(bndDiff));
for i = 1:length(bndDiff)
    if bndDiff(i) > movThr
        tempVec(i) = 1;
    end
end

% Find start and end times for tremor epochs
a = 1;
for j = 2:length(tempVec)
    if tempVec(j) == 1 && tempVec(j-1) == 0
        tempTimepts(a,1) = j;
        
        % Find start timept
        for k = 1:j
            if bndDiff(k) < bndDiffThr
                startTemp = k;
            end
        end
        
        if exist('startTemp')
            segStart(a) = startTemp;
        else
            continue
        end

        % Find stop timept
        for ii = j:length(bndDiff)
            if bndDiff(ii) < (bndDiffThr)
                segStop(a) = ii;
                break
            end
        end

        a = a + 1;
    end
end

if exist('segStart')
    segEndPts(:,1) = unique(segStart)';
    if length(unique(segStart)) ~= length(unique(segStop))
        segEndPts(end) = [];
    end
    segEndPts(:,2) = unique(segStop)';
    
    % Create cell array structure to hold tremor epochs
    for m = 1:size(segEndPts,1)
        allEpochs{m} = signal(segEndPts(m,1):segEndPts(m,2));
    end
    
    newSig = signal;
    
    for i = 1:length(allEpochs)
        newSig(segEndPts(i,1):segEndPts(i,2)) = NaN;
    end
else
    newSig = signal;
end

% Convert difference of bounds values that are NaN or less than 5 to signal
a = 1;
for i = 1:length(bndDiff)
    if isnan(bndDiff(i)) || bndDiff(i) <= 5
        newSig(i) = NaN;
        tempSeqList(a) = i;
        a = a + 1;
    end
end

% Find timepoints of NaNs or values less than 5
if exist('tempSeqList','var')
    endpts = [];
    n = length(tempSeqList);
    start = tempSeqList(1);
    ending = tempSeqList(1);
    
    for i = 2:n
        if tempSeqList(i) == ending + 1
            ending = tempSeqList(i);
        else
            endpts = [endpts; start, ending];
            start = tempSeqList(i);
            ending = tempSeqList(i);
        end
    end
end

%%%%%%%%    REMOVE DATA WITHOUT ENOUGH SEQUENTIAL TIMEPOINTS %%%%%%%%%%%%%%

% Combine start and end timepoints for data that was removed
if exist('tempSeqList','var')
    endpts = [endpts; start, ending];
    if exist('segEndPts')
        segEndPts = sort(cat(1,endpts,segEndPts));
    else
        segEndPts = endpts;
    end
end

figure('Name','Filtered neurophys trace'); hold on
    plot(signal)
    plot(newSig)
    ylim([-6000 6000])
    xlabel('Timepoints')
    ylabel('\muVs')
    legend('Raw','Filtered')
