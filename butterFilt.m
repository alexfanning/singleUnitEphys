%   Filter ephys data using Butterworth 4th order
%
%   Alex Fanning, 03/16/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outSig = butterFilt(parameters,data,range,passSet,figNum)

[y, x] = butter(4, range/(parameters(1).sf/2),passSet);
outSig = filter(y, x, data);

figure(figNum); hold on
plot(data)
plot(outSig)
set(gca,'FontSize',16)
saveas(figNum,[parameters(1).recName 'FilteredSig.tif'])
