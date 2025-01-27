

%% Import sorted data

% Get sorted unit data
% sortedStruct = load(uigetfile('*.mat'));
% for i = 1:length(sortedStruct.channels.data)
%     sortedTimepts{i} = double(sortedStruct.channels.data(i).timestamps);
% end

sortedData = xlsread(uigetfile('*.xls'));
sortedData = sortrows(sortedData,2);
for i = 1:length(sortedData)
    indSpks(i,:) = sortedData(i,7:end);
end
indSpks = indSpks';

plot(indSpks)