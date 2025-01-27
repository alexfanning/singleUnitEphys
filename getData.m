%   Extract data from Blackrock NS files
%
%   Written by Alex Fanning, 5/4/23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dataOut,parameters] = getData(parameters,dataOut,num)

NS = openNSx;
dataOut{num} = double(NS.Data');
parameters(num).chIDs = NS.MetaTags.ChannelID;
parameters(num).sf = NS.MetaTags.SamplingFreq;
parameters(1).recName = NS.MetaTags.Filename;