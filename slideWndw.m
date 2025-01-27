%   Creates upper and lower bounds of data using percentiles of a sliding
%   window
%
%   Alex Fanning, August 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f0,g0] = slideWndw(parameters,data)

N = length(data);

lowRange(1) = parameters(1).diffRange(1) - 5;
lowRange(2) = parameters(1).diffRange(1) + 5;
highRange(1) = parameters(1).diffRange(2) - 5;
highRange(2) = parameters(1).diffRange(2) + 5;

for i = 1:N
    dom = i-parameters(1).wndwSize/2:i+parameters(1).wndwSize/2;
    dom = dom(find(dom>=1 & dom<=N));
    pc = data(dom);
    mi = prctile(pc,lowRange(1));
    ma = prctile(pc,lowRange(2));
    ti = prctile(pc,highRange(1));
    ta = prctile(pc,highRange(2));
    id = find(pc>mi & pc<ma);
    tid = find(pc>ti & pc<ta);
    f0(i) = median(pc(id));
    g0(i) = median(pc(tid));
end