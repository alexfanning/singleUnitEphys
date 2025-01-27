%   Extract data around stimulus
%
%   Alex Fanning, 07/24/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [stimSrndTimepts,stimPause] = peristimData(waveTimepts,stimIdx,sig,parameters,type)

if type == 1
    for i = 1:length(stimIdx)
    
        stimPause(i) = (waveTimepts(stimIdx(i,2) + 1) - waveTimepts(stimIdx(i,2))) / (parameters(1).sf / 1000);
    
        if waveTimepts(stimIdx(i,1)) < 3000 || waveTimepts(stimIdx(i,2)) + 6000 > length(sig)
            continue
        else
            timeWndw(1) = waveTimepts(stimIdx(i,1)) - 3000;
            temp1 = waveTimepts > timeWndw(1);
            temp2 = waveTimepts < waveTimepts(stimIdx(i,1));
            temp3 = find(temp1==1);
            temp4 = find(temp2==1);
            temp5 = temp3(ismember(temp3,temp4));
            
            timeWndw(2) = waveTimepts(stimIdx(i,2)) + 5999;
            temp6 = waveTimepts < timeWndw(2);
            temp7 = waveTimepts > waveTimepts(stimIdx(i,2));
            temp8 = find(temp6==1);
            temp9 = find(temp7==1);
            temp10 = temp8(ismember(temp8,temp9));
            
            tempTimepts2 = cat(1,waveTimepts(temp5),waveTimepts(temp10));
            
            for ii = 1:length(tempTimepts2)
    
                if tempTimepts2(ii) <= waveTimepts(stimIdx(i,1))
                    stimSrndIdxs(ii) = tempTimepts2(ii) - waveTimepts(stimIdx(i,1));
                elseif tempTimepts2(ii) > waveTimepts(stimIdx(i,1)) && tempTimepts2(ii) < waveTimepts(stimIdx(i,2))
                    stimSrndIdxs(ii) = [];
                elseif tempTimepts2(ii) >= waveTimepts(stimIdx(i,2))
                    stimSrndIdxs(ii) = tempTimepts2(ii) - waveTimepts(stimIdx(i,2));
                end
    
            end
    
        end
        
        stimSrndTimepts{i} = stimSrndIdxs;
    end
elseif type == 2

    for i = 1:length(stimIdx)
        if stimIdx(i,1) < 15000 || stimIdx(i,2) + 30000 > length(sig)
            continue
        else
            timeWndw(1) = stimIdx(i,1) - 15000;
            temp1 = waveTimepts > timeWndw(1);
            temp2 = waveTimepts < stimIdx(i,1);
            temp3 = find(temp1==1);
            temp4 = find(temp2==1);
            temp5 = temp3(ismember(temp3,temp4));
            
            timeWndw(2) = stimIdx(i,1) + 29999;
            temp6 = waveTimepts < timeWndw(2);
            temp7 = waveTimepts > stimIdx(i,1);
            temp8 = find(temp6==1);
            temp9 = find(temp7==1);
            temp10 = temp8(ismember(temp8,temp9));
            
            tempTimepts2 = cat(1,waveTimepts(temp5),waveTimepts(temp10));

            for ii = 1:length(tempTimepts2)
    
                    stimSrndIdxs(ii) = tempTimepts2(ii) - stimIdx(i,1);
    
            end
        end
        stimSrndTimepts{i} = stimSrndIdxs;

    end
    stimPause = 0;
elseif type == 3
    for i = 1:length(waveTimepts)
        if waveTimepts(i,1) < 7500 || waveTimepts(i,2) + 14999 > length(sig)
            continue
        else
            timeWndw(1) = waveTimepts(i,1) - 3000;
            timeWndw(2) = waveTimepts(i,2) + 14999;
            a = 1;
            for t = timeWndw(1):timeWndw(2)
                temp(a) = stimIdx(t);
                a = a + 1;
            end

            stimSrndTimepts{i} = temp;
            clear temp
        end

    end
    stimPause = 0;
elseif type == 4
    for i = 1:length(stimIdx)
        temp1 = waveTimepts >= stimIdx(i,1);
        temp2 = waveTimepts <= stimIdx(i,2);
        temp3 = find(temp1==1);
        temp4 = find(temp2==1);
        temp5 = temp3(ismember(temp3,temp4));
            
        stimSrndTimepts{i} = temp5;

    end
    stimPause = 0;
end
