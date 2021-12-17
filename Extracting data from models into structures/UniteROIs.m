mapNames=["TTO", "TLO", "TPO", "TPC1", "TPC2", "TPC3", "TLS", "TF1", "TF2"];
modelProgName=["Log", "Lin"];%, "LinLog"];
modelSpace=["DurationPeriod", "OccupancyPeriod", "OnTimeOffTime"];
DTnames=["All", "Odd", "Even"];
hemispheres=["Left", "Right"];
shapes=["Oval", "Circ", "OvalFreeExponent", "OvalFreeExponentHRF"];

for whichDT=1:length(DTname
    for whichSpace=1:length(modelSpace)
        for whichProg=1:length(modelProgName)
            if whichSpace==3 && whichProg==3
                %Do nothing
            else
                if whichProg==3
%                     if whichSpace==1 %Remove when models are done
%                         nShapes=[1 3];
%                     elseif whichSpace==2 %Remove
                        nShapes=1; 
%                     end                
                elseif whichSpace==3 && whichProg==2 
                    nShapes=[1:3];
                elseif whichSpace==3
                    nShapes=[1 2];
                elseif whichSpace==2 && whichProg==2 
                    nShapes=[1:3];
                elseif whichSpace==1 && whichProg==1 
                    nShapes=[1:3];
                elseif whichSpace==1 && whichProg==2  
                    nShapes=[1:4];
                else
                    nShapes=1:2;    %1:4
                end
                for whichShape=nShapes
                    
                    for whichMap=1:length(mapNames)
                        if eval(char(strcat(['isfield(',thisData, ', ''',char(mapNames(whichMap)),''')'])))
                            for whichHemi=1:length(hemispheres)
                                targetData=char(strcat(thisData, '.All.',hemispheres(whichHemi),'.', modelSpace(whichSpace), '.', modelProgName(whichProg),'.', shapes(whichShape), '.', DTnames(whichDT)));
                                %                             if ~exist('dataAllSub', 'var')
                                %                                 eval([targetData, '=[];'])
                                %                             end
                                
                                tmp=char(strcat(thisData, '.', mapNames(whichMap)));
                                if eval(char(strcat(['isfield(',tmp, ', ''',char(hemispheres(whichHemi)),''')'])))
                                    sourceData=char(strcat(thisData, '.', mapNames(whichMap), '.',hemispheres(whichHemi),'.', modelSpace(whichSpace), '.', modelProgName(whichProg),'.', shapes(whichShape), '.', DTnames(whichDT)));
                                    
                                    %Check and initialise variables
                                    try
                                        eval(['~isfield(',targetData, ', ''mapOrder'');'])
                                    catch
                                        eval([targetData, '=[];'])
                                    end
                                    
                                    if eval(['~isfield(',targetData, ', ''mapOrder'')'])
                                        eval(char(strcat(targetData, '.mapOrder= "', mapNames(whichMap),'";')));
                                        eval([targetData, '.mapIndices=', num2str(length(eval([sourceData '.ratio']))), ';']);
                                        eval([targetData, '.meanDist=25;']);%,sourceData, '.meanDist;']);
                                    else
                                        %Fill variables
                                        eval(char(strcat(targetData, '.mapOrder = [', targetData, '.mapOrder, "', mapNames(whichMap), '"];')))
                                        try
                                        eval([targetData, '.mapIndices = [', targetData, '.mapIndices, ', num2str(length(eval([sourceData '.ratio']))), '];'])
                                        catch
                                            targetData;
                                        end
                                        %eval([targetData, '.meanDist = [', targetData, '.meanDist, ', sourceData, '.meanDist];'])
                                    end
                                    
                                    if eval(['~isfield(',thisData, ', ''All'')']) || eval(['~isfield(',targetData, ', ''ratio'')'])
                                        eval([targetData, '.ratio=',sourceData, '.ratio;']);
                                        eval([targetData, '.x0s=',sourceData, '.x0s;']);
                                        eval([targetData, '.y0s=',sourceData, '.y0s;']);
                                        eval([targetData, '.sigmas=',sourceData, '.sigmas;']);
                                        eval([targetData, '.sigmaMinor=',sourceData, '.sigmaMinor;']);
                                        eval([targetData, '.sigmaTheta=',sourceData, '.sigmaTheta;']);
                                        eval([targetData, '.meanSignal=',sourceData, '.meanSignal;']);
                                        eval([targetData, '.roiIndices=',sourceData, '.roiIndices;']);
                                        eval([targetData, '.ves=',sourceData, '.ves;']);
                                        if eval(['isfield(',sourceData, ', ''exp'')'])
                                            eval([targetData, '.exp=',sourceData, '.exp;']);
                                        end
                                    else
                                        eval([targetData, '.ratio = [', targetData, '.ratio, ', sourceData, '.ratio];'])
                                        eval([targetData, '.x0s = [', targetData, '.x0s, ', sourceData, '.x0s];'])
                                        eval([targetData, '.y0s = [', targetData, '.y0s, ', sourceData, '.y0s];'])
                                        eval([targetData, '.sigmas = [', targetData, '.sigmas, ', sourceData, '.sigmas];'])
                                        eval([targetData, '.sigmaMinor = [', targetData, '.sigmaMinor, ', sourceData, '.sigmaMinor];'])
                                        eval([targetData, '.sigmaTheta = [', targetData, '.sigmaTheta, ', sourceData, '.sigmaTheta];'])
                                        eval([targetData, '.meanSignal = [', targetData, '.meanSignal, ', sourceData, '.meanSignal];'])
                                        eval([targetData, '.roiIndices = [', targetData, '.roiIndices, ', sourceData, '.roiIndices];'])
                                        eval([targetData, '.ves = [', targetData, '.ves, ', sourceData, '.ves];'])
                                        if eval(['isfield(',sourceData, ', ''exp'')'])
                                            eval([targetData, '.exp=[' targetData, '.exp, ',sourceData, '.exp];']);
                                        end
                                    end
                                    
                                    if whichDT>1
                                        if eval(['~isfield(',targetData, ', ''vesXval'')'])
                                            eval([targetData, '.vesXval=',sourceData, '.vesXval;']);
                                        else
                                            eval([targetData, '.vesXval = [', targetData, '.vesXval, ', sourceData, '.vesXval];'])
                                        end
                                    end
                                    
                                else
                                    try
                                        eval(['~isfield(',targetData, ', ''mapOrder'');'])
                                    catch
                                        eval([targetData, '=[];'])
                                    end
                                    
                                    if eval(['~isfield(',targetData, ', ''mapOrder'')'])
                                        eval(char(strcat(targetData, '.mapOrder= "', mapNames(whichMap),'";')));
                                        eval([targetData, '.mapIndices=', num2str(0),';'])
                                        eval([targetData, '.meanDist=nan;']);
                                    else
                                        %Fill variables
                                        eval(char(strcat(targetData, '.mapOrder = [', targetData, '.mapOrder, "', mapNames(whichMap), '"];')))
                                        eval([targetData, '.mapIndices = [', targetData, '.mapIndices, ', num2str(0), '];'])
                                        eval([targetData, '.meanDist = [', targetData, '.meanDist, nan];'])
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

%Monotonic
modelNames=["Flat", "LinDF", "LinDCompF", "CompDCompF", "GausDLinF", "GausDCompF", "GausOLinF", "GausOCompF", "TemporalFreq","TemporalFreqPerEvent"]; 
for whichDT=1:length(DTnames)
    for whichModel=1:length(modelNames)
        for whichMap=1:length(mapNames)
            if eval(char(strcat(['isfield(',thisData, ', ''',char(mapNames(whichMap)),''')'])))
                for whichHemi=1:length(hemispheres)
                    targetData=char(strcat(thisData, '.All.',hemispheres(whichHemi),'.Monotonic.', modelNames(whichModel), '.', DTnames(whichDT)));
                    %                             if ~exist('dataAllSub', 'var')
                    %                                 eval([targetData, '=[];'])
                    %                             end
                    
                    tmp=char(strcat(thisData, '.', mapNames(whichMap)));
                    if eval(char(strcat(['isfield(',tmp, ', ''',char(hemispheres(whichHemi)),''')'])))
                        sourceData=char(strcat(thisData, '.', mapNames(whichMap), '.',hemispheres(whichHemi),'.Monotonic.', modelNames(whichModel), '.', DTnames(whichDT)));
                        
                        %Check and initialise variables
                        try
                            eval(['~isfield(',targetData, ', ''mapOrder'');'])
                        catch
                            eval([targetData, '=[];'])
                        end
                        
                        if eval(['~isfield(',targetData, ', ''mapOrder'')'])
                            eval(char(strcat(targetData, '.mapOrder= "', mapNames(whichMap),'";')));
                            eval([targetData, '.mapIndices=', num2str(length(eval([sourceData '.ratio']))), ';']);
                            eval([targetData, '.meanDist=25;']);%,sourceData, '.meanDist;']);
                        else
                            %Fill variables
                            eval(char(strcat(targetData, '.mapOrder = [', targetData, '.mapOrder, "', mapNames(whichMap), '"];')))
                            try
                                eval([targetData, '.mapIndices = [', targetData, '.mapIndices, ', num2str(length(eval([sourceData '.ratio']))), '];'])
                            catch
                                targetData;
                            end
                            %eval([targetData, '.meanDist = [', targetData, '.meanDist, ', sourceData, '.meanDist];'])
                        end
                        
                        if eval(['~isfield(',thisData, ', ''All'')']) || eval(['~isfield(',targetData, ', ''ratio'')'])
                            eval([targetData, '.ratio=',sourceData, '.ratio;']);
                            eval([targetData, '.x0s=',sourceData, '.x0s;']);
                            eval([targetData, '.y0s=',sourceData, '.y0s;']);
                            eval([targetData, '.sigmas=',sourceData, '.sigmas;']);
                            eval([targetData, '.sigmaMinor=',sourceData, '.sigmaMinor;']);
                            eval([targetData, '.sigmaTheta=',sourceData, '.sigmaTheta;']);
                            eval([targetData, '.meanSignal=',sourceData, '.meanSignal;']);
                            eval([targetData, '.roiIndices=',sourceData, '.roiIndices;']);
                            eval([targetData, '.ves=',sourceData, '.ves;']);
                            if eval(['isfield(',sourceData, ', ''exp'')'])
                                eval([targetData, '.exp=',sourceData, '.exp;']);
                            end
                        else
                            eval([targetData, '.ratio = [', targetData, '.ratio, ', sourceData, '.ratio];'])
                            eval([targetData, '.x0s = [', targetData, '.x0s, ', sourceData, '.x0s];'])
                            eval([targetData, '.y0s = [', targetData, '.y0s, ', sourceData, '.y0s];'])
                            eval([targetData, '.sigmas = [', targetData, '.sigmas, ', sourceData, '.sigmas];'])
                            eval([targetData, '.sigmaMinor = [', targetData, '.sigmaMinor, ', sourceData, '.sigmaMinor];'])
                            eval([targetData, '.sigmaTheta = [', targetData, '.sigmaTheta, ', sourceData, '.sigmaTheta];'])
                            eval([targetData, '.meanSignal = [', targetData, '.meanSignal, ', sourceData, '.meanSignal];'])
                            eval([targetData, '.roiIndices = [', targetData, '.roiIndices, ', sourceData, '.roiIndices];'])
                            eval([targetData, '.ves = [', targetData, '.ves, ', sourceData, '.ves];'])
                            if eval(['isfield(',sourceData, ', ''exp'')'])
                                eval([targetData, '.exp=[' targetData, '.exp, ',sourceData, '.exp];']);
                            end
                        end
                        
                        if whichDT>1
                            if eval(['~isfield(',targetData, ', ''vesXval'')'])
                                eval([targetData, '.vesXval=',sourceData, '.vesXval;']);
                            else
                                eval([targetData, '.vesXval = [', targetData, '.vesXval, ', sourceData, '.vesXval];'])
                            end
                        end
                        
                    else
                        try
                            eval(['~isfield(',targetData, ', ''mapOrder'');'])
                        catch
                            eval([targetData, '=[];'])
                        end
                        
                        if eval(['~isfield(',targetData, ', ''mapOrder'')'])
                            eval(char(strcat(targetData, '.mapOrder= "', mapNames(whichMap),'";')));
                            eval([targetData, '.mapIndices=', num2str(0),';'])
                            eval([targetData, '.meanDist=nan;']);
                        else
                            %Fill variables
                            eval(char(strcat(targetData, '.mapOrder = [', targetData, '.mapOrder, "', mapNames(whichMap), '"];')))
                            eval([targetData, '.mapIndices = [', targetData, '.mapIndices, ', num2str(0), '];'])
                            eval([targetData, '.meanDist = [', targetData, '.meanDist, nan];'])
                        end
                    end
                end
            end
        end
        
    end
end