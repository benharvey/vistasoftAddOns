mapNames=["TLO", "TTOP", "TTOA", "TPO", "TLS", "TPCI", "TPCM", "TPCS", "TFI" "TFS"];
DTnames=["All", "Odd", "Even"];
hemispheres=["Left", "Right", "Both"];
modelNames=["DurationPeriod", "OccupancyPeriod", "OnTimeOffTime", "LogDurationPeriod","DurationPeriodHRF","DurationPeriodNoExponent","TemporalFreq", "Flat", "LinDF", "LinDCompF", "TemporalFreqPerEvent","GausOLinF", "GausDLinF", "GausOCompF", "GausDCompF",];
for whichDT=1:length(DTnames)
    for whichModel=1:length(modelNames)
        
        for whichMap=1:length(mapNames)
            if eval(char(strcat(['isfield(',thisData, ', ''',char(mapNames(whichMap)),''')'])))
                for whichHemi=1:2%length(hemispheres)
                    targetData=char(strcat('dataAllSub.', mapNames(whichMap), '.',hemispheres(whichHemi),'.', modelNames(whichModel), '.', DTnames(whichDT)));
                    %                             if ~exist('dataAllSub', 'var')
                    %                                 eval([targetData, '=[];'])
                    %                             end
                    
                    tmp=char(strcat(thisData, '.', mapNames(whichMap)));
                    if eval(char(strcat(['isfield(',tmp, ', ''',char(hemispheres(whichHemi)),''')'])))
                        sourceData=char(strcat(thisData, '.', mapNames(whichMap), '.',hemispheres(whichHemi),'.', modelNames(whichModel), '.', DTnames(whichDT)));
                        
                        %Check and initialise variables
                        try
                            eval(['~isfield(',targetData, ', ''subjectOrder'');'])
                        catch
                            eval([targetData, '=[];'])
                        end
                        
                        if eval(['~isfield(',targetData, ', ''subjectOrder'')'])
                            eval([targetData, '.subjectOrder=thisData;']);
                            eval([targetData, '.subjectIndices=', num2str(length(eval([sourceData '.ratio']))), ';']);
                            eval([targetData, '.meanDist=',sourceData, '.meanDist;']);
                        else
                            %Fill variables
                            eval([targetData, '.subjectOrder = [', targetData, '.subjectOrder, "', thisData, '"];'])
                            eval([targetData, '.subjectIndices = [', targetData, '.subjectIndices, ', num2str(length(eval([sourceData '.ratio']))), '];'])
                            eval([targetData, '.meanDist = [', targetData, '.meanDist, ', sourceData, '.meanDist];'])
                        end
                        
                        if eval(['~isfield(',targetData, ', ''ratio'')'])
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
                        
                        if whichDT>1 && whichModel~=5
                            if eval(['~isfield(',targetData, ', ''vesXval'')'])
                                eval([targetData, '.vesXval=',sourceData, '.vesXval;']);
                            else
                                eval([targetData, '.vesXval = [', targetData, '.vesXval, ', sourceData, '.vesXval];'])
                            end
                        end
                        
                    else
                        try
                            eval(['~isfield(',targetData, ', ''subjectOrder'');'])
                        catch
                            eval([targetData, '=[];'])
                        end
                        
                        if eval(['~isfield(',targetData, ', ''subjectOrder'')'])
                            eval([targetData, '.subjectOrder=thisData;'])
                            eval([targetData, '.subjectIndices=', num2str(0),';'])
                            eval([targetData, '.meanDist=nan;']);
                        else
                            %Fill variables
                            eval([targetData, '.subjectOrder = [', targetData, '.subjectOrder "', thisData, '"];'])
                            eval([targetData, '.subjectIndices = [', targetData, '.subjectIndices ', num2str(0), '];'])
                            eval([targetData, '.meanDist = [', targetData, '.meanDist, nan];'])
                        end
                    end
                end
            end
        end
    end
end