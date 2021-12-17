DTs=DTsAll;

hemispheres=["Left", "Right"];
for whichDT=1:length(DTs)
    if DTs(whichDT)>0
        for whichFolder=1:length(modelFolders)
            whichModelNames=modelNames{whichFolder};
            %whichModelFiles=modelFileNames{whichFolder};
            for whichModel=1:length(whichModelNames)
                for whichMap=1:length(mapNames)
                    if eval(char(strcat(['isfield(',thisData, ', ''',char(mapNames(whichMap)),''')'])))
                        for whichHemi=1:length(hemispheres)
                            targetData=char(strcat(thisData, '.All.',hemispheres(whichHemi),'.', whichModelNames(whichModel), '.', DTnames(whichDT)));
                            %                             if ~exist('dataAllSub', 'var')
                            %                                 eval([targetData, '=[];'])
                            %                             end
                            
                            tmp=char(strcat(thisData, '.', mapNames(whichMap)));
                            if eval(char(strcat(['isfield(',tmp, ', ''',char(hemispheres(whichHemi)),''')'])))
                                sourceData=char(strcat(thisData, '.', mapNames(whichMap), '.',hemispheres(whichHemi),'.', whichModelNames(whichModel), '.', DTnames(whichDT)));
                                
                                %Check and initialise variables
                                try
                                    eval(['~isfield(',targetData, ', ''mapOrder'');'])
                                catch
                                    eval([targetData, '=[];'])
                                end
                                
                                if eval(['~isfield(',targetData, ', ''mapOrder'')'])
                                    eval(char(strcat(targetData, '.mapOrder= "', mapNames(whichMap),'";')));
                                    eval([targetData, '.mapIndices=', num2str(length(eval([sourceData '.x0s']))), ';']);
                                else
                                    %Fill variables
                                    eval(char(strcat(targetData, '.mapOrder = [', targetData, '.mapOrder, "', mapNames(whichMap), '"];')))
                                    try
                                        eval([targetData, '.mapIndices = [', targetData, '.mapIndices, ', num2str(length(eval([sourceData '.x0s']))), '];'])
                                    catch
                                        targetData;
                                    end
                                    %eval([targetData, '.meanDist = [', targetData, '.meanDist, ', sourceData, '.meanDist];'])
                                end
                                
                                if eval(['~isfield(',thisData, ', ''All'')']) || eval(['~isfield(',targetData, ', ''x0s'')'])
                                    eval([targetData, '.x0s=',sourceData, '.x0s;']);
                                    eval([targetData, '.y0s=',sourceData, '.y0s;']);
                                    eval([targetData, '.sigmas=',sourceData, '.sigmas;']);
                                    eval([targetData, '.sigmaMinor=',sourceData, '.sigmaMinor;']);
                                    eval([targetData, '.sigmaTheta=',sourceData, '.sigmaTheta;']);
                                    eval([targetData, '.meanSignal=',sourceData, '.meanSignal;']);
                                    eval([targetData, '.roiIndices=',sourceData, '.roiIndices;']);
                                    eval([targetData, '.ves=',sourceData, '.ves;']);
                                    eval([targetData, '.rss=',sourceData, '.rss;']);
                                    eval([targetData, '.rawrss=',sourceData, '.rawrss;']);
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
                                    eval([targetData, '.rss = [', targetData, '.rss, ', sourceData, '.rss];'])
                                    eval([targetData, '.rawrss = [', targetData, '.rawrss, ', sourceData, '.rawrss];'])
                                    if eval(['isfield(',sourceData, ', ''exp'')'])
                                        eval([targetData, '.exp=[' targetData, '.exp, ',sourceData, '.exp];']);
                                    end
                                end
                                
                                if ismember(DTsAll(whichDT), allXvalDTs)
                                    if eval(['~isfield(',targetData, ', ''vesXval'')'])
                                        eval([targetData, '.vesXval=',sourceData, '.vesXval;']);
                                        eval([targetData, '.rssXval=',sourceData, '.rssXval;']);
                                        eval([targetData, '.rawrssXval=',sourceData, '.rawrssXval;']);
                                    else
                                        eval([targetData, '.vesXval = [', targetData, '.vesXval, ', sourceData, '.vesXval];'])
                                        eval([targetData, '.rssXval = [', targetData, '.rssXval, ', sourceData, '.rssXval];'])
                                        eval([targetData, '.rawrssXval = [', targetData, '.rawrssXval, ', sourceData, '.rawrssXval];'])
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
                                else
                                    %Fill variables
                                    eval(char(strcat(targetData, '.mapOrder = [', targetData, '.mapOrder, "', mapNames(whichMap), '"];')))
                                    eval([targetData, '.mapIndices = [', targetData, '.mapIndices, ', num2str(0), '];'])
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
