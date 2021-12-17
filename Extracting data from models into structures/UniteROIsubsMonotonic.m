DTs=DTsAll;


hemispheres=["Left", "Right"];
for whichDT=1:length(DTnames)
    for whichModel=1:length(modelNames)
        
        for whichMap=1:length(mapNames)
            if eval(char(strcat(['isfield(',thisData, ', ''',char(mapNames(whichMap)),''')'])))
                for whichHemi=1:2
                    targetData=char(strcat('dataAllSub.', mapNames(whichMap), '.',hemispheres(whichHemi),'.', modelNames(whichModel), '.', DTnames(whichDT)));
                    
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
                        else
                            %Fill variables
                            eval([targetData, '.subjectOrder = [', targetData, '.subjectOrder, "', thisData, '"];'])
                            eval([targetData, '.subjectIndices = [', targetData, '.subjectIndices, ', num2str(length(eval([sourceData '.ratio']))), '];'])
                        end
                        
                        if eval(['~isfield(',targetData, ', ''x0s'')'])
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
                            eval(['~isfield(',targetData, ', ''subjectOrder'');'])
                        catch
                            eval([targetData, '=[];'])
                        end
                        
                        if eval(['~isfield(',targetData, ', ''subjectOrder'')'])
                            eval([targetData, '.subjectOrder=thisData;'])
                            eval([targetData, '.subjectIndices=', num2str(0),';'])
                        else
                            %Fill variables
                            eval([targetData, '.subjectOrder = [', targetData, '.subjectOrder "', thisData, '"];'])
                            eval([targetData, '.subjectIndices = [', targetData, '.subjectIndices ', num2str(0), '];'])
                        end
                    end
                end
            end
        end
    end
end