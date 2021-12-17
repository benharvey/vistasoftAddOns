
% modelProgName=["Log", "Lin"];%, "LinLog"];
% modelSpace=["DurationPeriod", "OccupancyPeriod", "OnTimeOffTime"];
DTnames=["All", "Odd", "Even"];
DTs=5:7;
% 
% shapes=["Oval", "Circ", "OvalFreeExponent", "OvalFreeExponentHRF"];
% 
% mapNames=["TLO", "TTOP", "TTOA", "TPO", "TLS", "TPCI", "TPCM", "TPCS", "TFI" "TFS"];
% %Folder names to ge models from
% modelFolders=["SearchFitFreeExponent", "HrfFitFreeExponent", "SearchFit", "Monotonic"];
% %Names for models in output structure (for each folder)
% modelNames{1}=["DurationPeriod", "OccupancyPeriod", "OnTimeOffTime", "LogDurationPeriod"];
% modelNames{2}=["DurationPeriodHRF"];
% modelNames{3}=["DurationPeriodNoExponent"];
% modelNames{4}=["TemporalFreq", "Flat", "LinDF", "LinDCompF", "TemporalFreqPerEvent","GausOLinF", "GausDLinF", "GausOCompF", "GausDCompF",];
% %Unique parts of corresponding model file names (for each folder)
% modelFileNames{1}=["Lin-2dOvalGaussian-DurationPeriod", "Lin-2dOvalGaussian-OccupancyPeriod", "Lin-2dOvalGaussian-OnTimeOffTime", "Log-2dOvalGaussian-DurationPeriod"];
% modelFileNames{2}=["Lin-2dOvalGaussian-DurationPeriod"];
% modelFileNames{3}=["Lin-2dOvalGaussian-DurationPeriod"];
% modelFileNames{4}=["TemporalFreq-", "Flat", "linearXlinearY", "linearXcompressiveY", "TemporalFreqPerEvent", "1dGaussianXlinearY-OccupancyFreq", "1dGaussianXlinearY-DurFreq", "1dGaussianXcompressiveY-OccupancyFreq", "1dGaussianXcompressiveY-DurFreq"];
% 
% paths{1}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/S08_Combined';
% paths{2}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/S09_Combined';
% paths{3}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/S10_Combined';
% paths{4}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/S11_Combined';
% paths{5}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/S12_Combined';
% paths{6}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/S13_Combined_Seg20';
% paths{7}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/HarvCombined';
% paths{8}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/NellCombined';
% subData=["dataS8", "dataS9", "dataS10", "dataS11", "dataS12", "dataS13", "dataHarv", "dataNell"];


hemispheres=["Left", "Right"];
for whichDT=1:length(DTs)
    if DTs(whichDT)>0
        for whichFolder=1:length(modelFolders)
            %folder=[path filesep 'Gray' filesep dataTYPES(DTs(whichDT)).name filesep char(modelFolders(whichFolder))];
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
                                
                                if whichFolder~=2 && whichDT>1 
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
