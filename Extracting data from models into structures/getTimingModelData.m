function data=getTimingModelData(path, modelFolders, modelNames, modelFileNames, mapNames, DTs);
%Get data from various candidate timing models

cd(path)
mrGlobals;
mrVista 3;
%load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsMix.mat')
allXvalDTs=DTs(2:3);
DTnames=["All", "Odd", "Even"];
path=char(path);

%Load all ROIs
VOLUME{1}=deleteAllROIs(VOLUME{1}); VOLUME{1}=refreshScreen(VOLUME{1},0);
okROIs=zeros(length(mapNames), 2);
for whichMap=1:length(mapNames)
    [VOLUME{1}, ok] = loadROI(VOLUME{1}, char(strcat('Left', mapNames(whichMap), 'Map')), [],[],[],1);
    if ok
        [VOLUME{1}] = loadROI(VOLUME{1}, char(strcat('Left', mapNames(whichMap), 'Low')), [],[],[],1);
        [VOLUME{1}] = loadROI(VOLUME{1}, char(strcat('Left', mapNames(whichMap), 'High')), [],[],[],1);
        okROIs(whichMap,1)=ok;
    else
        VOLUME{1} = makePointROI(VOLUME{1},[1 whichMap 1]);
        VOLUME{1} = makePointROI(VOLUME{1},[1 whichMap 2]);
        VOLUME{1} = makePointROI(VOLUME{1},[1 whichMap 3]);
    end
    
    [VOLUME{1}, ok] = loadROI(VOLUME{1}, char(strcat('Right', mapNames(whichMap), 'Map')), [],[],[],1);
    if ok
        [VOLUME{1}] = loadROI(VOLUME{1}, char(strcat('Right', mapNames(whichMap), 'Low')), [],[],[],1);
        [VOLUME{1}] = loadROI(VOLUME{1}, char(strcat('Right', mapNames(whichMap), 'High')), [],[],[],1);
        okROIs(whichMap,2)=ok;
    else
        VOLUME{1} = makePointROI(VOLUME{1},[1 whichMap 4]);
        VOLUME{1} = makePointROI(VOLUME{1},[1 whichMap 5]);
        VOLUME{1} = makePointROI(VOLUME{1},[1 whichMap 6]);
    end
end
VOLUME{1} = refreshScreen(VOLUME{1}, 0);


for whichDT=1:length(DTs)
    if DTs(whichDT)>0
        VOLUME{1}=viewSet(VOLUME{1}, 'curdt', DTs(whichDT));
        for whichFolder=1:length(modelFolders)
            folder=[path filesep 'Gray' filesep dataTYPES(DTs(whichDT)).name filesep char(modelFolders(whichFolder))];
            whichModelNames=modelNames{whichFolder};
            whichModelFiles=modelFileNames{whichFolder};
            for whichModel=1:length(whichModelNames)
                modelFile=dir([char(strcat(folder,filesep, '*', whichModelFiles(whichModel), '*'))]);
                
                %                     if isempty(modelFile)
                %                         VOLUME{1}=rmSelect(VOLUME{1}, 1,[], char(strcat('*', modelProg(whichProg), '2dOvalGaussian-', modelSpace(whichSpace), '*', whichModel)));
                %                     elseif length(modelFile)>1
                %                         VOLUME{1}=rmSelect(VOLUME{1}, 1,[], char(strcat('*', modelProg(whichProg), '2dOvalGaussian-', modelSpace(whichSpace), '*', whichModel)), length(modelFile));
                %                     else
                VOLUME{1}=rmSelect(VOLUME{1}, 1,[folder, filesep, modelFile.name]);
                %                     end
                %                     %cd(oldFolder);
                try
                VOLUME{1} = rmLoadDefault(VOLUME{1});
                catch
                    VOLUME{1}
                end
                VOLUME{1}=refreshScreen(VOLUME{1},0);
                
                if whichFolder~=2 && (whichDT==2 || whichDT==3)
                    xvalFolder=[path filesep 'Gray' filesep dataTYPES(DTs(5-whichDT)).name filesep char(modelFolders(whichFolder)) filesep 'xvalRefit'];
                    xvalModelFile=dir([char(strcat(xvalFolder,filesep, '*', whichModelFiles(whichModel), '*'))]);
                    xvalModel=load([xvalFolder, filesep, xvalModelFile.name], 'model');
                    xvalVes=rmGet(xvalModel.model{1}, 've');
                end
                
                for whichMap=1:length(mapNames)
                    dataName=char(strcat('data.', mapNames(whichMap), '.Left.', whichModelNames(whichModel), '.', DTnames(whichDT)));
                    
                    %Left hemisphere first entries
                    if okROIs(whichMap,1)
                        if whichDT==1 && whichFolder==1 && whichModel==1
                            dataTmp = RoiDistanceRatio(VOLUME{1}, 2+(whichMap-1)*6, 3+(whichMap-1)*6, 1+(whichMap-1)*6);
                        else
                            sourceDataName=char(strcat('data.', mapNames(whichMap), '.Left.', modelNames{1}(1), '.', DTnames(1)));
                            try
                            eval(['dataTmp=', sourceDataName, ';'])
                            catch
                                sourceDataName
                               % 
                            end
                            dataTmp.x0s=VOLUME{1}.rm.retinotopyModels{1}.x0(dataTmp.roiIndices);
                            dataTmp.y0s=VOLUME{1}.rm.retinotopyModels{1}.y0(dataTmp.roiIndices);
                            dataTmp.sigmas=VOLUME{1}.rm.retinotopyModels{1}.sigma.major(dataTmp.roiIndices);
                            dataTmp.sigmaMinor=VOLUME{1}.rm.retinotopyModels{1}.sigma.minor(dataTmp.roiIndices);
                            dataTmp.sigmaTheta=VOLUME{1}.rm.retinotopyModels{1}.sigma.theta(dataTmp.roiIndices);
                            if isfield(VOLUME{1}.rm.retinotopyModels{1}, 'exp')
                                dataTmp.exp=VOLUME{1}.rm.retinotopyModels{1}.exp(dataTmp.roiIndices);
                            end
                            
                            ves=1-(VOLUME{1}.rm.retinotopyModels{1}.rss(dataTmp.roiIndices)./VOLUME{1}.rm.retinotopyModels{1}.rawrss(dataTmp.roiIndices));
                            ves(~isfinite(ves)) = 0;
                            ves = max(ves, 0);
                            ves = min(ves, 1);
                            dataTmp.ves=ves;
                        end
                        dataTmp.folder=folder;
                        dataTmp.modelFile=modelFile.name;
                        if whichFolder~=2 && (whichDT==2 || whichDT==3)
                            dataTmp.vesXval=xvalVes(dataTmp.roiIndices);
                            dataTmp.folderXval=xvalFolder;
                            dataTmp.modelFileXval=xvalModelFile.name;
                        end
                        eval([dataName, '=dataTmp;'])
                        
                    end
                    
                    dataName=char(strcat('data.', mapNames(whichMap), '.Right.', whichModelNames(whichModel), '.', DTnames(whichDT)));
                    %Right hemisphere
                    if okROIs(whichMap,2)
                        if whichDT==1 && whichFolder==1 && whichModel==1
                            dataTmp = RoiDistanceRatio(VOLUME{1}, 5+(whichMap-1)*6, 6+(whichMap-1)*6, 4+(whichMap-1)*6);
                        else
                            sourceDataName=char(strcat('data.', mapNames(whichMap), '.Right.', modelNames{1}(1), '.', DTnames(1)));
                            eval(['dataTmp=', sourceDataName, ';'])
                            dataTmp.x0s=VOLUME{1}.rm.retinotopyModels{1}.x0(dataTmp.roiIndices);
                            dataTmp.y0s=VOLUME{1}.rm.retinotopyModels{1}.y0(dataTmp.roiIndices);
                            dataTmp.sigmas=VOLUME{1}.rm.retinotopyModels{1}.sigma.major(dataTmp.roiIndices);
                            dataTmp.sigmaMinor=VOLUME{1}.rm.retinotopyModels{1}.sigma.minor(dataTmp.roiIndices);
                            dataTmp.sigmaTheta=VOLUME{1}.rm.retinotopyModels{1}.sigma.theta(dataTmp.roiIndices);
                            if isfield(VOLUME{1}.rm.retinotopyModels{1}, 'exp')
                                dataTmp.exp=VOLUME{1}.rm.retinotopyModels{1}.exp(dataTmp.roiIndices);
                            end
                            
                            ves=1-(VOLUME{1}.rm.retinotopyModels{1}.rss(dataTmp.roiIndices)./VOLUME{1}.rm.retinotopyModels{1}.rawrss(dataTmp.roiIndices));
                            ves(~isfinite(ves)) = 0;
                            ves = max(ves, 0);
                            ves = min(ves, 1);
                            dataTmp.ves=ves;
                        end
                        dataTmp.folder=folder;
                        dataTmp.modelFile=modelFile.name;
                        if whichFolder~=2 && (whichDT==2 || whichDT==3)
                            dataTmp.vesXval=xvalVes(dataTmp.roiIndices);
                            dataTmp.folderXval=xvalFolder;
                            dataTmp.modelFileXval=xvalModelFile.name;
                        end
                        eval([dataName, '=dataTmp;'])
                    end
                    
                end
            end
        end
    end
end

close(1); mrvCleanWorkspace;
end

