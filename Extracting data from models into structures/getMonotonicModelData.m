function data=getMonotonicModelData(path, modelFolders, modelNames, modelFileNames, mapNames, DTsAll, allXvalDTs, DTnames);
%Get data from various candidate timing models

cd(path)
mrGlobals;
mrVista 3;
%load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsMix.mat')

DTnames=["Area","AreaOdd", "AreaEven", "Size","SizeOdd","SizeEven",  "Circ","CircOdd", "CircEven","Dense","DenseOdd","DenseEven"];
path=char(path);

%Load all ROIs
VOLUME{1}=deleteAllROIs(VOLUME{1}); VOLUME{1}=refreshScreen(VOLUME{1},0);
okROIs=zeros(length(mapNames), 2);
for whichMap=1:length(mapNames)
    [VOLUME{1}, ok] = loadROI(VOLUME{1}, char(strcat('Left', mapNames(whichMap))), [],[],[],1);
    [VOLUME{1}, ok] = loadROI(VOLUME{1}, char(strcat('Right', mapNames(whichMap))), [],[],[],1);
end
VOLUME{1} = refreshScreen(VOLUME{1}, 0);


for whichDT=1:length(DTsAll)
    if DTsAll(whichDT)>0
        VOLUME{1}=viewSet(VOLUME{1}, 'curdt', DTsAll(whichDT));
        for whichFolder=1:length(modelFolders)
            folder=[path filesep 'Gray' filesep dataTYPES(DTsAll(whichDT)).name filesep char(modelFolders(whichFolder))];
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
                
                if ismember(DTsAll(whichDT), allXvalDTs)
                    %NOTE: Looks in the other DT of the cross-validatin
                    %pair
                    modul=mod(whichDT, 3);
                    if modul==0;
                        modul=3;
                    end
                    otherDT=floor((whichDT-1)/3)*3+5-modul
                    xvalFolder=[path filesep 'Gray' filesep dataTYPES(otherDT).name filesep char(modelFolders(whichFolder)) filesep 'xvalRefit'];
                    xvalModelFile=dir([char(strcat(xvalFolder,filesep, '*', whichModelFiles(whichModel), '*'))]);
                    xvalModel=load([xvalFolder, filesep, xvalModelFile.name], 'model');
                    xvalVes=rmGet(xvalModel.model{1}, 've');
                    xvalRSS=rmGet(xvalModel.model{1}, 'rss');
                    xvalRawRSS=rmGet(xvalModel.model{1}, 'rawrss');
                end
                
                for whichMap=1:length(mapNames)
                    dataName=char(strcat('data.', mapNames(whichMap), '.Left.', whichModelNames(whichModel), '.', DTnames(whichDT)));
                    
                    %Left hemisphere first entries
                    [tmp, dataTmp.roiIndices]=intersectCols(VOLUME{1}.coords, VOLUME{1}.ROIs((whichMap-1)*2+1).coords);
                    dataTmp.x0s=VOLUME{1}.rm.retinotopyModels{1}.x0(dataTmp.roiIndices);
                    dataTmp.y0s=VOLUME{1}.rm.retinotopyModels{1}.y0(dataTmp.roiIndices);
                    dataTmp.sigmas=VOLUME{1}.rm.retinotopyModels{1}.sigma.major(dataTmp.roiIndices);
                    dataTmp.sigmaMinor=VOLUME{1}.rm.retinotopyModels{1}.sigma.minor(dataTmp.roiIndices);
                    dataTmp.sigmaTheta=VOLUME{1}.rm.retinotopyModels{1}.sigma.theta(dataTmp.roiIndices);
                    if isfield(VOLUME{1}.rm.retinotopyModels{1}, 'exp')
                        dataTmp.exp=VOLUME{1}.rm.retinotopyModels{1}.exp(dataTmp.roiIndices);
                    end
                    
                    dataTmp.rss=VOLUME{1}.rm.retinotopyModels{1}.rss(dataTmp.roiIndices);
                    dataTmp.rawrss=VOLUME{1}.rm.retinotopyModels{1}.rawrss(dataTmp.roiIndices);
                    
                    %Converts rss to ve
                    ves=1-(dataTmp.rss./dataTmp.rawrss;
                    ves(~isfinite(ves)) = 0;
                    ves = max(ves, 0);
                    ves = min(ves, 1);
                    dataTmp.ves=ves;
                    
                    dataTmp.folder=folder;
                    dataTmp.modelFile=modelFile.name;
                    if ismember(DTsAll(whichDT), allXvalDTs)
                        dataTmp.vesXval=xvalVes(dataTmp.roiIndices);
                        dataTmp.rssXval=xvalRSS(dataTmp.roiIndices);
                        dataTmp.rawrssXval=xvalRawRSS(dataTmp.roiIndices);
                        dataTmp.folderXval=xvalFolder;
                        dataTmp.modelFileXval=xvalModelFile.name;
                    end
                    eval([dataName, '=dataTmp;'])
                    
                    
                    dataName=char(strcat('data.', mapNames(whichMap), '.Right.', whichModelNames(whichModel), '.', DTnames(whichDT)));
                    
                    %Left hemisphere first entries
                    [tmp, dataTmp.roiIndices]=intersectCols(VOLUME{1}.coords, VOLUME{1}.ROIs((whichMap-1)*2+2).coords);
                    dataTmp.x0s=VOLUME{1}.rm.retinotopyModels{1}.x0(dataTmp.roiIndices);
                    dataTmp.y0s=VOLUME{1}.rm.retinotopyModels{1}.y0(dataTmp.roiIndices);
                    dataTmp.sigmas=VOLUME{1}.rm.retinotopyModels{1}.sigma.major(dataTmp.roiIndices);
                    dataTmp.sigmaMinor=VOLUME{1}.rm.retinotopyModels{1}.sigma.minor(dataTmp.roiIndices);
                    dataTmp.sigmaTheta=VOLUME{1}.rm.retinotopyModels{1}.sigma.theta(dataTmp.roiIndices);
                    if isfield(VOLUME{1}.rm.retinotopyModels{1}, 'exp')
                        dataTmp.exp=VOLUME{1}.rm.retinotopyModels{1}.exp(dataTmp.roiIndices);
                    end
                    
                    dataTmp.rss=VOLUME{1}.rm.retinotopyModels{1}.rss(dataTmp.roiIndices);
                    dataTmp.rawrss=VOLUME{1}.rm.retinotopyModels{1}.rawrss(dataTmp.roiIndices);
                    
                    %Converts rss to ve
                    ves=1-(dataTmp.rss./dataTmp.rawrss;
                    ves(~isfinite(ves)) = 0;
                    ves = max(ves, 0);
                    ves = min(ves, 1);
                    dataTmp.ves=ves;
                    
                    dataTmp.folder=folder;
                    dataTmp.modelFile=modelFile.name;
                    if ismember(DTsAll(whichDT), allXvalDTs)
                        dataTmp.vesXval=xvalVes(dataTmp.roiIndices);
                        dataTmp.rssXval=xvalRSS(dataTmp.roiIndices);
                        dataTmp.rawrssXval=xvalRawRSS(dataTmp.roiIndices);
                        dataTmp.folderXval=xvalFolder;
                        dataTmp.modelFileXval=xvalModelFile.name;
                    end
                    eval([dataName, '=dataTmp;'])
                    
                end
            end
        end
    end
end

close(1); mrvCleanWorkspace;
end

