%BMH, 01/2020

%First, group all scans into a single mrVista gray view data type
%Then reorder scan blocks, where all stimulus configurations are mixed in
%a scan run in pseudo-random order
RemixTimingScanBlocks;

%Create an ROI of the layer 1 voxels and save to local file
VOLUME{1}=makeGrayROI(VOLUME{1}); VOLUME{1}=refreshScreen(VOLUME{1},0);
VOLUME{1} = roiRestricttoLayer1(VOLUME{1},VOLUME{1}.selectedROI); 
VOLUME{1} = refreshScreen(VOLUME{1},0);
saveROI(VOLUME{1}, 'selected', 1);

%Make parameter definitions for models
MakeDurationVsPeriod;

%Define file paths for model params
%For all stimulus configurations mixed in a scan run
ParamsPathMix='/media/Storage3/WholeBrainNewSeg/TimingModelParamsMix.mat';
%For stimulus configurations in separate runs or sessions
ParamsPathSessions='/media/Storage3/WholeBrainNewSeg/TimingModelParamsSessions.mat';

%Define file paths for each subject's data
paths{1}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/S08_Combined';
paths{2}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/S09_Combined';
paths{3}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/S10_Combined';
paths{4}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/S11_Combined';
paths{5}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/S12_Combined';
paths{6}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/S13_Combined';
paths{7}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/HarvCombined';
paths{8}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/NellCombined';

%Define which mrVista dataTYPES contain your data for all, odd and even
%runs (odd and even for cross-validation)
DT=5:7;

%Set maximum number of cores to run models on simultaneously, typically the
%same as the number of dataTYPES
maxCores=length(DTs);

%Define which subject's data you want to model. 
%Numbers correspond to 'paths' list above.
whichSubs=[1:8];

%Run competing models for these subjects. Each will take a few days.
RunCandidateTimingModels(paths, whichSubs, DTs, maxCores, ParamsPathMix, ParamsPathSessions, 1);

%IF YOU DON'T WANT TO COMPARE CANDIDATE MODELS, JUST RUN THE BEST ONE
RunCandidateTimingModels(paths, whichSubs, DTs, maxCores, ParamsPathMix, ParamsPathSessions, 0);

%NOW MANUALLY MOVE THE FINAL MODELS YOU WANT TO COMPARE INTO SUB-FOLDERS,
%WITH THE SAME FOLDER NAME, WITHIN YOUR THREE DATATYPE FOLDERS

%Specify this folder's name
FinalModelFolder='FinalModels';

%Run cross validation
CrossValidateCandidateTimingModels(paths, whichSubs, DTs, FinalModelFolder, ParamsPathMix, ParamsPathSessions);

%%
%DEFINE ROIS BY SURFACE-BASED CLUSTERING. Manual ROI definition works fine, but this is
%less subjective and is faster when there are many, non-contiguous maps.
%Most of the steps are done in AFNI, so first we need to got our surfaces
%into a standardized nifti format.
%AFNI AND R NEED TO BE INSTALLED AND ACCESSIBLE FROM A BASH TERMINAL.
%To do this, we first need the WHOLE surface, while mrVista's gray view
%typically limits the surface to the nodes that overlap with the scan
%volume. We also need this limited to layer 1 nodes, touching the white
%matter border, so the starting segmentation is not suitable.
%SO FIRST MAKE A COPY OF ONE OF YOUR SUBJECT'S MRVISTA SESSIONS THAT WILL
%BE IRRETIEVABLY BROKEN. THEN GO INTO THIS SESSION FOLDER.
%Place the subject's T1 and segmentation into the folder.
T1name='t1_1mm.nii';
segName='t1_classFS_43.nii';
VistaSurfaceToNifti(T1name, segName, [], 4); %

%Now go back to the session folder in which your model was run, from 'paths'
%Load the ROI in which the model was run, as the first in the ROI list.
smoothLevel=5;
task=1;
% expression='x>0.09 & x<1.01 & y>0.09'; % Set to different values to make sure the 0.1 and 1 are actually included
expression='x<100 & y<100'; % Set to different values to make sure the 0.1 and 1 are actually included

% Model in question HAS to be loaded
VistaModelDataToNifti(VOLUME,T1name,expression,smoothLevel);  %Third input is 'task', 1 for timing, 2 for numerosity


%COPY THE RESULTING NIFTIS INTO THE PATH FOR YOUR SUBJECT FROM THE LIST
%BELOW

paths{1}='/media/Storage3/WholeBrainNewSeg/ForAFNIclustering/FromAlessioMatlab/S08_Data';
paths{2}='/media/Storage3/WholeBrainNewSeg/ForAFNIclustering/FromAlessioMatlab/S09_Data';
paths{3}='/media/Storage3/WholeBrainNewSeg/ForAFNIclustering/FromAlessioMatlab/S10_Data';
paths{4}='/media/Storage3/WholeBrainNewSeg/ForAFNIclustering/FromAlessioMatlab/S11ForAFNI';
paths{5}='/media/Storage3/WholeBrainNewSeg/ForAFNIclustering/FromAlessioMatlab/S12_Data';
paths{6}='/media/Storage3/WholeBrainNewSeg/ForAFNIclustering/FromAlessioMatlab/S13_Data';
paths{7}='/media/Storage3/WholeBrainNewSeg/ForAFNIclustering/FromAlessioMatlab/S01ForAFNI';
paths{8}='/media/Storage3/WholeBrainNewSeg/ForAFNIclustering/FromAlessioMatlab/S03ForAFNI';

% Why is the copying across necessary? Can't we just use the same directory
% to avoid even more directories and confusion?
% Just try I guess.
paths{9}='/home/benharvey/Inventory_WhatIsDone_WhatisOutdated/WhereWeCanBreakthem/BH_FullGray/mrVistaSession';
paths{10}='/home/benharvey/Inventory_WhatIsDone_WhatisOutdated/WhereWeCanBreakthem/YC_FullGray/mrVistaSession';
paths{11}='/home/benharvey/Inventory_WhatIsDone_WhatisOutdated/WhereWeCanBreakthem/SH_FullGray/mrVistaSession';
paths{12}='/home/benharvey/Inventory_WhatIsDone_WhatisOutdated/WhereWeCanBreakthem/JP_FullGray/mrVistaSession';
paths{13}='/home/benharvey/Inventory_WhatIsDone_WhatisOutdated/WhereWeCanBreakthem/MV_FullGray/mrVistaSession';
paths{14}='/home/benharvey/Inventory_WhatIsDone_WhatisOutdated/WhereWeCanBreakthem/MvA_FullGray/mrVistaSession';

whichSubs=9;
selSizes=56; %Cluster size (or set of sizes) in mm^2 of cortical surface
% selSizes=86; %Cluster size (or set of sizes) in mm^2 of cortical surface
selVarExps= 0.76:0.002:0.998;  %For variance explained percentiles (less useful, allows easier comparison across subjects)
selVarExps=0.01:0.01:0.40; %For actual variance explained
actualVE=1; %True if using actual variance explained, false if using percentiles.

% ALSO THE LOCATION FOR THE OUTPUT!
mainDir='/media/Storage3/WholeBrainNewSeg/ForAFNIclustering/FromAlessioMatlab'; %IMPORTANT: LOCATION OF R SCRIPTS 
AFNIclusterMultiThreshold(paths, whichSubs, task, smoothLevel, selSizes, selVarExps, actualVE, mainDir); %Surface-based clustering at a range of variance explained thresholds

%For each subject, choose the threshold where the cluster count flattens out from this plot
%of cluster count vs threshold

% For three test subjects, the resulting number of clusters is 1. Changed the
% expression values to be more inclusive, and got normal results back. Redo
% BH, YC and SH to check if this was the right solution. Check for BH,
% normal again. YC and SH still shit, MUST fix t1's 
subjectNum=9;
cd(mainDir);
AFNIclusterCountVsThreshold(selSizes, selVarExps, paths, subjectNum);

%Add this threshold here, either as quatile or actual variance explained
subjectNum=9;
selectedQuantile=0.9;
actualVE=1
selectedVarExp=0.2;
AFNIclustersToNifti(selSizes, selectedQuantile, selectedVarExp, actualVE);
 
%Convert one hemispheres cluster niftis to mrVista. There are too many
%clusters to do both hemispheres together
%Choose hemisphere
LeftOrRight='Left'; %Acceptable values are 'Left' and 'Right'
ClusterNifti=readFileNifti(strcat('Clusters', LeftOrRight, 'Layer1_zzz.nii.gz'));

%Now CD to the model session
mrVista 3;
VOLUME = NiftiClustersToVistaROIs(VOLUME, ClusterNifti,LeftOrRight)%, now automatically saves as left or right

%The resulting ROIs are numbered and should be named.
%Some may need to be joined or split into multiple maps.

%%
%Extract data from these models in the defined map ROIs
%These are the base names of the maps. Each should have an ROI for the map
%(called Left*Map), a line ROI for the low timing border (Left*Low), a line  
%ROI for the low timing border (Left*High), and a point ROI for the center
%(Left*Center)
mapNames=["TLO", "TTOP", "TTOA", "TPO", "TLS", "TPCI", "TPCM", "TPCS", "TFI" "TFS"];
%Folder names to get models from
modelFolders=["SearchFitFreeExponent", "HrfFitFreeExponent", "SearchFit", "Monotonic"];
%Names for models in output structure (for each folder)
modelNames{1}=["DurationPeriod", "OccupancyPeriod", "OnTimeOffTime", "LogDurationPeriod"];
modelNames{2}=["DurationPeriodHRF"];
modelNames{3}=["DurationPeriodNoExponent"];
modelNames{4}=["TemporalFreq", "Flat", "LinDF", "LinDCompF", "TemporalFreqPerEvent","GausOLinF", "GausDLinF", "GausOCompF", "GausDCompF",];
%Unique parts of corresponding model file names (for each folder)
modelFileNames{1}=["Lin-2dOvalGaussian-DurationPeriod", "Lin-2dOvalGaussian-OccupancyPeriod", "Lin-2dOvalGaussian-OnTimeOffTime", "Log-2dOvalGaussian-DurationPeriod"];
modelFileNames{2}=["Lin-2dOvalGaussian-DurationPeriod"];
modelFileNames{3}=["Lin-2dOvalGaussian-DurationPeriod"];
modelFileNames{4}=["TemporalFreq-", "Flat", "linearXlinearY", "linearXcompressiveY", "TemporalFreqPerEvent", "1dGaussianXlinearY-OccupancyFreq", "1dGaussianXlinearY-DurFreq", "1dGaussianXcompressiveY-OccupancyFreq", "1dGaussianXcompressiveY-DurFreq"];

%Names of each subject's resulting data structure
subjectOrder=["dataS8", "dataS9", "dataS10", "dataS11", "dataS12", "dataS13", "dataS1", "dataS3"];

%Which subjects' data to get
whichSubs=[1:8];

%Get that data. MOST SUBSEQUENT STEPS RELY ON THE DATA BEING STRUCTURED IN
%A SPECIFIC WAY, FROM THIS STEP.
for thisSub=whichSubs
    dataTmp=getTimingModelData(paths{thisSub}, modelFolders, modelNames, modelFileNames, mapNames, DT);
    eval([char(subjectOrder(thisSub)), '=dataTmp;'])
end

%Join all ROIs within subject into a field called 'All'
for thisSub=whichSubs
    thisData=char(subjectOrder(thisSub));
    UniteROIsRevision;
end

%Group ROI data across subject into a structure called 'dataAllSub'
for thisSub=whichSubs
    thisData=char(subjectOrder(thisSub));
    UniteROIsubsRevision;
end

%Names of models to be compared, in order of bars in resulting plot.
%Follows modelNames cells above, i.e. fields in subData structures.
modelNamesAll={'TemporalFreq', 'Flat', 'LinDF', 'LinDCompF', 'TemporalFreqPerEvent','GausOLinF', 'GausDLinF', 'GausOCompF', 'GausDCompF', 'DurationPeriodNoExponent', 'DurationPeriod', 'OccupancyPeriod', 'OnTimeOffTime', 'LogDurationPeriod'};

%Optionally, for comparisons that treat each subject (rather than each map)
%as an independent measure
mapNames=["All"];

%Now compare timig model fits
CompareTimingModelFits;

%Revert to full list of map names
mapNames=["TLO", "TTOP", "TTOA", "TPO", "TLS", "TPCI", "TPCM", "TPCS", "TFI" "TFS"];

%Histogram of timing preferences in each maps
TimingPreferenceHistograms;

%Make (very many) plots of timing vs distance
whichSubs=1:8;
whichPlots=[0 1 0 0 0 0 0];
TimingPlots;

%Make plots of tuning extent vs preferred timing, from group data
subjectOrder=[subjectOrder, "dataAllSubs"];
whichSubs=9;
whichPlots=[0 0 0 0 0 1 1];
TimingPlots;

%Get counts from resulting stats structure
countDurOddEven=zeros(size(stats{1}.pDurOddEven));
countPerOddEven=zeros(size(stats{1}.pDurOddEven));
countDistDur=zeros(size(stats{1}.pDistDurCorr));
countDistPer=zeros(size(stats{1}.pDistDurCorr));
countDurPer=zeros(size(stats{1}.pDurPerCorr));
for n=1:8
    countDurOddEven=countDurOddEven+(stats{n}.pDurOddEven<0.05 & stats{n}.rDurOddEven>0);
    countPerOddEven=countPerOddEven+(stats{n}.pPerOddEven<0.05 & stats{n}.rPerOddEven>0);
    
    countDistDur=countDistDur+(stats{n}.pDistDurCorr<0.05 & stats{n}.rDistDurCorr>0);
    countDistPer=countDistPer+(stats{n}.pDistPerCorr<0.05 & stats{n}.rDistPerCorr>0);
    
    countDurPer=countDurPer+(stats{n}.pDurPerCorr<0.05 & stats{n}.rDurPerCorr>0);
end

%ANOVAs of changes through hierarchy
TimingANOVAs(dataS8, dataS9, dataS10, dataS11, dataS12, dataS13, dataS1, dataS3);
































%Group MNI coords
subOrder={'S8', 'S9', 'S10', 'S11', 'S12', 'S13', 'S14', 'S15'};
mapNames={'TLO', 'TTOP', 'TTOA', 'TPO', 'TLS', 'TPCI', 'TPCM', 'TPCS', 'TFI' 'TFS'};%Revision
hemiNames={'Left', 'Right'};
okROI=zeros([length(subOrder), length(mapNames), length(hemiNames)]);
coordsROI=nan([length(subOrder), length(mapNames), length(hemiNames), 3]);
for sub=1:length(subOrder)
    for map=1:length(mapNames)
        for hemi=1:length(hemiNames)
            for n=1:eval(['length(names', subOrder{sub}, ')'])
                if strcmp(eval(['names', subOrder{sub}, '{', num2str(n), '}']), [hemiNames{hemi}, mapNames{map}, 'Center'])
                    okROI(sub, map, hemi)=1;
                    coordsROI(sub, map, hemi, :)=eval(['points', subOrder{sub}, '(', num2str(n), ',:)']);
                end
            end
        end
    end
end

%Populate table, alternating means and SDs
[round(nanmean(coordsROI(:,:,1,1),1)); round(nanstd(coordsROI(:,:,1,1),1)); round(nanmean(coordsROI(:,:,1,2),1)); round(nanstd(coordsROI(:,:,1,2),1)); round(nanmean(coordsROI(:,:,1,3),1)); round(nanstd(coordsROI(:,:,1,3),1));]'
[round(nanmean(coordsROI(:,:,2,1),1)); round(nanstd(coordsROI(:,:,2,1),1)); round(nanmean(coordsROI(:,:,2,2),1)); round(nanstd(coordsROI(:,:,2,2),1)); round(nanmean(coordsROI(:,:,2,3),1)); round(nanstd(coordsROI(:,:,2,3),1));]'

%AFNI results
origin=[128 145 140];
origin8=repmat(origin, [8,1])

coordsROIAFNI=nan([length(subOrder), length(mapNames), length(hemiNames), 3]);
okROIAFNI=zeros([length(subOrder), length(mapNames), length(hemiNames)]);
for map=1:length(mapNames)
    for hemi=1:length(hemiNames)
        eval(['thisData=CenterCoords', hemiNames{hemi}, mapNames{map},';']);
        %Cutting missing data
        for line=1:8
            if thisData(line,1)==85 && thisData(line,2)==95 && thisData(line,3)==83
                thisData(line,:)=nan;
            else
                okROIAFNI(line, map, hemi)=1;
            end
        end            
        tmp=(thisData-origin8);
        MNI=[-tmp(:,1) tmp(:,3) -tmp(:,2)];
        coordsROIAFNI(:, map, hemi,:)=MNI;
    end
end

%Populate table, alternating means and SDs
[round(nanmean(coordsROIAFNI(:,:,1,1),1)); round(nanstd(coordsROIAFNI(:,:,1,1),1)); round(nanmean(coordsROIAFNI(:,:,1,2),1)); round(nanstd(coordsROIAFNI(:,:,1,2),1)); round(nanmean(coordsROIAFNI(:,:,1,3),1)); round(nanstd(coordsROIAFNI(:,:,1,3),1));]'
[round(nanmean(coordsROIAFNI(:,:,2,1),1)); round(nanstd(coordsROIAFNI(:,:,2,1),1)); round(nanmean(coordsROIAFNI(:,:,2,2),1)); round(nanstd(coordsROIAFNI(:,:,2,2),1)); round(nanmean(coordsROIAFNI(:,:,2,3),1)); round(nanstd(coordsROIAFNI(:,:,2,3),1));]'



