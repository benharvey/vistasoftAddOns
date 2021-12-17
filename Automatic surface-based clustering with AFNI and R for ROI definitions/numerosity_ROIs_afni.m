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

task = 2; % task: 1 = timing, 2 = numerosity

ROIpaths{1} = '/mnt/WholeBrainNewSeg-afni/numerosity_ROIs_afni/Dumo_Combined_smooth';modelpaths{1} = '/mnt/WholeBrainNewSeg-afni/numerosity_ROIs_afni/Dumo_Combined_orig';
ROIpaths{2} = '/mnt/WholeBrainNewSeg-afni/numerosity_ROIs_afni/Frac_Combined_smooth';modelpaths{2} = '/mnt/WholeBrainNewSeg-afni/numerosity_ROIs_afni/Frac_Combined_orig';
ROIpaths{3} = '/mnt/WholeBrainNewSeg-afni/numerosity_ROIs_afni/Harv_Combined_smooth';modelpaths{3} = '/mnt/WholeBrainNewSeg-afni/numerosity_ROIs_afni/Harv_Combined_orig';
ROIpaths{4} = '/mnt/WholeBrainNewSeg-afni/numerosity_ROIs_afni/Klei_Combined_smooth';modelpaths{4} = '/mnt/WholeBrainNewSeg-afni/numerosity_ROIs_afni/Klei_Combined_orig';
ROIpaths{5} = '/mnt/WholeBrainNewSeg-afni/numerosity_ROIs_afni/Nell_Combined_smooth';modelpaths{5} = '/mnt/WholeBrainNewSeg-afni/numerosity_ROIs_afni/Nell_Combined_orig';
ROIpaths{6} = '/mnt/WholeBrainNewSeg-afni/numerosity_ROIs_afni/S8_Combined_smooth';modelpaths{6} = '/mnt/WholeBrainNewSeg-afni/numerosity_ROIs_afni/S8_Combined_orig';
ROIpaths{7} = '/mnt/WholeBrainNewSeg-afni/numerosity_ROIs_afni/S9_Combined_smooth';modelpaths{7} = '/mnt/WholeBrainNewSeg-afni/numerosity_ROIs_afni/S9_Combined_orig';
ROIpaths{8} = '/mnt/WholeBrainNewSeg-afni/numerosity_ROIs_afni/S10_Combined_smooth';modelpaths{8} = '/mnt/WholeBrainNewSeg-afni/numerosity_ROIs_afni/S10_Combined_orig';
ROIpaths{9} = '/mnt/WholeBrainNewSeg-afni/numerosity_ROIs_afni/S11_Combined_smooth';modelpaths{9} = '/mnt/WholeBrainNewSeg-afni/numerosity_ROIs_afni/S11_Combined_orig';
ROIpaths{10} = '/mnt/WholeBrainNewSeg-afni/numerosity_ROIs_afni/S12_Combined_smooth';modelpaths{10} = '/mnt/WholeBrainNewSeg-afni/numerosity_ROIs_afni/S12_Combined_orig';
ROIpaths{11} = '/mnt/WholeBrainNewSeg-afni/numerosity_ROIs_afni/S13_Combined_smooth';modelpaths{11} = '/mnt/WholeBrainNewSeg-afni/numerosity_ROIs_afni/S13_Combined_orig';

clusterpaths{1} = '/mnt/WholeBrainNewSeg-afni/Dumo_Combined';
clusterpaths{2} = '/mnt/WholeBrainNewSeg-afni/Frac_Combined';
clusterpaths{3} = '/mnt/WholeBrainNewSeg-afni/Harv_Combined';
clusterpaths{4} = '/mnt/WholeBrainNewSeg-afni/Klei_Combined';
clusterpaths{5} = '/mnt/WholeBrainNewSeg-afni/Nell_Combined';
clusterpaths{6} = '/mnt/WholeBrainNewSeg-afni/S8_Combined';
clusterpaths{7} = '/mnt/WholeBrainNewSeg-afni/S9_Combined';
clusterpaths{8} = '/mnt/WholeBrainNewSeg-afni/S10_Combined';
clusterpaths{9} = '/mnt/WholeBrainNewSeg-afni/S11_Combined';
clusterpaths{10} = '/mnt/WholeBrainNewSeg-afni/S12_Combined';
clusterpaths{11} = '/mnt/WholeBrainNewSeg-afni/S13_Combined';

whichSubs = 1;
thisSub = 1;
cd(ROIpaths{whichSubs}(1:end-7))

T1_path = dir('t1_1mm.nii.gz');
T1name = T1_path.name;
T1_path = [T1_path.folder,'/',T1_path.name];

seg_path = dir('*1mm_*.nii.gz');if isempty(seg_path),seg_path = dir('*class*.nii.gz');end
segName = seg_path.name;
seg_path = {[seg_path.folder,'/',seg_path.name]};

numGrayLayers = 4;

startup vn

VOLUME{1} = VistaSurfaceToNifti(T1name, segName, seg_path, numGrayLayers);

close(1);close(2);mrvCleanWorkspace

system(['cp ',ROIpaths{whichSubs},'/allTmp.mat ',modelpaths{whichSubs}]);
system(['cp ',ROIpaths{whichSubs},'/SegmentationLayer1.nii.gz ',modelpaths{whichSubs}]);

cd(modelpaths{whichSubs})
mrVista 3

modelName = {'*Log-FullBlanks-DT0.5-Fineret*gFit-gFit*'};
load('/mnt/WholeBrainNewSeg-afni/allModelParamsMix.mat');
stim = NumerosityDTs{whichSubs(thisSub)}(1); % NumerosityAll (L1)
stimName = regexprep((dataTYPES(1,stim).name), '[^\w'']','');
stimFolder = fullfile(modelpaths{whichSubs},'Gray',dataTYPES(1,stim).name);
cd(stimFolder)
modelFile = [stimFolder,'/',strcat(ls(modelName{1}))];
cd ../..

VOLUME{1}.curDataType = NumerosityDTs{whichSubs(thisSub)}(1); % NumerosityAll (L1)
VOLUME{1} = refreshScreen(VOLUME{1},0);
VOLUME{1} = rmSelect(VOLUME{1}, 2, modelFile); VOLUME{1} = rmLoadDefault(VOLUME{1});
ROIFolder = fullfile(modelpaths{whichSubs},'Gray/ROIs/gray-Layer1.mat');
VOLUME{1} = loadROI(VOLUME{1}, ROIFolder, [], [], 1, 1); VOLUME{1} = selectCurROISlice(VOLUME{1}); VOLUME{1} = refreshScreen(VOLUME{1});

%Now go back to the session folder in which your model was run, from 'paths'
%Load the ROI in which the model was run, as the first in the ROI list.
smoothLevel = 10;
VistaModelDataToNifti(VOLUME, T1name, task, smoothLevel);

selSizes = 86; %Cluster size (or set of sizes) in mm^2 of cortical surface
selVarExps =  0.76:0.002:0.998;  %#ok<NASGU> %For variance explained percentiles (less useful, allows easier comparison across subjects)
selVarExps = 0.01:0.01:0.40; %For actual variance explained
actualVE = 1; %True if using actual variance explained, false if using percentiles.
mainDir = '/mnt/WholeBrainNewSeg-afni/matlab/Clustering_for_ROI_definitions'; %IMPORTANT: LOCATION OF R SCRIPTS 
AFNIclusterMultiThreshold(ROIpaths, whichSubs, task, smoothLevel, selSizes, selVarExps, actualVE, mainDir); %Surface-based clusterig at a range of variance explained thresholds

%For each subject, choose the threshold where the cluster count flattens out from this plot
%of cluster count vs threshold
subjectNum = 6;
cd(ROIpaths{subjectNum});
AFNIclusterCountVsThreshold(selSizes, selVarExps, ROIpaths, whichSubs);

%Add this threshold here, either as quatile or actual variance explained
subjectNum = 6;
selectedQuantile = 0.90;
actualVE = 1;
selectedVarExp = 0.2;
AFNIclustersToNifti; %Use that clustering threshold to get clusters as nifti

startup vn

%Convert one hemispheres cluster niftis to mrVista. There are too many
%clusters to do both hemispheres together
%Choose hemisphere
LeftOrRight = 'Left'; %Acceptable values are 'Left' and 'Right'
ClusterNifti = readFileNifti(strcat('Clusters', LeftOrRight, 'Layer1_zzz.nii.gz'));
cd(clusterpaths{subjectNum})
mrVista 3;
VOLUME = NiftiClustersToVistaROIs(VOLUME, ClusterNifti);
close all
mrvCleanWorkspace;

cd(ROIpaths{subjectNum})
LeftOrRight = 'Right'; %Acceptable values are 'Left' and 'Right'
ClusterNifti = readFileNifti(strcat('Clusters', LeftOrRight, 'Layer1_zzz.nii.gz'));
cd(clusterpaths{subjectNum})
mrVista 3;
VOLUME = NiftiClustersToVistaROIs(VOLUME, ClusterNifti);
close all
mrvCleanWorkspace;

%The resulting ROIs are numbered and should be named.
%Some may need to be joined or split into multiple maps.