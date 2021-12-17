startup vo
mrVista
mrVista 3
load('/media/Storage2/WholeBrainNewSeg/Nell/Nell030713Area/Inplane/anat.mat')
INPLANE{1}
INPLANE{1}.ROIs
INPLANE{1}.ROIs.coords(:,1:10)
for n=1:size(INPLANE{1}.ROIs.coords,2); anat(INPLANE{1}.ROIs.coords(1,n), INPLANE{1}.ROIs.coords(2,n), INPLANE{1}.ROIs.coords(3,n))=0; end
figure; imagesc(anat(:,:,21)
figure; imagesc(anat(:,:,21))
load('/media/Storage2/WholeBrainNewSeg/Nell/Nell030713Area/Inplane/anat.mat')
anat1=anat([5:128 1:4],:,:);
anat2=anat([125:128 1:124],:,:);
anat3=anat(:,[5:128 1:4],:);
anat4=anat(:,[125:128 1:124],:);
anat5=anat(:,:,[40:43 1:39]);
anat6=anat(:,:,[5:43 1:4]);
anatall=min(anat1,anat2);
anatall=min(anatall,anat3);
anatall=min(anatall,anat4);
anatall=min(anatall,anat5);
anatall=min(anatall,anat6);
figure; imagesc(anatall(:,:,21))
for n=1:size(INPLANE{1}.ROIs.coords,2); anatall2(INPLANE{1}.ROIs.coords(1,n), INPLANE{1}.ROIs.coords(2,n), INPLANE{1}.ROIs.coords(3,n))=0; end
figure; imagesc(anatall2(:,:,21))
anat1=anat([7:128 1:6],:,:);
anat2=anat([123:128 1:122],:,:);
anat3=anat(:,[7:128 1:6],:);
anat4=anat(:,[123:128 1:122],:);
anat5=anat(:,:,[38:43 1:37]);
anat6=anat(:,:,[7:43 1:6]);
anatall=min(anat1,anat2);
anatall=min(anatall,anat3);
anatall=min(anatall,anat4);
anatall=min(anatall,anat5);
anatall=min(anatall,anat6);
for n=1:size(INPLANE{1}.ROIs.coords,2); anatall2(INPLANE{1}.ROIs.coords(1,n), INPLANE{1}.ROIs.coords(2,n), INPLANE{1}.ROIs.coords(3,n))=0; end
figure; imagesc(anatall2(:,:,21))
anatall2(58:72,:,:)=0;
figure; imagesc(anatall2(:,:,21))
for n=1:size(INPLANE{1}.ROIs.coords,2); anatall2(INPLANE{1}.ROIs.coords(1,n), INPLANE{1}.ROIs.coords(2,n), INPLANE{1}.ROIs.coords(3,n))=0; end
anatall2(:,58:72,:)=0;
figure; imagesc(anatall2(:,:,21))
for n=1:size(INPLANE{1}.ROIs.coords,2); anatall2(INPLANE{1}.ROIs.coords(1,n), INPLANE{1}.ROIs.coords(2,n), INPLANE{1}.ROIs.coords(3,n))=0; end
anatall2(:,58:72,:)=0;
figure; imagesc(anatall2(:,:,21))
anatall2(:,58:72,:)=0;
figure; imagesc(anatall2(:,:,10))
mask=anatall2>300;
sum(mask(:))
1288128
128*128

%For each subject
fname='T1_1mm_class_20';

thinSegmentation([fname, '.nii.gz'], 1);
thinSegmentation([fname, 'thinned.nii.gz'], 1);
thinSegmentation([fname, 'thinnedthinned.nii.gz'], 1);
segmentation=readFileNifti([fname, 'thinnedthinnedthinned.nii.gz']);
[xs, ys, zs]=ind2sub(size(segmentation.data), find(segmentation.data==3 | segmentation.data==4));

%For each acquired session, in order
mrVista;
mrVista 3;
xform = inv(mrSESSION.alignment);
ipVoxSize = mrSESSION.inplanes.voxelSize;
volVoxSize = readVolAnatHeader(vANATOMYPATH);

VOLUME{1}=makeGrayROI(VOLUME{1}); VOLUME{1}=refreshScreen(VOLUME{1},0);
 volume=checkSelectedVolume; INPLANE{1}=vol2ipAllROIs(volume,INPLANE{1}); INPLANE{1}=refreshScreen(INPLANE{1},0); clear volume
 INPLANE{1}=refreshScreen(INPLANE{1},0);
 
coords = xformROIcoords([size(VOLUME{1}.anat, 1)-zs, size(VOLUME{1}.anat, 2)-ys, xs]', xform, volVoxSize, ipVoxSize);

% Toss coords outside the inplanes
ipSize = viewSize(INPLANE{1});
indices = ((coords(1,:) >= 1) & (coords(1,:) <= ipSize(1)) & ...
    (coords(2,:) >= 1) & (coords(2,:) <= ipSize(2)) & ...
    (coords(3,:) >= 1) & (coords(3,:) <= ipSize(3)));
coords = coords(:,indices);
INPLANE{1}.ROIs(1).coords=coords;
 INPLANE{1}=refreshScreen(INPLANE{1},0);

load(['Inplane/MotionComp/TSeries/Scan', num2str(1),'/tSeries', num2str(1), '.mat']);
nScans=size(dir('Inplane/MotionComp/TSeries'),1)-2;
nFrames=size(tSeries,1);
tSeriesall=zeros([size(INPLANE{1}.anat,1), size(INPLANE{1}.anat,2), size(INPLANE{1}.anat,3), nFrames, nScans]); 

for scanCounter=1:nScans
    for n=1:size(INPLANE{1}.anat,3);
        load(['Inplane/MotionComp/TSeries/Scan', num2str(scanCounter),'/tSeries', num2str(n), '.mat']);
        tSeries=reshape(tSeries, [nFrames size(INPLANE{1}.anat,1) size(INPLANE{1}.anat,2)]);
        tSeries=permute(tSeries, [2 3 1]);
        tSeriesall(:,:,n,:, scanCounter)=tSeries;
    end
end

tSeriesOut=zeros([nFrames, size(coords, 2), nScans]);
for counter=1:size(coords, 2)
    tSeriesOut(:,counter, :)=tSeriesall(coords(1,counter),coords(2,counter),coords(3, counter),:,:);
end

%After all 4 sessions are saved as tSeriesOut1, 2, 3 and 4, join session
%together
minLength=min([size(tSeriesOut1,2) size(tSeriesOut2,2) size(tSeriesOut3,2) size(tSeriesOut4,2)]);
tSeriesOutAll=cat(3, tSeriesOut1(:,1:minLength, :), tSeriesOut2(:,1:minLength, :), tSeriesOut3(:,1:minLength, :), tSeriesOut4(:,1:minLength, :));

% minLength=min([size(tSeriesOut1,2) size(tSeriesOut2,2) size(tSeriesOut3,2) size(tSeriesOut4,2) size(tSeriesOut5,2)]);
% tSeriesOutAll=cat(3, tSeriesOut1(:,1:minLength, :), tSeriesOut2(:,1:minLength, :), tSeriesOut3(:,1:minLength, :), tSeriesOut4(:,1:minLength, :), tSeriesOut5(:,1:minLength, :));


%Rearrange scans for timing data
S8order=[22 10 6 16 17 3, 9 20 14 23 13 11 15, 4 19 24 2 12 1 7, 21 18 5 8 1]; %Added 1 (1234) for extra scan
S9order=[7 16 10 11 17 20 6, 15 8 21 5 9 1 2, 19 22 12 18 24 23 14, 13 3 4 1 24 11 15];
S10order=[23 9 24 1 13 3, 2 4 11 15 12 10, 8 7 18 19 16, 17 14 6 20 5 21 22];
S11order=[8 16 10 11 6 21, 9 14 18 17 23, 22 24 5 1 13 7, 3 2 19 20 12 15 8]; %Repeated number 8 due to excessive motion. Dropped scan 4 (number 4 in session 2) from bad data
S12order=[22 6 3 16 11 7, 17 14 8 5 21 19 15, 1 23 2 4 18, 24 13 9 20 10 12 1];

S13order=[7 20 21 6 12 24 14, 22 3 10 8 15 19 11, 5 4 16 1 13 2 23 17, 18 9];

ThisSubjectOrder=S13order;


%Convert scan order to stimulus sweep
Scan2Sweep(1,:)=[1 2 3 4];
Scan2Sweep(2,:)=[1 3 2 4];
Scan2Sweep(3,:)=[2 1 3 4];
Scan2Sweep(4,:)=[2 3 1 4];
Scan2Sweep(5,:)=[3 1 2 4];
Scan2Sweep(6,:)=[3 2 1 4];
Scan2Sweep(7,:)=[1 2 4 3];
Scan2Sweep(8,:)=[1 3 4 2];
Scan2Sweep(9,:)=[2 1 4 3];
Scan2Sweep(10,:)=[2 3 4 1];
Scan2Sweep(11,:)=[3 1 4 2];
Scan2Sweep(12,:)=[3 2 4 1];
Scan2Sweep(13,:)=[1 4 2 3];
Scan2Sweep(14,:)=[1 4 3 2];
Scan2Sweep(15,:)=[2 4 1 3];
Scan2Sweep(16,:)=[2 4 3 1];
Scan2Sweep(17,:)=[3 4 1 2];
Scan2Sweep(18,:)=[3 4 2 1];
Scan2Sweep(19,:)=[4 1 2 3];
Scan2Sweep(20,:)=[4 1 3 2];
Scan2Sweep(21,:)=[4 2 1 3];
Scan2Sweep(22,:)=[4 2 3 1];
Scan2Sweep(23,:)=[4 3 1 2];
Scan2Sweep(24,:)=[4 3 2 1];

%Initialize scan averaging structure
Sweeps{1}.data=zeros([56, size(tSeriesOutAll,2), size(tSeriesOutAll,3)]);
Sweeps{2}.data=zeros([56, size(tSeriesOutAll,2), size(tSeriesOutAll,3)]);
Sweeps{3}.data=zeros([56, size(tSeriesOutAll,2), size(tSeriesOutAll,3)]);
Sweeps{4}.data=zeros([56, size(tSeriesOutAll,2), size(tSeriesOutAll,3)]);

%For each scan
for thisScan=1:size(tSeriesOutAll,3)
   
    %Extract the 4 sweeps
    tSeriesSweeps= reshape(tSeriesOutAll(:,:,thisScan), [56, 4, size(tSeriesOutAll, 2)]);
    for thisSweep=1:4
        Sweeps{Scan2Sweep(ThisSubjectOrder(thisScan),thisSweep)}.data(:, :, thisScan)=squeeze(tSeriesSweeps(:,thisSweep,:));
    end
end

%Create new dataTYPES in combined session(if FIRST SUBJECT)
% duplicateDataType(VOLUME{1}, 'TimingSweepsNull', 4);
% duplicateDataType(VOLUME{1}, 'OddScansNull', 4);
% duplicateDataType(VOLUME{1}, 'EvenScansNull', 4);
% duplicateDataType(VOLUME{1}, 'FirstHalfNull', 4);
% duplicateDataType(VOLUME{1}, 'SecondHalfNull', 4);

%Save sweeps to file as separate scans, for all scans together and
%cross-validation dataTYPES.
%After doing this, duplicate dataTYPE and collapse tSeries to layer 1
%voxels.
nScans=size(tSeriesOutAll,3);

S8indices=1:24487;
S9indices=24488:52685;
S10indices=52686:78854;
S11indices=78855:103572;
S12indices=103573:129695;
S13indices=129696:163131;

thisSubInds=S13indices;

for thisSweep=1:4
    load(['Gray/TimingSweepsNull/TSeries/Scan', num2str(thisSweep),'/tSeries', num2str(1), '.mat']);
    tSeries(:,thisSubInds)=mean(Sweeps{thisSweep}.data,3);
    save(['Gray/TimingSweepsNull/TSeries/Scan', num2str(thisSweep),'/tSeries', num2str(1), '.mat'], 'tSeries')
end

for thisSweep=1:4
    load(['Gray/OddScansNull/TSeries/Scan', num2str(thisSweep),'/tSeries', num2str(1), '.mat']);
    tSeries(:,thisSubInds)=mean(Sweeps{thisSweep}.data(:,:,(1:2:nScans)),3);
    save(['Gray/OddScansNull/TSeries/Scan', num2str(thisSweep),'/tSeries', num2str(1), '.mat'], 'tSeries')
end

for thisSweep=1:4
    load(['Gray/EvenScansNull/TSeries/Scan', num2str(thisSweep),'/tSeries', num2str(1), '.mat']);
    tSeries(:,thisSubInds)=mean(Sweeps{thisSweep}.data(:,:,(2:2:nScans)),3);
    save(['Gray/EvenScans/TSeries/Scan', num2str(thisSweep),'/tSeries', num2str(1), '.mat'], 'tSeries')
end

for thisSweep=1:4
    load(['Gray/FirstHalfNull/TSeries/Scan', num2str(thisSweep),'/tSeries', num2str(1), '.mat']);
    tSeries(:,thisSubInds)=mean(Sweeps{thisSweep}.data(:,:,1:round(nScans/2)),3);
    save(['Gray/FirstHalfNull/TSeries/Scan', num2str(thisSweep),'/tSeries', num2str(1), '.mat'], 'tSeries')
end

for thisSweep=1:4
    load(['Gray/SecondHalfNull/TSeries/Scan', num2str(thisSweep),'/tSeries', num2str(1), '.mat']);
    tSeries(:,thisSubInds)=mean(Sweeps{thisSweep}.data(:,:,(round(nScans/2)+1):nScans),3);
    save(['Gray/SecondHalfNull/TSeries/Scan', num2str(thisSweep),'/tSeries', num2str(1), '.mat'], 'tSeries')
end

%After running model, extract all of each subject's VE values (for all
%modelled voxel, real data)
NellVeAllSmall=rmGet(model{1}, 've');
NellVeAllLarge=rmGet(model{1}, 've');
NellAllSmall=NellAllSmall(roiIndices);
NellVeAllSmall=NellVeAllSmall(roiIndices);
NellVeAllLarge=NellVeAllLarge(roiIndices);
[tmp roiIndices]=intersectCols(VOLUME{1}.coords, VOLUME{1}.ROIs(1).coords);
HarvVeAllSmall=rmGet(model{1}, 've');
HarvVeAllLarge=rmGet(model{1}, 've');
HarvVeAllSmall=HarvVeAllSmall(roiIndices);
HarvVeAllLarge=HarvVeAllLarge(roiIndices);
[tmp roiIndices]=intersectCols(VOLUME{1}.coords, VOLUME{1}.ROIs(1).coords);
DumoVeAllSmall=rmGet(model{1}, 've');
DumoVeAllLarge=rmGet(model{1}, 've');
DumoVeAllSmall=DumoVeAllSmall(roiIndices);
DumoVeAllLarge=DumoVeAllLarge(roiIndices);
[tmp roiIndices]=intersectCols(VOLUME{1}.coords, VOLUME{1}.ROIs(1).coords);
FracVeAllSmall=rmGet(model{1}, 've');
FracVeAllLarge=rmGet(model{1}, 've');
FracVeAllSmall=FracVeAllSmall(roiIndices);
FracVeAllLarge=FracVeAllLarge(roiIndices);
AllVeAllSmall=[NellVeAllSmall HarvVeAllSmall DumoVeAllSmall FracVeAllSmall];
AllVeAllLarge=[NellVeAllLarge HarvVeAllLarge DumoVeAllLarge FracVeAllLarge];

%Now get null VEs
NullVeAll=rmGet(model{1}, 've');
[tmp roiIndices]=intersectCols(VOLUME{1}.coords, VOLUME{1}.ROIs(1).coords);
NullVeAll=NullVeAll(roiIndices);

%And compare
NullVeAll=sort(NullVeAll);
setToTest=AllVeAllLarge;
setToTest=sort(setToTest);
pVals=zeros(size(setToTest));
for n=1:length(setToTest)
    pVals(n)=sum(NullVeAll>setToTest(n))/length(NullVeAll);
end
pVals(pVals<=(1/length(NullVeAll)))=1/length(NullVeAll);

[h, crit_p, adj_ci_cvrg, pValFDR]=fdr_bh(pVals);
whichVE=0.38;
critIndex=find(setToTest>whichVE, 1, 'first');
critpValFDR=pValFDR(critIndex);
