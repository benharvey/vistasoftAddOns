%First, group all scans into a single mrVista gray view data type

%Subject-specific scan orders, randomly shuffled. Commas separate sessions.
S8order=[22 10 6 16 17 3, 9 20 14 23 13 11 15, 4 19 24 2 12 1 7, 21 18 5 8 1]; %Added 1 (1234) for extra scan
S10order=[23 9 24 1 13 3, 2 4 11 15 12 10, 8 7 18 19 16, 17 14 6 20 5 21 22];
S11order=[8 16 10 11 6 21, 9 14 18 17 23, 22 24 5 1 13 7, 3 2 19 20 12 15 8]; %Repeated number 8 due to excessive motion. Dropped scan 4 (number 4 in session 2) from bad data
S12order=[22 6 3 16 11 7, 17 14 8 5 21 19 15, 1 23 2 4 18, 24 13 9 20 10 12 1];
S13order=[7 20 21 6 12 24 14 22 3 10 8 15 19 11 5 4 16 1 13 2 23 17 18 9];
S9order=[7 16 10 11 17 20 6, 15 8 21 5 9 1 2, 19 22 12 18 24 23 14, 13 3 4 1 24 11 15];

%SPECIFY WHICH SUBJECT
ThisSubjectOrder=S9order;


% 1='Timing pRF TR=2.1, Duration Constant Luminance';
% 2='Timing pRF TR=2.1, Duration Constant Event, Variable Interval';
% 3='Timing pRF TR=2.1, Duration Constant Frequency';
% 4='Timing pRF TR=2.1, Duration Gaps';

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

%Make a new dataTYPE for average sweep data, and cross-validation dataTYPES
groupScans(VOLUME{1}, 1:4, 'TimingSweeps');
VOLUME{1}=viewSet(VOLUME{1},'curdt',4);
tSeriesClipFrames(VOLUME{1},1:4,0,56);
duplicateDataType(VOLUME{1}, 'OddScans', 4);
duplicateDataType(VOLUME{1}, 'EvenScans', 4);
duplicateDataType(VOLUME{1}, 'FirstHalf', 4);
duplicateDataType(VOLUME{1}, 'SecondHalf', 4);

VOLUME{1}=viewSet(VOLUME{1},'curdt','MotionComp');
nScans=size(dataTYPES(3).scanParams, 2);

%Initialize scan averaging structure
Sweeps{1}.data=zeros([56, size(VOLUME{1}.coords, 2), nScans]);
Sweeps{2}.data=zeros([56, size(VOLUME{1}.coords, 2), nScans]);
Sweeps{3}.data=zeros([56, size(VOLUME{1}.coords, 2), nScans]);
Sweeps{4}.data=zeros([56, size(VOLUME{1}.coords, 2), nScans]);

%For each scan
for thisScan=1:nScans
    %Load time series
    tSeries = loadtSeries(VOLUME{1}, thisScan, 1);

    %Extract the 4 sweeps
    tSeriesSweeps= reshape(tSeries, [56, 4, size(VOLUME{1}.coords, 2)]);
    for thisSweep=1:4
        Sweeps{Scan2Sweep(ThisSubjectOrder(thisScan),thisSweep)}.data(:, :, thisScan)=squeeze(tSeriesSweeps(:,thisSweep,:));
    end
end

%Save sweeps to file as separate scans, for all scans together and
%cross-validation dataTYPES.
%After doing this, duplicate dataTYPE and collapse tSeries to layer 1
%voxels.
VOLUME{1}=viewSet(VOLUME{1},'curdt','TimingSweeps');
for thisSweep=1:4
    savetSeries(mean(Sweeps{thisSweep}.data,3),VOLUME{1},thisSweep,1);   
end
collapseOverLayersTseries(VOLUME{1});

VOLUME{1}=viewSet(VOLUME{1},'curdt','OddScans');
for thisSweep=1:4
    savetSeries(mean(Sweeps{thisSweep}.data(:,:,(1:2:nScans)),3),VOLUME{1},thisSweep,1);   
end
collapseOverLayersTseries(VOLUME{1});

VOLUME{1}=viewSet(VOLUME{1},'curdt','EvenScans');
for thisSweep=1:4
    savetSeries(mean(Sweeps{thisSweep}.data(:,:,(2:2:nScans)),3),VOLUME{1},thisSweep,1);   
end
collapseOverLayersTseries(VOLUME{1});

VOLUME{1}=viewSet(VOLUME{1},'curdt','FirstHalf');
for thisSweep=1:4
    savetSeries(mean(Sweeps{thisSweep}.data(:,:,1:round(nScans/2)),3),VOLUME{1},thisSweep,1);   
end
collapseOverLayersTseries(VOLUME{1});

VOLUME{1}=viewSet(VOLUME{1},'curdt','SecondHalf');
for thisSweep=1:4
    savetSeries(mean(Sweeps{thisSweep}.data(:,:,(round(nScans/2)+1):nScans),3),VOLUME{1},thisSweep,1);   
end
collapseOverLayersTseries(VOLUME{1});