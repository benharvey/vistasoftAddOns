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

%Create an ROI of the layer 1 voxels and save to local file
VOLUME{1}=makeGrayROI(VOLUME{1}); VOLUME{1}=refreshScreen(VOLUME{1},0);
VOLUME{1} = roiRestricttoLayer1(VOLUME{1},VOLUME{1}.selectedROI); 
VOLUME{1} = refreshScreen(VOLUME{1},0);
saveROI(VOLUME{1}, 'selected', 1);

%Make parameter definitions for models
MakeDurationVsPeriod;

%2-dimensional duration vs period models (mixed scans)
load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsMix.mat')

combinedDT=9
setAllRetParams(paramsDurationPeriod, combinedDT);
rmRunDuration2d(VOLUME{1},combinedDT, 'gray-Layer1',8,{'1g'},[],[],[],0);
rmRunDuration2dOval(VOLUME{1},combinedDT, 'gray-Layer1',8,{'1g'},[],[],[],0);
rmRunLogDuration2d(VOLUME{1},combinedDT{1}, 'gray-Layer1',8,{'1g'},[],[],[],0);
rmRunLogDuration2dOval(VOLUME{1},combinedDT, 'gray-Layer1',8,{'1g'},[],[],[],0);

combinedDT=[10 11]
setAllRetParams(paramsDurationPeriod, combinedDT);
rmRunDuration2d(VOLUME{1},combinedDT, 'gray-Layer1',8,{'1g'},[],[],[],0);
rmRunLogDuration2d(VOLUME{1},combinedDT, 'gray-Layer1',8,{'1g'},[],[],[],0);
rmRunDuration2dOval(VOLUME{1},combinedDT, 'gray-Layer1',8,{'1g'},[],[],[],0);
rmRunLogDuration2dOval(VOLUME{1},combinedDT, 'gray-Layer1',8,{'1g'},[],[],[],0);

combinedDT=[9:11]
setAllRetParams(paramsDurationPeriod, combinedDT);
rmRunLogDuration2d(VOLUME{1},combinedDT, 'gray-Layer1',1,{'1g'},[],[],[],0);
rmRunDuration2d(VOLUME{1},combinedDT, 'gray-Layer1',1,{'1g'},[],[],[],0);
rmRunDuration2dOval(VOLUME{1},combinedDT, 'gray-Layer1',1,{'1g'},[],[],[],0);
rmRunLogDuration2dOval(VOLUME{1},combinedDT, 'gray-Layer1',1,{'1g'},[],[],[],0);

combinedDT=[9:11]
setAllRetParams(paramsDurationPeriod, combinedDT);
rmRunDuration2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0);
rmRunDuration2d(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0);
rmRunLogDuration2d(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0);
rmRunLogDuration2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0);

%2-dimensional duration vs off time models
combinedDT=9;
setAllRetParams(paramsOnTimeOffTime, combinedDT);
rmRunOnTimeOffTime2d(VOLUME{1},combinedDT, 'gray-Layer1',1,{'1g'},[],[],[],0);
rmRunLogOnTimeOffTime2d(VOLUME{1},combinedDT, 'gray-Layer1',1,{'1g'},[],[],[],0);
rmRunOnTimeOffTime2dOval(VOLUME{1},combinedDT, 'gray-Layer1',1,{'1g'},[],[],[],0);
rmRunLogOnTimeOffTime2dOval(VOLUME{1},combinedDT, 'gray-Layer1',1,{'1g'},[],[],[],0);

combinedDT=[9 10 11]
setAllRetParams(paramsOnTimeOffTime, combinedDT);
rmRunLogOnTimeOffTime2d(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0);
rmRunLogOnTimeOffTime2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0);
rmRunOnTimeOffTime2d(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0);
rmRunOnTimeOffTime2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0);

%2-dimensional Occupancy vs period models
combinedDT=[9:11];
setAllRetParams(paramsOccupancyPeriod, combinedDT)
rmRunOccupancyPeriod2d(VOLUME{1},combinedDT, 'gray-Layer1',1,{'1g'},[],[],[],0);
rmRunLogOccupancyPeriod2d(VOLUME{1},combinedDT, 'gray-Layer1',1,{'1g'},[],[],[],0);
rmRunOccupancyPeriod2dOval(VOLUME{1},combinedDT, 'gray-Layer1',1,{'1g'},[],[],[],0);
rmRunLogOccupancyPeriod2dOval(VOLUME{1},combinedDT, 'gray-Layer1',1,{'1g'},[],[],[],0);

combinedDT=[9:11]
setAllRetParams(paramsOccupancyPeriod, combinedDT)
rmRunLogOccupancyPeriod2d(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0);
rmRunLogOccupancyPeriod2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0);
rmRunOccupancyPeriod2d(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0);
rmRunOccupancyPeriod2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0);

%With log intensity of stimulus
setAllRetParams(paramsDurationPeriod, combinedDT);
rmRunDuration2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0,1);
rmRunDuration2d(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0,1);
rmRunLogDuration2d(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0,1);
rmRunLogDuration2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0,1);
rmRunLinDurNormLogPer2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0, 1);

setAllRetParams(paramsOccupancyPeriod, combinedDT)
rmRunLogOccupancyPeriod2d(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0,1);
rmRunLogOccupancyPeriod2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0,1);
rmRunOccupancyPeriod2d(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0,1);
rmRunOccupancyPeriod2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0,1);
rmRunOccupancyNormLogPeriod2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0,1);

%Further possibilities
combinedDT=9:11
setAllRetParams(paramsDurationPeriod, combinedDT);
%rmRunInvDurFreq2dOval(VOLUME{1},combinedDT, 'TPC1LeftMap',4,{'1g'},[],[],[],0);
%rmRunLogInvDurFreq2dOval(VOLUME{1},combinedDT, 'TPC1LeftMap',4,{'1g'},[],[],[],0);
%rmRunLogDurFreq2dOval(VOLUME{1},combinedDT, 'TPC1LeftMap',4,{'1g'},[],[],[],0);
rmRunLinDurLogPer2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0);
rmRunLinDurNormLogPer2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0);

rmRunDurFreq2d(VOLUME{1},combinedDT(1), 'TPC1LeftMap',12,{'1g'},[],[],[],0); %Duration tuning with compressive monotonic frequency

setAllRetParams(paramsOccupancyPeriod, combinedDT)
rmRunOccupancyLogPeriod2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0); %Good 
rmRunOccupancyNormLogPeriod2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0);

%For testing
rmRunLinDurNormLogPer2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0);
setAllRetParams(paramsDurationPeriod, combinedDT);
rmRunDuration2dOval(VOLUME{1},combinedDT(1), 'TPC1LeftMap',4,{'1g'},[],[],[],0);
rmRunLinDurNormLogPer2dOval(VOLUME{1},combinedDT(1), 'TPC1LeftMap',4,{'1g'},[],[],[],0);
rmRunDuration2dOval(VOLUME{1},combinedDT(1), 'TPC1LeftMap',4,{'1g'},[],[],[],0,'free');

%Monotonic
combinedDT=9:11;
setAllRetParams(paramsDurationPeriod, combinedDT);
rmRunDuration2d(VOLUME{1},combinedDT, 'gray-Layer1',9,{'1g'},[],[],[],0);
rmRunLinDurNormLogPer2dOval(VOLUME{1},combinedDT, 'gray-Layer1',9,{'1g'},[],[],[],0);

%Loop to run for all subjects
paths{1}='/media/Storage3/WholeBrainNewSeg/Subject 8/S8_Combined';
paths{2}='/media/Storage3/WholeBrainNewSeg/Subject 9/S9_Combined';
paths{3}='/media/Storage3/WholeBrainNewSeg/Subject 10/S10_Combined';
paths{4}='/media/Storage3/WholeBrainNewSeg/Subject 11/S11_Combined';
paths{5}='/media/Storage3/WholeBrainNewSeg/Subject 12/S12_Combined';
paths{6}='/media/Storage3/WholeBrainNewSeg/Subject 13/S13_Combined';
paths{7}='/media/Storage3/WholeBrainNewSeg/Harv/HarvCombined';
paths{8}='/media/Storage3/WholeBrainNewSeg/Nell/NellCombined';

load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsMix.mat')
combinedDT=9:11;
whichSubs=[1 2 3 4 5 6];
maxCores=3;
for thisSub=1[8, 7, 1];%1:length(whichSubs)
   cd(paths{thisSub})
   mrVista 3;
   
   %The code to run for each
   if thisSub<=6
       combinedDT=9:11;
       load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsMix.mat')
       setAllRetParams(paramsDurationPeriod, combinedDT);
       rmRunDuration2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],maxCores,0);
       rmRunDuration2d(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],maxCores,0);
       rmRunLogDuration2d(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],maxCores,0);
       rmRunLogDuration2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],maxCores,0);
   else
       if thisSub==7
           combinedDT=[67 69 70]; %For Harv
       elseif thisSub==8
           combinedDT=[72 84 85];
       end
       load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsSessions.mat')
       setAllRetParams(paramsDurationPeriod, combinedDT);
       rmRunDuration2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],maxCores,1);
       rmRunDuration2d(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],maxCores,1);
       rmRunLogDuration2d(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],maxCores,1);
       rmRunLogDuration2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],maxCores,1);
   end
   close(1); mrvCleanWorkspace; 
end



%With log intensity of stimulus
whichSubs=[4 5 8]%2 3 7]% 4 5 8];
for thisSub=whichSubs
    cd(paths{thisSub})
    mrVista 3;
    
    %The code to run for each
    if thisSub<=6
        combinedDT=9:11;
        load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsMix.mat')
        setAllRetParams(paramsDurationPeriod, combinedDT);
        rmRunDuration2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0,1);
        rmRunDuration2d(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0,1);
        rmRunLogDuration2d(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0,1);
        rmRunLogDuration2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0,1);
        rmRunLinDurNormLogPer2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0, 1);
    else
        if thisSub==7
            combinedDT=[67 69 70]; %For Harv
        elseif thisSub==8
            combinedDT=[72 84 85];
        end
        load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsSessions.mat')
        setAllRetParams(paramsDurationPeriod, combinedDT);
        rmRunDuration2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],1,1);
        rmRunDuration2d(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],1,1);
        rmRunLogDuration2d(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],1,1);
        rmRunLogDuration2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],1,1);
        rmRunLinDurNormLogPer2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],1, 1);
    end
    close(1); mrvCleanWorkspace;
end

whichSubs=[3 5 7 8]
for thisSub=whichSubs
    cd(paths{thisSub})
    mrVista 3;
    
    %The code to run for each
    if thisSub<=6
        combinedDT=9:11;
        load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsMix.mat')
        setAllRetParams(paramsOccupancyPeriod, combinedDT)
        rmRunLogOccupancyPeriod2d(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0,1);
        rmRunLogOccupancyPeriod2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0,1);
        rmRunOccupancyPeriod2d(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0,1);
        rmRunOccupancyPeriod2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0,1);
        rmRunOccupancyNormLogPeriod2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0,1);
    else
        if thisSub==7
            combinedDT=[67 69 70]; %For Harv
        elseif thisSub==8
            combinedDT=[72 84 85];
        end
        load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsSessions.mat')
        setAllRetParams(paramsOccupancyPeriod, combinedDT)
        rmRunLogOccupancyPeriod2d(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],1,1);
        rmRunLogOccupancyPeriod2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],1,1);
        rmRunOccupancyPeriod2d(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],1,1);
        rmRunOccupancyPeriod2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],1,1);
        rmRunOccupancyNormLogPeriod2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],1,1);
    end
    close(1); mrvCleanWorkspace;
end

%With free compressive exponent.
%Duration period
whichSubs=[4:6 8]%[6 5 4 7]%2 3 7]% 4 5 8];
for thisSub=whichSubs
    cd(paths{thisSub})
    mrVista 3;
    
    %The code to run for each
    if thisSub<=6
        combinedDT=9:11;
        load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsMix.mat')
        setAllRetParams(paramsDurationPeriod, combinedDT);
        rmRunDuration2dOval(VOLUME{1},combinedDT, 'gray-Layer1',5,{'1g'},[],[],[],0,'free');

    else
        if thisSub==7
            combinedDT=[67 69 70]; %For Harv
        elseif thisSub==8
            combinedDT=[72 84 85];
        end
        load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsSessions.mat')
        setAllRetParams(paramsDurationPeriod, combinedDT);
        rmRunDuration2dOval(VOLUME{1},combinedDT, 'gray-Layer1',5,{'1g'},[],[],[],1,'free');
    end
    close(1); mrvCleanWorkspace;
end

%logDuration logPeriod
whichSubs=[7]%[6 5 4 7]%2 3 7]% 4 5 8];
for thisSub=whichSubs
    cd(paths{thisSub})
    mrVista 3;
    
    %The code to run for each
    if thisSub<=6
        combinedDT=9:11;
        load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsMix.mat')
        setAllRetParams(paramsDurationPeriod, combinedDT);
        rmRunLogDuration2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0,'free');

    else
        if thisSub==7
            combinedDT=[67 69 70]; %For Harv
        elseif thisSub==8
            combinedDT=[72 84 85];
        end
        load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsSessions.mat')
        setAllRetParams(paramsDurationPeriod, combinedDT);
        rmRunLogDuration2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],1,'free');
    end
    close(1); mrvCleanWorkspace;
end

%Occupancy period, free exponent
whichSubs=[7]%[6 5 4 7]%2 3 7]% 4 5 8];
for thisSub=whichSubs
    cd(paths{thisSub})
    mrVista 3;
    
    %The code to run for each
    if thisSub<=6
        combinedDT=9:11;
        load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsMix.mat')
        setAllRetParams(paramsOccupancyPeriod, combinedDT);
        rmRunOccupancyPeriod2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0,'free');

    else
        if thisSub==7
            combinedDT=[67 69 70]; %For Harv
        elseif thisSub==8
            combinedDT=[72 84 85];
        end
        load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsSessions.mat')
        setAllRetParams(paramsOccupancyPeriod, combinedDT);
        rmRunOccupancyPeriod2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],1,'free');
    end
    close(1); mrvCleanWorkspace;
end

%OnTime OffTime, free exponent
whichSubs=[8]%[6 5 4 7]%2 3 7]% 4 5 8];
for thisSub=whichSubs
    cd(paths{thisSub})
    mrVista 3;
    
    %The code to run for each
    if thisSub<=6
        combinedDT=9:11;
        load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsMix.mat')
        setAllRetParams(paramsOnTimeOffTime, combinedDT);
        rmRunOnTimeOffTime2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0,'free');

    else
        if thisSub==7
            combinedDT=[67 69 70]; %For Harv
        elseif thisSub==8
            combinedDT=[72 84 85];
        end
        load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsSessions.mat')
        setAllRetParams(paramsOnTimeOffTime, combinedDT);
        rmRunOnTimeOffTime2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],1,'free');
    end
    close(1); mrvCleanWorkspace;
end


%Duration period, models with monotonic components
whichSubs=[4:6 8]%[6 5 4 7]%2 3 7]% 4 5 8];
for thisSub=whichSubs
    cd(paths{thisSub})
    mrVista 3;
    
    %The code to run for each
    if thisSub<=6
        combinedDT=9:11;
        load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsMix.mat')
        setAllRetParams(paramsDurationPeriod, combinedDT);
        for whichModel=14:16
            rmRunDurFreq2d(VOLUME{1},combinedDT, 'gray-Layer1',whichModel,{'1g'},[],[],[],0);
        end

    else
        if thisSub==7
            combinedDT=[67 69 70]; %For Harv
        elseif thisSub==8
            combinedDT=[72 84 85];
        end
        load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsSessions.mat')
        setAllRetParams(paramsDurationPeriod, combinedDT);
        for whichModel=14:16
            rmRunDurFreq2d(VOLUME{1},combinedDT, 'gray-Layer1',whichModel,{'1g'},[],[],[],1);
        end
    end
    close(1); mrvCleanWorkspace;
end

%Occupancy period, models with monotonic components
whichSubs=[2:8]%[6 5 4 7]%2 3 7]% 4 5 8];
for thisSub=whichSubs
    cd(paths{thisSub})
    mrVista 3;
    
    %The code to run for each
    if thisSub<=6
        combinedDT=9:11;
        load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsMix.mat')
        setAllRetParams(paramsOccupancyPeriod, combinedDT);
        for whichModel=[14]
            rmRunOccupancyFreq2d(VOLUME{1},combinedDT, 'gray-Layer1',whichModel,{'1g'},[],[],[],0);
        end

    else
        if thisSub==7
            combinedDT=[67 69 70]; %For Harv
        elseif thisSub==8
            combinedDT=[72 84 85];
        end
        load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsSessions.mat')
        setAllRetParams(paramsOccupancyPeriod, combinedDT);
        for whichModel=[14]
            rmRunOccupancyFreq2d(VOLUME{1},combinedDT, 'gray-Layer1',whichModel,{'1g'},[],[],[],1);
        end
    end
    close(1); mrvCleanWorkspace;
end

%Temporal frequency tuned
whichSubs=[1:3 7]%[6 5 4 7]%2 3 7]% 4 5 8];
for thisSub=whichSubs
    cd(paths{thisSub})
    mrVista 3;
    
    %The code to run for each
    if thisSub<=6
        combinedDT=9:11;
        load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsMix.mat')
        setAllRetParams(paramsTemporalFrequency, combinedDT);
        for whichModel=[14]
            rmRunTemporalFreq1d(VOLUME{1},combinedDT, 'gray-Layer1',whichModel,{'1g'},[],[],[],0);
        end

    else
        if thisSub==7
            combinedDT=[67 69 70]; %For Harv
        elseif thisSub==8
            combinedDT=[72 84 85];
        end
        load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsSessions.mat')
        setAllRetParams(paramsTemporalFrequency, combinedDT);
        for whichModel=[14]
            rmRunTemporalFreq1d(VOLUME{1},combinedDT, 'gray-Layer1',whichModel,{'1g'},[],[],[],1);
        end
    end
    close(1); mrvCleanWorkspace;
end

%After correcting even scan split
load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsMix.mat')

combinedDT=11;
setAllRetParams(paramsDurationPeriod, combinedDT);
rmRunDuration2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0);
rmRunLogDuration2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0);
rmRunLinDurNormLogPer2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0);

setAllRetParams(paramsOnTimeOffTime, combinedDT);
rmRunLogOnTimeOffTime2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0);
rmRunOnTimeOffTime2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0);

setAllRetParams(paramsOccupancyPeriod, combinedDT)
rmRunLogOccupancyPeriod2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0);
rmRunOccupancyPeriod2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0);
rmRunOccupancyNormLogPeriod2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0);

setAllRetParams(paramsDurationPeriod, combinedDT);
rmRunDuration2d(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0);
rmRunLogDuration2d(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0);

setAllRetParams(paramsOnTimeOffTime, combinedDT);
rmRunLogOnTimeOffTime2d(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0);
rmRunOnTimeOffTime2d(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0);

setAllRetParams(paramsOccupancyPeriod, combinedDT)
rmRunLogOccupancyPeriod2d(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0);
rmRunOccupancyPeriod2d(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0);

%For null distributions
combinedDT=14:16;
setAllRetParams(paramsDurationPeriod, combinedDT);
rmRunDuration2dOval(VOLUME{1},combinedDT, 'NullDataROI',4,{'1g'},[],[],[],0);
rmRunLogDuration2dOval(VOLUME{1},combinedDT, 'NullDataROI',4,{'1g'},[],[],[],0);
rmRunLinDurNormLogPer2dOval(VOLUME{1},combinedDT, 'NullDataROI',4,{'1g'},[],[],[],0);

setAllRetParams(paramsOnTimeOffTime, combinedDT);
rmRunLogOnTimeOffTime2dOval(VOLUME{1},combinedDT, 'NullDataROI',4,{'1g'},[],[],[],0);
rmRunOnTimeOffTime2dOval(VOLUME{1},combinedDT, 'NullDataROI',4,{'1g'},[],[],[],0);

setAllRetParams(paramsOccupancyPeriod, combinedDT)
rmRunLogOccupancyPeriod2dOval(VOLUME{1},combinedDT, 'NullDataROI',4,{'1g'},[],[],[],0);
rmRunOccupancyPeriod2dOval(VOLUME{1},combinedDT, 'NullDataROI',4,{'1g'},[],[],[],0);
rmRunOccupancyNormLogPeriod2dOval(VOLUME{1},combinedDT, 'NullDataROI',4,{'1g'},[],[],[],0);

setAllRetParams(paramsDurationPeriod, combinedDT);
rmRunDuration2d(VOLUME{1},combinedDT, 'NullDataROI',4,{'1g'},[],[],[],0);
rmRunLogDuration2d(VOLUME{1},combinedDT, 'NullDataROI',4,{'1g'},[],[],[],0);

setAllRetParams(paramsOnTimeOffTime, combinedDT);
rmRunLogOnTimeOffTime2d(VOLUME{1},combinedDT, 'NullDataROI',4,{'1g'},[],[],[],0);
rmRunOnTimeOffTime2d(VOLUME{1},combinedDT, 'NullDataROI',4,{'1g'},[],[],[],0);

setAllRetParams(paramsOccupancyPeriod, combinedDT)
rmRunLogOccupancyPeriod2d(VOLUME{1},combinedDT, 'NullDataROI',4,{'1g'},[],[],[],0);
rmRunOccupancyPeriod2d(VOLUME{1},combinedDT, 'NullDataROI',4,{'1g'},[],[],[],0);

setAllRetParams(paramsDurationPeriod, combinedDT);
rmRunDuration2dOval(VOLUME{1},combinedDT, 'NullDataROI',4,{'1g'},[],[],[],0,'free');

folderName=[pwd '/Gray/' dataTYPES(15).name];
modelFiles=dir([folderName '/*free-gFit.mat']);
rmMainPostGrid([1 15],'NullDataROI',4, [folderName '/' modelFiles(1).name]);

%For seperate sessions
for n=1:4
    paramsDurationPeriod(n).nCycles=3;
    paramsDurationPeriod(n).nUniqueRep=3;
    paramsDurationPeriod(n).nFrames=168;
    paramsOccupancyPeriod(n).nCycles=3;
    paramsOccupancyPeriod(n).nUniqueRep=3;
    paramsOccupancyPeriod(n).nFrames=168;
    paramsOnTimeOffTime(n).nCycles=3;
    paramsOnTimeOffTime(n).nUniqueRep=3;
    paramsOnTimeOffTime(n).nFrames=168;
end

paramsDurationPeriod(1).paramsFile='/media/Storage3/WholeBrainNewSeg/Nell/NellCombined/Stimuli/params_ConLumOff2D.mat';
paramsOccupancyPeriod(1).paramsFile='/media/Storage3/WholeBrainNewSeg/Nell/NellCombined/Stimuli/params_ConLumOff2D.mat';
paramsOnTimeOffTime(1).paramsFile='/media/Storage3/WholeBrainNewSeg/Nell/NellCombined/Stimuli/params_ConLumOff2D.mat';

paramsDurationPeriod(2).paramsFile='/media/Storage3/WholeBrainNewSeg/Nell/NellCombined/Stimuli/params_ConDurOff2D.mat';
paramsOccupancyPeriod(2).paramsFile='/media/Storage3/WholeBrainNewSeg/Nell/NellCombined/Stimuli/params_ConDurOff2D.mat';
paramsOnTimeOffTime(2).paramsFile='/media/Storage3/WholeBrainNewSeg/Nell/NellCombined/Stimuli/params_ConDurOff2D.mat';

paramsDurationPeriod(3).paramsFile='/media/Storage3/WholeBrainNewSeg/Nell/NellCombined/Stimuli/params_ConPerOff2D.mat';
paramsOccupancyPeriod(3).paramsFile='/media/Storage3/WholeBrainNewSeg/Nell/NellCombined/Stimuli/params_ConPerOff2D.mat';
paramsOnTimeOffTime(3).paramsFile='/media/Storage3/WholeBrainNewSeg/Nell/NellCombined/Stimuli/params_ConPerOff2D.mat';

paramsDurationPeriod(4).paramsFile='/media/Storage3/WholeBrainNewSeg/Nell/NellCombined/Stimuli/params_TimeGapsOff2D.mat';
paramsOccupancyPeriod(4).paramsFile='/media/Storage3/WholeBrainNewSeg/Nell/NellCombined/Stimuli/params_TimeGapsOff2D.mat';
paramsOnTimeOffTime(4).paramsFile='/media/Storage3/WholeBrainNewSeg/Nell/NellCombined/Stimuli/params_TimeGapsOff2D.mat';
   
load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsSessions.mat')

combinedDT=[67 69 70] %For Harv
combinedDT=[72 84 85] %For Nell
setAllRetParams(paramsDurationPeriod, combinedDT);
rmRunLogDuration2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],1);
rmRunDuration2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],1);
rmRunDuration2d(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],1);
rmRunLogDuration2d(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],1);

combinedDT=[84 85]
setAllRetParams(paramsOccupancyPeriod, combinedDT)
rmRunOccupancyPeriod2d(VOLUME{1},combinedDT, 'gray-Layer1',1,{'1g'},[],[],[],1);
rmRunOccupancyPeriod2dOval(VOLUME{1},combinedDT, 'gray-Layer1',1,{'1g'},[],[],[],1);
rmRunLogOccupancyPeriod2d(VOLUME{1},combinedDT, 'gray-Layer1',1,{'1g'},[],[],[],1);
rmRunLogOccupancyPeriod2dOval(VOLUME{1},combinedDT, 'gray-Layer1',1,{'1g'},[],[],[],1);

combinedDT=[72 84 85]
setAllRetParams(paramsOnTimeOffTime, combinedDT);
rmRunLogOnTimeOffTime2d(VOLUME{1},combinedDT, 'gray-Layer1',1,{'1g'},[],[],[],1);
rmRunLogOnTimeOffTime2dOval(VOLUME{1},combinedDT, 'gray-Layer1',1,{'1g'},[],[],[],1);
rmRunOnTimeOffTime2d(VOLUME{1},combinedDT, 'gray-Layer1',1,{'1g'},[],[],[],1);
rmRunOnTimeOffTime2dOval(VOLUME{1},combinedDT, 'gray-Layer1',1,{'1g'},[],[],[],1);

setAllRetParams(paramsDurationPeriod, combinedDT);
%rmRunInvDurFreq2dOval(VOLUME{1},combinedDT, 'TPC1LeftMap',4,{'1g'},[],[],[],0);
%rmRunLogInvDurFreq2dOval(VOLUME{1},combinedDT, 'TPC1LeftMap',4,{'1g'},[],[],[],0);
%rmRunLogDurFreq2dOval(VOLUME{1},combinedDT, 'TPC1LeftMap',4,{'1g'},[],[],[],0);
rmRunLinDurLogPer2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],1);
rmRunLinDurNormLogPer2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],1);
rmRunDurFreq2d(VOLUME{1},combinedDT(1), 'TPC1LeftMap',12,{'1g'},[],[],[],1); %Duration tuning with compressive monotonic frequency


setAllRetParams(paramsOccupancyPeriod, combinedDT)
rmRunOccupancyLogPeriod2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],1); %Good fit
rmRunOccupancyNormLogPeriod2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],1); %Good fit


%Monotonic
setAllRetParams(paramsDurationPeriod, combinedDT);
rmRunDuration2d(VOLUME{1},combinedDT, 'gray-Layer1',9,{'1g'},[],[],[],1);

%To build from existing gfit models to sfit and HRF
combinedDT=[9:11];
if length(combinedDT)>1
    parpool(length(combinedDT))
end
parfor whichDT=1:length(combinedDT)
    folderName=[pwd '/Gray/' dataTYPES(combinedDT(whichDT)).name];
    modelFiles=dir([folderName '/*-gFit.mat']);
    for whichModel=1:length(modelFiles)
        rmMainPostGrid([1 combinedDT(whichDT)],'gray-Layer1',4, [folderName '/' modelFiles(whichModel).name]);
    end
end


whichDT=1
folderName=[pwd '/Gray/' dataTYPES(combinedDT(whichDT)).name];
modelFiles=dir([folderName '/*-fFit-fFit.mat']);
if length(modelFiles)>1 && length(modelFiles)<5
    parpool(length(modelFiles))
elseif length(modelFiles)>1
    parpool(3)
end
parfor whichModel=1:length(modelFiles)
    rmMainPostSearch([1 combinedDT(whichDT)],'gray-Layer1',5, [folderName '/' modelFiles(whichModel).name]);
end


%HRF fit for best models only
for whichDT=1:length(combinedDT)
    clear('modelFiles');
    folderName=[pwd '/Gray/' dataTYPES(combinedDT(whichDT)).name, '/SearchFit'];
    modelFiles(1)=dir([char(strcat(folderName, '/*-LinNormLog-2dOvalGaussian-DurationPeriod*'))]);
    modelFiles(2)=dir([char(strcat(folderName, '/*-LinNormLog-2dOvalGaussian-OccupancyPeriod*'))]);
    if whichDT>1
        modelFiles(3)=dir([char(strcat(folderName, '/*-Lin-2dOvalGaussian-OccupancyPeriod*'))]);
        modelFiles(4)=dir([char(strcat(folderName, '/*-Lin-2dOvalGaussian-DurationPeriod*'))]);
    end
    
    parpool(length(modelFiles))
    parfor whichModel=1:length(modelFiles)
        rmMainPostSearch([1 combinedDT(whichDT)],'gray-Layer1',5, [folderName '/' modelFiles(whichModel).name]);
    end
    delete(gcp('nocreate'))
end

for whichDT=1:length(combinedDT)
    clear('modelFiles');
    folderName=[pwd '/Gray/' dataTYPES(combinedDT(whichDT)).name, '/SearchFit'];
    modelFiles(1)=dir([char(strcat(folderName, '/*-Log-2dCircGaussian-OccupancyPeriod*'))]);
    modelFiles(2)=dir([char(strcat(folderName, '/*-Lin-2dCircGaussian-DurationPeriod*'))]);
    modelFiles(3)=dir([char(strcat(folderName, '/*-Log-2dCircGaussian-DurationPeriod*'))]);
    %modelFiles(3)=dir([char(strcat(folderName, '/*-Lin-2dCircGaussian-OccupancyPeriod*'))]);
    modelFiles(4)=dir([char(strcat(folderName, '/*-Lin-2dOvalGaussian-DurationPeriod*'))]);
    
    
    parpool(4)
    parfor whichModel=1:length(modelFiles)
        rmMainPostSearch([1 combinedDT(whichDT)],'gray-Layer1',5, [folderName '/' modelFiles(whichModel).name]);
    end
    delete(gcp('nocreate'))
end

paths{1}='/media/Storage3/WholeBrainNewSeg/Subject 8/S8_Combined';
paths{2}='/media/Storage3/WholeBrainNewSeg/Subject 9/S9_Combined';
paths{3}='/media/Storage3/WholeBrainNewSeg/Subject 10/S10_Combined';
paths{4}='/media/Storage3/WholeBrainNewSeg/Subject 11/S11_Combined';
paths{5}='/media/Storage3/WholeBrainNewSeg/Subject 12/S12_Combined';
paths{6}='/media/Storage3/WholeBrainNewSeg/Subject 13/S13_Combined';
paths{7}='/media/Storage3/WholeBrainNewSeg/Harv/HarvCombined';
paths{8}='/media/Storage3/WholeBrainNewSeg/Nell/NellCombined';

for thisSub=[5];%1:length(whichSubs)
    cd(paths{thisSub})
    mrVista 3;
    %For a single model
    
    if thisSub==7
        combinedDT=[67 69];% 70]; %For Harv
    elseif thisSub==8
        combinedDT=[72 84 85];
    else
        combinedDT=[9];% 10 11];
    end
    
    parfor whichDT=1:length(combinedDT)
        folderName=[pwd '/Gray/' dataTYPES(combinedDT(whichDT)).name '/SearchFitFreeExponent'];
        modelFiles=dir([char(strcat(folderName, '/*-Lin-2dOvalGaussian-DurationPeriod*expIntensity-free*.mat'))]);
        rmMainPostSearch([1 combinedDT(whichDT)],'gray-Layer1',5, [folderName '/' modelFiles.name]);
    end
    delete(gcp('nocreate'))
    close(1); mrvCleanWorkspace;
end

rmMainPostGrid([1 combinedDT],'ROI10PC',4, [folderName '/' modelFiles.name]);

for thisSub=[3 2 1 7 8];%1:length(whichSubs)
    cd(paths{thisSub})
    mrVista 3;
    if thisSub==7
        combinedDT=[67 69 70]; %For Harv
    elseif thisSub==8
        combinedDT=[72 84 85];
    else
        combinedDT=[9 10 11];
    end
    
    for whichDT= 1:length(combinedDT)
        folderName=[pwd '/Gray/' dataTYPES(combinedDT(whichDT)).name];
        eval(['!mkdir ',  '"',folderName, '/OldOval"']);
        eval(['!mkdir ',  '"',folderName, '/OvalGFit"']);
        
        modelFilesLin=dir([char(strcat(folderName, '/SearchFit/*Oval*'))]);
        modelFilesLog=dir([char(strcat(folderName, '/SearchFitLogIntensity/*Oval*'))]);
        
        for n=1:length(modelFilesLin)
            eval(['!mv ', '"', modelFilesLin(n).folder, '/',  modelFilesLin(n).name, '" "',folderName, '/OldOval/', modelFilesLin(n).name,'"']);
        end
        for n=1:length(modelFilesLog)
            eval(['!mv ', '"', modelFilesLog(n).folder, '/',  modelFilesLog(n).name, '" "',folderName, '/OldOval/', modelFilesLog(n).name,'"']);
        end
        
        for n=1:length(modelFilesLin)
            bit=modelFilesLin(n).name(1:25);
            mf=dir(char(strcat(folderName, '/', bit, '*gFit.mat')));
            eval(['!mv ', '"', mf(1).folder, '/',  mf(1).name, '" "',folderName, '/OvalGFit/', mf(1).name,'"']);
            mf=dir(char(strcat(folderName, '/', bit, '*')));
            for m=1:length(mf)
                eval(['!mv ', '"', mf(m).folder, '/',  mf(m).name, '" "',folderName, '/OldOval/', mf(m).name,'"']);
            end
        end
        
        for n=1:length(modelFilesLog)
            bit=modelFilesLog(n).name(1:25);
            mf=dir(char(strcat(folderName, '/', bit, '*gFit.mat')));
            eval(['!mv ', '"', mf(1).folder, '/',  mf(1).name, '" "',folderName, '/OvalGFit/', mf(1).name,'"']);
            mf=dir(char(strcat(folderName, '/', bit, '*')));
            for m=1:length(mf)
                eval(['!mv ', '"', mf(m).folder, '/',  mf(m).name, '" "',folderName, '/OldOval/', mf(m).name,'"']);
            end
        end
    end
    close(1); mrvCleanWorkspace;
end

for thisSub=[3 2 1];%1:length(whichSubs)
    cd(paths{thisSub})
    mrVista 3;
    %For a single model
    
    if thisSub==7
        combinedDT=[67 69 70]; %For Harv
    elseif thisSub==8
        combinedDT=[72 84 85];
    else
        combinedDT=[9 10 11];
    end
    
    for whichDT=1:length(combinedDT)
        folderName=[pwd '/Gray/' dataTYPES(combinedDT(whichDT)).name, '/OvalGFit'];
        modelFiles=dir([char(strcat(folderName, '/*.mat'))]);
        thisDT=combinedDT(whichDT);
        parfor n=1:length(modelFiles)
            rmMainPostGrid([1 thisDT],'gray-Layer1',4, [folderName '/' modelFiles(n).name]);
        end
    end
    delete(gcp('nocreate'))
    close(1); mrvCleanWorkspace;
end

for thisSub=[7 2];%1:length(whichSubs)
    cd(paths{thisSub})
    mrVista 3;
    %For a single model
    
    if thisSub==7
        combinedDT=[67 69 70]; %For Harv
    elseif thisSub==8
        combinedDT=[72 84 85];
    else
        combinedDT=[9 10 11];
    end
    
    parfor whichDT=1:length(combinedDT)
        folderName=[pwd '/Gray/' dataTYPES(combinedDT(whichDT)).name, '/SearchFitFreeExponent'];
        modelFiles=dir([char(strcat(folderName, '/*OvalGaussian*-fFit.mat'))]);
        rmMainPostSearch([1 combinedDT(whichDT)],'ROI10PC',5, [folderName '/' modelFiles.name]);
    end
    delete(gcp('nocreate'))
    close(1); mrvCleanWorkspace;
end


%FOR AFNI PREPROCESS, RUNNING ALTERNATIVE MODELS
%CB revision
paths{1}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/S08_Combined';
paths{2}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/S09_Combined';
paths{3}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/S10_Combined';
paths{4}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/S11_Combined';
paths{5}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/S12_Combined';
paths{6}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/S13_Combined_Seg20';
paths{7}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/HarvCombined';
paths{8}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/NellCombined';
%Data types for all, odd and even runs
combinedDT=5:7;

%Whole brain combined with numersoity and VFM
paths{1}='/home/benharvey/WholeBrainNewSeg_afni/S8_Combined';
paths{2}='/home/benharvey/WholeBrainNewSeg_afni/S9_Combined';
paths{3}='/home/benharvey/WholeBrainNewSeg_afni/S10_Combined';
paths{4}='/home/benharvey/WholeBrainNewSeg_afni/S11_Combined';
paths{5}='/home/benharvey/WholeBrainNewSeg_afni/S12_Combined';
paths{6}='/media/Storage4/8TB_HDD/WholeBrainNewSeg_afni/S13_Combined_seg40/mrVistaSession';
paths{7}='/media/Storage4/8TB_HDD/WholeBrainNewSeg_afni/Harv_Combined/mrVistaSession';
paths{8}='/home/benharvey/WholeBrainNewSeg_afni/Nell_Combined';

%Data types for all, odd and even runs
combinedDTsAll{1}=19:21;
combinedDTsAll{2}=19:21;
combinedDTsAll{3}=19:21;
combinedDTsAll{4}=19:21;
combinedDTsAll{5}=19:21;
combinedDTsAll{6}=19:21;
combinedDTsAll{7}=23:25;
combinedDTsAll{8}=23:25;


whichSubs=[7]; %To do: 

maxCores=3;
for thisSub=1:length(whichSubs)
   combinedDT=combinedDTsAll{whichSubs(thisSub)};
   cd(paths{whichSubs(thisSub)})
   mrVista 3;
   
   %The code to run for each
   if whichSubs(thisSub)<=6
       load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsMix.mat')
       setAllRetParams(paramsDurationPeriod, combinedDT);
       rmRunDuration2dOval(VOLUME{1},combinedDT, 'gray-Layer1',5,{'1g'},[],[],[],0,'free'); %This is the best fitting model. Setting the 4th input arguement to '5' fits the HRF parameters too.
%        setAllRetParams(paramsDurationPeriod, combinedDT);
%        rmRunLogDuration2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0,'free');
%        setAllRetParams(paramsDurationPeriod, combinedDT);
%        rmRunDuration2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],maxCores,0);
%        
%        for whichModel=[12 14 15 16]
%            setAllRetParams(paramsDurationPeriod, combinedDT);
%            rmRunDurFreq2d(VOLUME{1},combinedDT, 'gray-Layer1',whichModel,{'1g'},[],[],[],0);
%        end
%        
%        setAllRetParams(paramsOccupancyPeriod, combinedDT);
%        rmRunOccupancyPeriod2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0,'free');
%        for whichModel=[12 14]
%            setAllRetParams(paramsOccupancyPeriod, combinedDT);
%            rmRunOccupancyFreq2d(VOLUME{1},combinedDT, 'gray-Layer1',whichModel,{'1g'},[],[],[],0);
%        end
%        
%        setAllRetParams(paramsOnTimeOffTime, combinedDT);
%        rmRunOnTimeOffTime2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0,'free');
%        
%        setAllRetParams(paramsTemporalFrequency, combinedDT);
%        rmRunTemporalFreq1d(VOLUME{1},combinedDT, 'gray-Layer1',14,{'1g'},[],[],[],0);
%        
%        load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsMix.mat')
%        for n=1:4
%            paramsTemporalFrequency(n).paramsFile=[paramsTemporalFrequency(n).paramsFile(1:(end-9)), 'TF', paramsTemporalFrequency(n).paramsFile((end-5):end)];
%        end
%        setAllRetParams(paramsTemporalFrequency, combinedDT);
%        rmRunTemporalFreq1d(VOLUME{1},combinedDT, 'gray-Layer1',14,{'1g'},[],sprintf('retModel-%s-Lin-1dGaussianXnoY-TemporalFreq-%s',datestr(now,'yyyymmdd-HHMMSS'), num2str(20,2)),[],0);
%        load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsMix.mat')
       
   else
       load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsSessions.mat')
       setAllRetParams(paramsDurationPeriod, combinedDT);
       rmRunDuration2dOval(VOLUME{1},combinedDT, 'gray-Layer1',5,{'1g'},[],[],[],1,'free'); %This is the best fitting model. Setting the 4th input arguement to '5' fits the HRF parameters too.
%        setAllRetParams(paramsDurationPeriod, combinedDT);
%        rmRunDuration2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],maxCores,1);
%        setAllRetParams(paramsDurationPeriod, combinedDT);
%        rmRunLogDuration2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],1,'free');
%        for whichModel=[12 14 15 16]
%            setAllRetParams(paramsDurationPeriod, combinedDT);
%            rmRunDurFreq2d(VOLUME{1},combinedDT, 'gray-Layer1',whichModel,{'1g'},[],[],[],1);
%        end
%        
%        setAllRetParams(paramsOccupancyPeriod, combinedDT);
%        rmRunOccupancyPeriod2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],1,'free');
%        for whichModel=[12 14]
%            setAllRetParams(paramsOccupancyPeriod, combinedDT);
%            rmRunOccupancyFreq2d(VOLUME{1},combinedDT, 'gray-Layer1',whichModel,{'1g'},[],[],[],1);
%        end
%        
%        setAllRetParams(paramsOnTimeOffTime, combinedDT);
%        rmRunOnTimeOffTime2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],1,'free');
%        
%        setAllRetParams(paramsTemporalFrequency, combinedDT);
%        rmRunTemporalFreq1d(VOLUME{1},combinedDT, 'gray-Layer1',14,{'1g'},[],[],[],1);
%        
%        load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsSessions.mat')
%        for n=1:4
%            paramsTemporalFrequency(n).paramsFile=[paramsTemporalFrequency(n).paramsFile(1:(end-9)), 'TF', paramsTemporalFrequency(n).paramsFile((end-5):end)];
%        end
%        setAllRetParams(paramsTemporalFrequency, combinedDT);
%        rmRunTemporalFreq1d(VOLUME{1},combinedDT, 'gray-Layer1',14,{'1g'},[],sprintf('retModel-%s-Lin-1dGaussianXnoY-TemporalFreq-%s',datestr(now,'yyyymmdd-HHMMSS'), num2str(20,2)),[],1);
%        load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsSessions.mat')
   end
   close(1); mrvCleanWorkspace; 
end

%When searchFit or HRFfit goes wrong
whichDTs=1:3;
parpool(length(whichDTs))
parfor whichDT=whichDTs
    folderName=[pwd '/Gray/' dataTYPES(combinedDT(whichDT)).name]
    %S10
    %modelFiles='retModel-20210709-173216-Lin-2dOvalGaussian-DurationPeriod-DT0.5-maxValue-2-expIntensity-free-fFit.mat';
     %S11
    %modelFiles='retModel-20210709-141736-Lin-2dOvalGaussian-DurationPeriod-DT0.5-maxValue-2-expIntensity-free-fFit.mat';
     %S8
    modelFiles='retModel-20210709-141202-Lin-2dOvalGaussian-DurationPeriod-DT0.5-maxValue-2-expIntensity-free-fFit.mat';
%rmMainPostGrid([1 combinedDT(whichDT)],'gray-Layer1',5, [folderName '/' modelFiles]);
    rmMainPostSearch([1 combinedDT(whichDT)],'gray-Layer1',5, [folderName '/' modelFiles]);
end
delete(gcp('nocreate'))
close(1); mrvCleanWorkspace; 

%Small test ROI
rmRunDuration2dOval(VOLUME{1},combinedDT(whichDT), 'LeftTPOMap',5,{'1g'},[],[],[],0,'free');

whichSubs=[8];
%Redoing models with a searchfit problem (December 2nd)
maxCores=3;
for thisSub=1:length(whichSubs)
   cd(paths{whichSubs(thisSub)})
   mrVista 3;
   
   %The code to run for each
   if whichSubs(thisSub)<=6
       combinedDT=5:7;
       load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsMix.mat')
       setAllRetParams(paramsDurationPeriod, combinedDT);
       
       for whichModel=[12]
           setAllRetParams(paramsDurationPeriod, combinedDT);
           rmRunDurFreq2d(VOLUME{1},combinedDT, 'gray-Layer1',whichModel,{'1g'},[],[],[],0);
       end
       
       for whichModel=[12]
           setAllRetParams(paramsOccupancyPeriod, combinedDT);
           rmRunOccupancyFreq2d(VOLUME{1},combinedDT, 'gray-Layer1',whichModel,{'1g'},[],[],[],0);
       end

   else
       if whichSubs(thisSub)==7
           combinedDT=HarvDTs; %For Harv
       elseif whichSubs(thisSub)==8
           combinedDT=NellDTs;
       end
       load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsSessions.mat')
       
       for whichModel=[12]
           setAllRetParams(paramsDurationPeriod, combinedDT);
           rmRunDurFreq2d(VOLUME{1},combinedDT, 'gray-Layer1',whichModel,{'1g'},[],[],[],1);
       end
       
       for whichModel=[12]
           setAllRetParams(paramsOccupancyPeriod, combinedDT);
           rmRunOccupancyFreq2d(VOLUME{1},combinedDT, 'gray-Layer1',whichModel,{'1g'},[],[],[],1);
       end
       
   end
   close(1); mrvCleanWorkspace; 
end


%Fitting HRF (Only really needed for all data, free exponent)
whichSubs=[6 4 5]; %Done: 1 2 3 5 6
for thisSub=1:length(whichSubs);
    cd(paths{whichSubs(thisSub)})
    mrVista 3;
    if whichSubs(thisSub)<=6
        combinedDT=5:7;
    elseif whichSubs(thisSub)==7
        combinedDT=HarvDTs; %For Harv
    elseif whichSubs(thisSub)==8
        combinedDT=NellDTs;
    end
    whichDT=2;%:length(combinedDT)
    folderName=[pwd '/Gray/' dataTYPES(combinedDT(whichDT)).name, '/SearchFitFreeExponent'];
    modelFiles=dir([char(strcat(folderName, '/*-Lin-2dOvalGaussian-DurationPeriod*-free-*'))]);
    rmMainPostSearch([1 combinedDT(whichDT)],'gray-Layer1',5, [folderName '/' modelFiles.name]);
    close(1); mrvCleanWorkspace;
end
 
whichSubs=[6 4 5]; %Done: 1 2 3 5 6
for thisSub=1:length(whichSubs);
    cd(paths{whichSubs(thisSub)})
    mrVista 3;
    if whichSubs(thisSub)<=6
        combinedDT=5:7;
    elseif whichSubs(thisSub)==7
        combinedDT=HarvDTs; %For Harv
    elseif whichSubs(thisSub)==8
        combinedDT=NellDTs;
    end
    whichDT=3;%:length(combinedDT)
    folderName=[pwd '/Gray/' dataTYPES(combinedDT(whichDT)).name, '/SearchFitFreeExponent'];
    modelFiles=dir([char(strcat(folderName, '/*-Lin-2dOvalGaussian-DurationPeriod*-free-*'))]);
    rmMainPostSearch([1 combinedDT(whichDT)],'gray-Layer1',5, [folderName '/' modelFiles.name]);
    close(1); mrvCleanWorkspace;
end

%Fitting xval models with HRF params from all. Move models to HRF fit
%folder first. (THIS APPROACH CAUSES SOME PROBLEMS, SOMETHING IS WRONG)
whichSubs=[6]; %Done: 1 2 3 5 6
for thisSub=1:length(whichSubs);
    cd(paths{whichSubs(thisSub)})
    mrVista 3;
    if whichSubs(thisSub)<=6
        combinedDT=5:7;
        folderName=[pwd '/Gray/' dataTYPES(combinedDT(1)).name, '/HrfFitFreeExponent'];
        rmFile=dir([char(strcat(folderName, '/*-Lin-2dOvalGaussian-DurationPeriod*-free-*'))]);
        load([folderName '/' rmFile.name]);
        hrfParams{1}=params.stim(1).hrfType;
        hrfParams{2}=params.stim(1).hrfParams{2};
        load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsMix.mat')
        setAllRetParams(paramsDurationPeriod, combinedDT);
        rmRunDuration2dOval(VOLUME{1},combinedDT(2:3), 'gray-Layer1',4,{'1g'},hrfParams,sprintf('retModel-%s-2dDGaussian-lin-FullBlanks-DT0.5-HRFfromAll',datestr(now,'yyyymmdd-HHMMSS')),[],0,'free');
    else
        if whichSubs(thisSub)==7
            combinedDT=HarvDTs; %For Harv
        elseif whichSubs(thisSub)==8
            combinedDT=NellDTs;
        end
        folderName=[pwd '/Gray/' dataTYPES(combinedDT(1)).name, '/HrfFitFreeExponent'];
        rmFile=dir([char(strcat(folderName, '/*-Lin-2dOvalGaussian-DurationPeriod*-free-*'))]);
        load([folderName '/' rmFile.name]);
        hrfParams{1}=params.stim(1).hrfType;
        hrfParams{2}=params.stim(1).hrfParams{2};
        load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsSessions.mat')
        setAllRetParams(paramsDurationPeriod, combinedDT);
        rmRunDuration2dOval(VOLUME{1},combinedDT(2:3), 'gray-Layer1',4,{'1g'},hrfParams,sprintf('retModel-%s-2dDGaussian-lin-FullBlanks-DT0.5-HRFfromAll',datestr(now,'yyyymmdd-HHMMSS')),[],1,'free');

    end
    close(1); mrvCleanWorkspace;
end

%Cross validation
whichSubs=[7]; %To do: 5
for thisSub=whichSubs
   cd(paths{thisSub})
   mrVista 3;
   
   %The code to run for each
   if thisSub<=6
       combinedDT=5:7;
       load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsMix.mat')
       
   else
       if thisSub==7
           combinedDT=HarvDTs; %For Harv
       elseif thisSub==8
           combinedDT=NellDTs;
       end
       load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsSessions.mat')
   end
   allXvalDTs=combinedDT(2:3);
   
   %Linear Gaussian models
   for n=1:length(allXvalDTs)
       files=dir(['Gray/' dataTYPES(allXvalDTs(n)).name, '/SearchFit/', '*.mat']);
       thisPath=['Gray/' dataTYPES(allXvalDTs(n)).name, '/SearchFit/'];
       otherPath=['Gray/' dataTYPES(allXvalDTs(3-n)).name, '/SearchFit/'];
%        eval(['!mkdir ',  '"',otherPath, 'xval"']);
%        eval(['!mkdir ',  '"',otherPath, 'xvalRefit"']);
       
       for whichFile=1:length(files)
           eval(['!cp ', '"', thisPath, files(whichFile).name, '" "', otherPath, 'xval/xval-', files(whichFile).name, '"']);
       end
   end
   parpool(2)
   parfor whichDT=1:length(allXvalDTs)
       folderName=[pwd '/Gray/' dataTYPES(allXvalDTs(whichDT)).name, '/SearchFit/xval'];
       modelFiles=dir([folderName '/*Fit.mat']);
       for whichModel=1:length(modelFiles)
           rmMainPostSearch([1 allXvalDTs(whichDT)],'gray-Layer1',4, [folderName '/' modelFiles(whichModel).name]);
       end
   end
   delete(gcp('nocreate'))
   
   %Exponential Gaussian models
   for n=1:length(allXvalDTs)
       files=dir(['Gray/' dataTYPES(allXvalDTs(n)).name, '/SearchFitFreeExponent/', '*.mat']);
       thisPath=['Gray/' dataTYPES(allXvalDTs(n)).name, '/SearchFitFreeExponent/'];
       otherPath=['Gray/' dataTYPES(allXvalDTs(3-n)).name, '/SearchFitFreeExponent/'];
%        eval(['!mkdir ',  '"',otherPath, 'xval"']);
%        eval(['!mkdir ',  '"',otherPath, 'xvalRefit"']);
       
       for whichFile=1:length(files)
           eval(['!cp ', '"', thisPath, files(whichFile).name, '" "', otherPath, 'xval/xval-', files(whichFile).name, '"']);
       end
   end
   parpool(2)
   parfor whichDT=1:length(allXvalDTs)
       folderName=[pwd '/Gray/' dataTYPES(allXvalDTs(whichDT)).name, '/SearchFitFreeExponent/xval'];
       modelFiles=dir([folderName '/*.mat']);
       for whichModel=4%1:length(modelFiles)
           rmMainPostSearch([1 allXvalDTs(whichDT)],'gray-Layer1',4, [folderName '/' modelFiles(whichModel).name]);
       end
   end
   delete(gcp('nocreate'))
   
   %Monotonic component models
   for n=1:length(allXvalDTs)
       files=dir(['Gray/' dataTYPES(allXvalDTs(n)).name, '/Monotonic/', '*.mat']);
       thisPath=['Gray/' dataTYPES(allXvalDTs(n)).name, '/Monotonic/'];
       otherPath=['Gray/' dataTYPES(allXvalDTs(3-n)).name, '/Monotonic/'];
%        eval(['!mkdir ',  '"',otherPath, 'xval"']);
%        eval(['!mkdir ',  '"',otherPath, 'xvalRefit"']);
       
       for whichFile=1:length(files)
           eval(['!cp ', '"', thisPath, files(whichFile).name, '" "', otherPath, 'xval/xval-', files(whichFile).name, '"']);
       end
   end
   
   parpool(2)
   parfor whichDT=1:length(allXvalDTs)
       folderName=[pwd '/Gray/' dataTYPES(allXvalDTs(whichDT)).name, '/Monotonic/xval'];
       modelFiles=dir([folderName '/*.mat']);
       for whichModel=1:length(modelFiles)
           rmMainPostSearch([1 allXvalDTs(whichDT)],'gray-Layer1',4, [folderName '/' modelFiles(whichModel).name]);
       end
   end
   delete(gcp('nocreate'))
   
      %For a set of models in a folder called Tmp
      for n=1:length(allXvalDTs)
          files=dir(['Gray/' dataTYPES(allXvalDTs(n)).name, '/Tmp/', '*.mat']);
          thisPath=['Gray/' dataTYPES(allXvalDTs(n)).name, '/Tmp/'];
          otherPath=['Gray/' dataTYPES(allXvalDTs(3-n)).name, '/Tmp/'];
          %        eval(['!mkdir ',  '"',otherPath, 'xval"']);
          %        eval(['!mkdir ',  '"',otherPath, 'xvalRefit"']);
          
          for whichFile=1:length(files)
              eval(['!cp ', '"', thisPath, files(whichFile).name, '" "', otherPath, 'xval/xval-', files(whichFile).name, '"']);
          end
      end
      
      parpool(2)
      parfor whichDT=1:length(allXvalDTs)
          folderName=[pwd '/Gray/' dataTYPES(allXvalDTs(whichDT)).name, '/Tmp/xval'];
          modelFiles=dir([folderName '/*Fit.mat']);
          for whichModel=1:length(modelFiles)
              rmMainPostSearch([1 allXvalDTs(whichDT)],'gray-Layer1',4, [folderName '/' modelFiles(whichModel).name]);
          end
      end
      delete(gcp('nocreate'))

   close(1); mrvCleanWorkspace; 
end

%Flat model, frequency only. Edit rmRunOccupancyFreq2d at line 60 for name
%Need to edit rmGridFit at line 210 to only test flat RF
%NOTE: This model has no free parameters, so does not need cross
%validation. Copy the model directly into the xvalRefit folder in the same
%data type.
whichSubs=[6];

maxCores=3;
for thisSub=1:length(whichSubs)
    cd(paths{whichSubs(thisSub)})
    mrVista 3;
    
    %The code to run for each
    if whichSubs(thisSub)<=6
        combinedDT=5:7;
        load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsMix.mat')
        setAllRetParams(paramsOccupancyPeriod, combinedDT);
        rmRunOccupancyFreq2d(VOLUME{1},combinedDT, 'gray-Layer1',14,{'1g'},[],[],[],0);
    else
       if thisSub==7
           combinedDT=HarvDTs; %For Harv
       elseif thisSub==8
           combinedDT=NellDTs;
       end
       load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsSessions.mat')
       setAllRetParams(paramsOccupancyPeriod, combinedDT);
        rmRunOccupancyFreq2d(VOLUME{1},combinedDT, 'gray-Layer1',14,{'1g'},[],[],[],1);
    end
    close(1); mrvCleanWorkspace;
end

%In resulting models, to make TSeries plots
colStr{1}='k'
colStr{2}='b'
colStr{3}='r'
colStr{4}='g'
figure; hold on;
for n=1:4
    plot(M.currTsData.x((1):(56))+1.05, M.currTsData.raw(((n-1)*56+1):(n*56)), [colStr{n} 'o']);
    plot(M.currTsData.x((1):(56))+1.05, M.currTsData.pred(((n-1)*56+1):(n*56)), [colStr{n} '-']);
end
axis square
set(gca, 'XTick', [0 29.4 58.8 88.2 117.6])
axis([0 117.6 10*-.75 10])

%Make pRF image
[X,Y]=meshgrid([0:0.01:1, 1.98:0.01:2.02],[0:0.01:1 2.08:0.01:2.12]);
%M.rfParams=[0.5    .6    0.75         0    0.25    pi/2    0.1474]
RF=rfGaussian2d(X,Y,M.rfParams(3),M.rfParams(5),M.rfParams(6), M.rfParams(1),M.rfParams(2));
% if length(M.rfParams)>=7
%     freq=5./Y;
%     scale=1./(freq'.^M.rfParams(7));
%     scale=scale.*freq';
%     RF=5*RF./scale;
% end
figure; imagesc(flipud(RF./max(RF(:)))); %colormap gray
axis image
colorbar

Yvals=0:0.01:20;
exps=M.rfParams(7);
figure; hold on;
for count=1:length(exps)
    plot(Yvals, ((Yvals.^exps(count))));
end
axis square
axis([0 20 0 10])
set(gca, 'XTick', 0:2:20)

%Make 1D pRF pass projections
X=0.05:0.05:1;
Y=0.05:0.05:1;
RF=rfGaussian2d(X,Y,M.rfParams(3),M.rfParams(5),M.rfParams(6), M.rfParams(1),M.rfParams(2));
if length(M.rfParams)>=7
    freq=5./Y;
    scale=1./(freq'.^M.rfParams(7));
    scale=scale.*freq';
    RF=5*RF./scale;
end
RF=[RF; RF];
RF=RF(:);
Xaxis=[X;X];
Xaxis=Xaxis(:);
Xaxis=[0; Xaxis(1:(end-1))];
figure; plot(Xaxis,RF, colStr{1})
RF=rfGaussian2d(0.05.*ones(size(X)),Y,M.rfParams(3),M.rfParams(5),M.rfParams(6), M.rfParams(1),M.rfParams(2));
if length(M.rfParams)>=7
    freq=5./Y;
    scale=1./(freq'.^M.rfParams(7));
    scale=scale.*freq';
    RF=5*RF./scale;
end
RF=[RF; RF];
RF=RF(:);
hold on; plot(Xaxis,RF, colStr{2})
RF=rfGaussian2d(X,ones(size(Y)),M.rfParams(3),M.rfParams(5),M.rfParams(6), M.rfParams(1),M.rfParams(2));
if length(M.rfParams)>=7
    freq=5./ones(size(Y));
    scale=1./(freq'.^M.rfParams(7));
    scale=scale.*freq';
    RF=5*RF./scale;
end
RF=[RF; RF];
RF=RF(:);
hold on; plot(Xaxis,RF, colStr{3})
X=[0.05:0.05:0.5, 0.05:0.05:0.5];
Y=[0.95:-0.05:0.5, 0.55:0.05:1];
RF=rfGaussian2d(X,Y,M.rfParams(3),M.rfParams(5),M.rfParams(6), M.rfParams(1),M.rfParams(2));
if length(M.rfParams)>=7
    freq=5./Y;
    scale=1./(freq'.^M.rfParams(7));
    scale=scale.*freq';
    RF=5*RF./scale;
end
RF=[RF; RF];
RF=RF(:);
hold on; plot(Xaxis,RF, colStr{4})

%Convert logOP RF into linDP space
[X,Y]=meshgrid([0:0.01:1 1.98:0.01:2.02],[0:0.01:1 2.08:0.01:2.12]);
X=log(X./Y)+4;
Y=log(Y)+4;
RF=rfGaussian2d(X,Y,MLogOP.rfParams(3),MLogOP.rfParams(5),MLogOP.rfParams(6), MLogOP.rfParams(1),MLogOP.rfParams(2));
figure; imagesc(flipud(RF))
[X,Y]=meshgrid([0:0.01:1 1.98:0.01:2.02],[0:0.01:1 2.08:0.01:2.12]);
RF=rfGaussian2d(X,Y,MLinDP.rfParams(3),MLinDP.rfParams(5),MLinDP.rfParams(6), MLinDP.rfParams(1),MLinDP.rfParams(2));
figure; imagesc(flipud(RF))

[X,Y]=meshgrid([0.05:0.05:1 2],[0.05:0.05:1 2.1]);
Xbin=zeros(size(X));
Xbin(X==0.05)=1;
Xbin(Y==1)=1;
Xbin(X==Y)=1;
figure; imagesc(flipud(Xbin));
Xs=X(Xbin==1);
Ys=Y(Xbin==1);
Os=Xs./Ys;
figure; plot(Os, Ys, '.')

X=[0.05:0.05:1 2];
Ylog=log(X);
Ylog=Ylog-min(Ylog);
Ylog=Ylog./max(Ylog);
Yexp=X.^0.2;
Yexp=Yexp-min(Yexp);
Yexp=Yexp./max(Yexp);

figure; plot(X, Ylog);
hold on; plot(X, Yexp);

%Convert compressive frequency parameters into parameters in real space
rfParams=[0.1734    0.3014    1.7317         0    0.3115    1.1579    0.1474]

Xgrid=-4:0.01:4;
Ygrid=0:0.01:2;
[Xgrid, Ygrid]=meshgrid(Xgrid, Ygrid);
[maxX, maxY, longAxis, shortAxis, longAxisOrientation] = quantifyCompressiveRF(Xgrid, Ygrid, rfParams(1), rfParams(2), rfParams(3), rfParams(5), rfParams(6), rfParams(7));



%Make stim timing plots
load('params_ConLumOff2D.mat')
figure; stem(stimulus.seq(253:1092)~=85, 'k-'); axis square; axis([0 840 0.1 0.9])
load('params_ConDurOff2D.mat')
figure; stem(stimulus.seq(253:1092)~=85, 'k-'); axis square; axis([0 840 0.1 0.9])   
load('params_ConPerOff2D.mat')
figure; stem(stimulus.seq(253:1092)~=85, 'k-'); axis square; axis([0 840 0.1 0.9])
load('params_TimeGapsOn2D.mat')
figure; stem(stimulus.seq([253:672 799:1218])~=85, 'k-'); axis square; axis([0 840 0.1 0.9])
hold on; plot([420 420], [0 1])


%Get drawn ROIs
%Line ROI
meshROI2Volume([], 2); VOLUME{1}.ROIs(VOLUME{1}.selectedROI) = roiSortLineCoords(VOLUME{1});VOLUME{1}=editROIFields(VOLUME{1}); setROIPopup(VOLUME{1});
%Filled ROI
meshROI2Volume([], 2); VOLUME{1}=editROIFields(VOLUME{1}); setROIPopup(VOLUME{1});



%Cross validation
basePath=pwd;
%Determine which data types to work on
allXvalDTs=combinedDT(2:3);

%CombineModels and copy to cross-validation folders
for n=1:length(allXvalDTs)
    files=dir(['Gray/' dataTYPES(allXvalDTs(n)).name, '/SearchFit/', '*.mat']);
    thisPath=['Gray/' dataTYPES(allXvalDTs(n)).name, '/SearchFit/'];
    otherPath=['Gray/' dataTYPES(allXvalDTs(3-n)).name, '/SearchFit/'];
    eval(['!mkdir ',  '"',otherPath, 'xval"']);
    eval(['!mkdir ',  '"',otherPath, 'xvalRefit"']);
    for whichFile=1:length(files)
        eval(['!cp ', '"', thisPath, files(whichFile).name, '" "', otherPath, 'xval/xval-', files(whichFile).name, '"']);
    end
end

%Determine fits of these models in cross-validated data
for whichDT=1:length(allXvalDTs)
    folderName=[pwd '/Gray/' dataTYPES(allXvalDTs(whichDT)).name, '/SearchFit/xval'];
    modelFiles=dir([folderName '/*-fFit.mat']);
    if length(modelFiles)>1 && length(modelFiles)<5
        parpool(length(modelFiles))
    elseif length(modelFiles)>1
        parpool(4)
    end
    parfor whichModel=1:length(modelFiles)
        rmMainPostSearch([1 allXvalDTs(whichDT)],'gray-Layer1',4, [folderName '/' modelFiles(whichModel).name]);
    end
    delete(gcp('nocreate'))
end

%Cross validation for log intensity models
paths{1}='/media/Storage3/WholeBrainNewSeg/Subject 8/S8_Combined';
paths{2}='/media/Storage3/WholeBrainNewSeg/Subject 9/S9_Combined';
paths{3}='/media/Storage3/WholeBrainNewSeg/Subject 10/S10_Combined';
paths{4}='/media/Storage3/WholeBrainNewSeg/Subject 11/S11_Combined';
paths{5}='/media/Storage3/WholeBrainNewSeg/Subject 12/S12_Combined';
paths{6}='/media/Storage3/WholeBrainNewSeg/Subject 13/S13_Combined';
paths{7}='/media/Storage3/WholeBrainNewSeg/Harv/HarvCombined';
paths{8}='/media/Storage3/WholeBrainNewSeg/Nell/NellCombined';

nCores=4;

%Redoing searchfits to get correct rss estimate
whichSubs=[7:8];
for thisSub=whichSubs
   cd(paths{thisSub})
   mrVista 3;
   
   %The code to run for each
   if thisSub<=6
       combinedDT=9:11;
       load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsMix.mat')
       
   else
       if thisSub==7
           combinedDT=[67 69 70]; %For Harv
       elseif thisSub==8
           combinedDT=[72 84 85];
       end
       load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsSessions.mat')
   end
   for whichDT=1:length(combinedDT)
       folderName=[pwd '/Gray/' dataTYPES(combinedDT(whichDT)).name, '/HrfFitFreeExponent'];
       modelFiles=dir([folderName '/*-fFit.mat']);
       if length(modelFiles)>1 && length(modelFiles)<5
           parpool(length(modelFiles))
       elseif length(modelFiles)>1
           parpool(nCores)
       end
       thisDT=combinedDT(whichDT);
       parfor whichModel=1:length(modelFiles)
           rmMainPostSearch([1 thisDT],'gray-Layer1',4, [folderName '/' modelFiles(whichModel).name]);
       end
       delete(gcp('nocreate'))
   end
   
   close(1); mrvCleanWorkspace; 
end

   
whichSubs=[7:8];
for thisSub=whichSubs
   cd(paths{thisSub})
   mrVista 3;
   
   %The code to run for each
   if thisSub<=6
       combinedDT=9:11;
       load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsMix.mat')
       
   else
       if thisSub==7
           combinedDT=[67 69 70]; %For Harv
       elseif thisSub==8
           combinedDT=[72 84 85];
       end
       load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsSessions.mat')
   end
   allXvalDTs=combinedDT(2:3);
   
   for n=1:length(allXvalDTs)
       files=dir(['Gray/' dataTYPES(allXvalDTs(n)).name, '/SearchFit/', '*.mat']);
       thisPath=['Gray/' dataTYPES(allXvalDTs(n)).name, '/SearchFit/'];
       otherPath=['Gray/' dataTYPES(allXvalDTs(3-n)).name, '/SearchFit/'];
       eval(['!mkdir ',  '"',otherPath, 'xval"']);
       eval(['!mkdir ',  '"',otherPath, 'xvalRefit"']);
       
       for whichFile=1:length(files)
           eval(['!cp ', '"', thisPath, files(whichFile).name, '" "', otherPath, 'xval/xval-', files(whichFile).name, '"']);
       end
   end
   for whichDT=1:length(allXvalDTs)
       folderName=[pwd '/Gray/' dataTYPES(allXvalDTs(whichDT)).name, '/SearchFit/xval'];
       modelFiles=dir([folderName '/*-fFit.mat']);
       if length(modelFiles)>1 && length(modelFiles)<5
           parpool(length(modelFiles))
       elseif length(modelFiles)>1
           parpool(nCores)
       end
       parfor whichModel=1:length(modelFiles)
           rmMainPostSearch([1 allXvalDTs(whichDT)],'gray-Layer1',4, [folderName '/' modelFiles(whichModel).name]);
       end
       delete(gcp('nocreate'))
   end
   
   close(1); mrvCleanWorkspace; 
end

whichSubs=[1:8];
for thisSub=whichSubs
   cd(paths{thisSub})
   mrVista 3;
   
   %The code to run for each
   if thisSub<=6
       combinedDT=9:11;
       load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsMix.mat')
       
   else
       if thisSub==7
           combinedDT=[67 69 70]; %For Harv
       elseif thisSub==8
           combinedDT=[72 84 85];
       end
       load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsSessions.mat')
   end
   allXvalDTs=combinedDT(2:3);
   
   for n=1:length(allXvalDTs)
       files=dir(['Gray/' dataTYPES(allXvalDTs(n)).name, '/SearchFitLogIntensity/', '*.mat']);
       thisPath=['Gray/' dataTYPES(allXvalDTs(n)).name, '/SearchFitLogIntensity/'];
       otherPath=['Gray/' dataTYPES(allXvalDTs(3-n)).name, '/SearchFitLogIntensity/'];
       eval(['!mkdir ',  '"',otherPath, 'xval"']);
       eval(['!mkdir ',  '"',otherPath, 'xvalRefit"']);
       
       for whichFile=1:length(files)
           eval(['!cp ', '"', thisPath, files(whichFile).name, '" "', otherPath, 'xval/xval-', files(whichFile).name, '"']);
       end
   end
   for whichDT=1:length(allXvalDTs)
       folderName=[pwd '/Gray/' dataTYPES(allXvalDTs(whichDT)).name, '/SearchFitLogIntensity/xval'];
       modelFiles=dir([folderName '/*-fFit.mat']);
       if length(modelFiles)>1 && length(modelFiles)<5
           parpool(length(modelFiles))
       elseif length(modelFiles)>1
           parpool(nCores)
       end
       parfor whichModel=1:length(modelFiles)
           rmMainPostSearch([1 allXvalDTs(whichDT)],'gray-Layer1',4, [folderName '/' modelFiles(whichModel).name]);
       end
       delete(gcp('nocreate'))
   end
   
   close(1); mrvCleanWorkspace; 
end

whichSubs=[8];
for thisSub=whichSubs
   cd(paths{thisSub})
   mrVista 3;
   
   %The code to run for each
   if thisSub<=6
       combinedDT=9:11;
       load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsMix.mat')
       
   else
       if thisSub==7
           combinedDT=[67 69 70]; %For Harv
       elseif thisSub==8
           combinedDT=[72 84 85];
       end
       load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsSessions.mat')
   end
   allXvalDTs=combinedDT(2:3);
   
   for n=1:length(allXvalDTs)
       files=dir(['Gray/' dataTYPES(allXvalDTs(n)).name, '/SearchFitFreeExponent/', '*.mat']);
       thisPath=['Gray/' dataTYPES(allXvalDTs(n)).name, '/SearchFitFreeExponent/'];
       otherPath=['Gray/' dataTYPES(allXvalDTs(3-n)).name, '/SearchFitFreeExponent/'];
       eval(['!mkdir ',  '"',otherPath, 'xval"']);
       eval(['!mkdir ',  '"',otherPath, 'xvalRefit"']);
       
       for whichFile=1:length(files)
           eval(['!cp ', '"', thisPath, files(whichFile).name, '" "', otherPath, 'xval/xval-', files(whichFile).name, '"']);
       end
   end
   for whichDT=1:length(allXvalDTs)
       folderName=[pwd '/Gray/' dataTYPES(allXvalDTs(whichDT)).name, '/SearchFitFreeExponent/xval'];
       modelFiles=dir([folderName '/*-fFit.mat']);
       if length(modelFiles)>1 && length(modelFiles)<5
           parpool(length(modelFiles))
       elseif length(modelFiles)>1
           parpool(nCores)
       end
       parfor whichModel=1:length(modelFiles)
           rmMainPostSearch([1 allXvalDTs(whichDT)],'gray-Layer1',4, [folderName '/' modelFiles(whichModel).name]);
       end
       delete(gcp('nocreate'))
   end
   
   close(1); mrvCleanWorkspace; 
end

whichSubs=[7:8];
for thisSub=whichSubs
   cd(paths{thisSub})
   mrVista 3;
   
   %The code to run for each
   if thisSub<=6
       combinedDT=9:11;
       load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsMix.mat')
       
   else
       if thisSub==7
           combinedDT=[67 69 70]; %For Harv
       elseif thisSub==8
           combinedDT=[72 84 85];
       end
       load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsSessions.mat')
   end
   allXvalDTs=combinedDT(2:3);
   
   for n=1:length(allXvalDTs)
       files=dir(['Gray/' dataTYPES(allXvalDTs(n)).name, '/HrfFitFreeExponent/', '*.mat']);
       thisPath=['Gray/' dataTYPES(allXvalDTs(n)).name, '/HrfFitFreeExponent/'];
       otherPath=['Gray/' dataTYPES(allXvalDTs(3-n)).name, '/HrfFitFreeExponent/'];
       eval(['!mkdir ',  '"',otherPath, 'xval"']);
       eval(['!mkdir ',  '"',otherPath, 'xvalRefit"']);
       
       for whichFile=1:length(files)
           eval(['!cp ', '"', thisPath, files(whichFile).name, '" "', otherPath, 'xval/xval-', files(whichFile).name, '"']);
       end
   end
   for whichDT=1:length(allXvalDTs)
       folderName=[pwd '/Gray/' dataTYPES(allXvalDTs(whichDT)).name, '/HrfFitFreeExponent/xval'];
       modelFiles=dir([folderName '/*-fFit.mat']);
       if length(modelFiles)>1 && length(modelFiles)<5
           parpool(length(modelFiles))
       elseif length(modelFiles)>1
           parpool(nCores)
       end
       parfor whichModel=1:length(modelFiles)
           rmMainPostSearch([1 allXvalDTs(whichDT)],'gray-Layer1',4, [folderName '/' modelFiles(whichModel).name]);
       end
       delete(gcp('nocreate'))
   end
   
   close(1); mrvCleanWorkspace; 
end

whichSubs=[1:8];
for thisSub=whichSubs
   cd(paths{thisSub})
   mrVista 3;
   
   %The code to run for each
   if thisSub<=6
       combinedDT=9:11;
       load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsMix.mat')
       
   else
       if thisSub==7
           combinedDT=[67 69 70]; %For Harv
       elseif thisSub==8
           combinedDT=[72 84 85];
       end
       load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsSessions.mat')
   end
   allXvalDTs=combinedDT(2:3);
   
   for n=1:length(allXvalDTs)
       files=dir(['Gray/' dataTYPES(allXvalDTs(n)).name, '/Monotonic/', '*TemporalFreqPerEvent*']);
       thisPath=['Gray/' dataTYPES(allXvalDTs(n)).name, '/Monotonic/'];
       otherPath=['Gray/' dataTYPES(allXvalDTs(3-n)).name, '/Monotonic/'];
%        eval(['!mkdir ',  '"',otherPath, 'xval"']);
%        eval(['!mkdir ',  '"',otherPath, 'xvalRefit"']);
       
       for whichFile=1:length(files)
           eval(['!cp ', '"', thisPath, files(whichFile).name, '" "', otherPath, 'xval/xval-', files(whichFile).name, '"']);
       end
   end
   
   for whichDT=1:length(allXvalDTs)
       folderName=[pwd '/Gray/' dataTYPES(allXvalDTs(whichDT)).name, '/Monotonic/xval'];
       modelFiles=dir([folderName '/*TemporalFreqPerEvent*']);
       if length(modelFiles)>1 && length(modelFiles)<5
           parpool(length(modelFiles))
       elseif length(modelFiles)>1
           parpool(nCores)
       end
       for whichModel=1:length(modelFiles)
           rmMainPostSearch([1 allXvalDTs(whichDT)],'gray-Layer1',4, [folderName '/' modelFiles(whichModel).name]);
       end
       delete(gcp('nocreate'))
   end
   %delete(gcp('nocreate'))
   close(1); mrvCleanWorkspace; 
end


%Cross validation for a sub-set of models
basePath=pwd;
%Determine which data types to work on
allXvalDTs=combinedDT(2:3);

%CombineModels and copy to cross-validation folders
for n=1:length(allXvalDTs)
    files=dir(['Gray/' dataTYPES(allXvalDTs(n)).name, '/SearchFit/', '*DurationPeriod*.mat']);
    thisPath=['Gray/' dataTYPES(allXvalDTs(n)).name, '/SearchFit/'];
    otherPath=['Gray/' dataTYPES(allXvalDTs(3-n)).name, '/SearchFit/'];
    
    for whichFile=1:length(files)
        eval(['!cp ', '"', thisPath, files(whichFile).name, '" "', otherPath, 'xval/xval-', files(whichFile).name, '"']);
    end
end

%Determine fits of these models in cross-validated data
for whichDT=1:length(allXvalDTs)
    folderName=[pwd '/Gray/' dataTYPES(allXvalDTs(whichDT)).name, '/SearchFit/xval'];
    modelFiles=dir([folderName '/*DurationPeriod*-fFit.mat']);
    if length(modelFiles)>1 && length(modelFiles)<5
        parpool(length(modelFiles))
    elseif length(modelFiles)>1
        parpool(4)
    end
    parfor whichModel=1:length(modelFiles)
        rmMainPostSearch([1 allXvalDTs(whichDT)],'gray-Layer1',4, [folderName '/' modelFiles(whichModel).name]);
    end
    delete(gcp('nocreate'))
end

%Get ROI data from models
%For AFNI preproc models, distances from first model only.
%For revision, more limited set of models, but more easily defined
mapNames=["TLO", "TTOP", "TTOA", "TPO", "TLS", "TPCI", "TPCM", "TPCS", "TFI" "TFS"];
%mapNames=["TPCI"];
%Folder names to ge models from
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

paths{1}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/S08_Combined';
paths{2}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/S09_Combined';
paths{3}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/S10_Combined';
paths{4}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/S11_Combined';
paths{5}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/S12_Combined';
paths{6}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/S13_Combined_Seg20';
paths{7}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/HarvCombined';
paths{8}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/NellCombined';
subData=["dataS8", "dataS9", "dataS10", "dataS11", "dataS12", "dataS13", "dataHarv", "dataNell"];

whichSubs=[1:8]
for thisSub=[7];%6:length(whichSubs)
    dataTmp=getTimingModelData(paths{whichSubs(thisSub)}, modelFolders, modelNames, modelFileNames, mapNames, [5:7]);
    %eval([char(subData(whichSubs(thisSub))), '=dataTmp;'])
end

dataS8=rmfield(dataS8, 'All');
dataS9=rmfield(dataS9, 'All');
dataS10=rmfield(dataS10, 'All');
dataS11=rmfield(dataS11, 'All');
dataS12=rmfield(dataS12, 'All');
dataS13=rmfield(dataS13, 'All');
dataHarv=rmfield(dataHarv, 'All');
dataNell=rmfield(dataNell, 'All');

%Join all ROIs within subject
thisData='dataS8';
UniteROIsRevision;
thisData='dataS9';
UniteROIsRevision;
thisData='dataS10';
UniteROIsRevision;
thisData='dataS11';
UniteROIsRevision;
thisData='dataS12';
UniteROIsRevision;
thisData='dataS13';
UniteROIsRevision;
thisData='dataHarv';
UniteROIsRevision;
thisData='dataNell';
UniteROIsRevision;

%The above code is also in the function 'UniteROIstruct.m'. This can be
%scripted for all subjects, as follows
clear('dataAllSub'); %First run through only
thisData='dataS8';
UniteROIsubsRevision;
thisData='dataS9';
UniteROIsubsRevision;
thisData='dataS10';
UniteROIsubsRevision;
thisData='dataS11';
UniteROIsubsRevision;
thisData='dataS12';
UniteROIsubsRevision;
thisData='dataS13';
UniteROIsubsRevision;
thisData='dataHarv';
UniteROIsubsRevision;
thisData='dataNell';
UniteROIsubsRevision;

%New model comparisons
subjectOrder={'dataS8', 'dataS9', 'dataS10', 'dataS11', 'dataS12', 'dataS13', 'dataHarv', 'dataNell'};
%subjectOrder={'dataAllSub'};
mapNames={'TLO', 'TTOP', 'TTOA', 'TPO', 'TLS', 'TPCI', 'TPCM', 'TPCS', 'TFI' 'TFS'};%Revision
%mapNames={'All'};
hemispheres={'Left', 'Right'};
modelNames={'TemporalFreq', 'Flat', 'LinDF', 'LinDCompF', 'TemporalFreqPerEvent','GausOLinF', 'GausDLinF', 'GausOCompF', 'GausDCompF', 'DurationPeriodNoExponent', 'DurationPeriod', 'OccupancyPeriod', 'OnTimeOffTime', 'LogDurationPeriod'};


clear('minTvals')
clear('tMatrix')
clear('sterrs')        
clear('barData')
clear('barDataSub')
%Get data locations
for whichSub=1:length(subjectOrder)
    for whichMap=1:length(mapNames)
        for whichHemi=1:2;

        for modelN=1:length(modelNames)
            targetDataOdd{whichSub, whichMap, whichHemi, modelN}=char(strcat(subjectOrder{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNames{modelN}, '.Odd.vesXval'));
            targetDataEven{whichSub, whichMap,whichHemi, modelN}=char(strcat(subjectOrder{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNames{modelN}, '.Even.vesXval'));
        end
        end
    end
end

%Get data from locations
barDataMeans=nan([length(subjectOrder), length(mapNames), 2,2,length(modelNames)]);
for whichSub=1:length(subjectOrder)
    for whichMap=1:length(mapNames)
        if isfield(eval(subjectOrder{whichSub}), char(mapNames{whichMap}))
            %mapOK(whichSub,whichMap)=1;
            for whichHemi=1:2;
                if isfield(eval([char(subjectOrder{whichSub}), '.', char(mapNames{whichMap})]), char(hemispheres{whichHemi}))
                    for modelN=1:length(modelNames);
                        barData{whichSub, whichMap, whichHemi, 1}(:,modelN)=eval(targetDataOdd{whichSub, whichMap, whichHemi, modelN});
                        barData{whichSub, whichMap, whichHemi, 2}(:,modelN)=eval(targetDataEven{whichSub, whichMap, whichHemi, modelN});
                    end
                    for oddEven=1:2
                        indices=max(barData{whichSub, whichMap, whichHemi, oddEven},[],2)>0.2;
                        barDataMeans(whichSub, whichMap, whichHemi, oddEven,:)=mean(barData{whichSub, whichMap, whichHemi, oddEven}(indices,:), 1);
                    end
                end
            end
%         else
%             mapOK(whichSub,whichMap)=0;
        end
    end
end

barPoints=[];
whichBars=1:length(modelNames);
for n=1:length(whichBars)
    tmp=barDataMeans(:, :, :, :,whichBars(n));
    tmp=(tmp(:));
    tmp=tmp(~isnan(tmp));
    barMeans(n)=mean(tmp);
    barStd(n)=std(tmp);
    barSerr(n)=std(tmp)/sqrt(length(tmp));
    barPoints(:,n)=tmp;
    CI95(n,:) = tinv([0.025 0.975], length(tmp)-1);
    CI95(n,:) =bsxfun(@times, barSerr(n), CI95(n,:));
end
figure; subplot(1,2,1);
bar(barMeans);
hold on; errorbar(1:size(barMeans,2), barMeans, CI95(:,1), CI95(:,2), 'k.')
axis square;
xlim([0.5 length(barMeans)+0.5])
ylim([0 0.35])

for x=1:length(whichBars)
    for y=1:length(whichBars)
        [~,p,ci,stats] = ttest(barPoints(:,x), barPoints(:,y), 'tail', 'both');
        pvals(x,y)=p;
        tMatrix(x,y)=stats.tstat;
    end
end
tMatrix(tMatrix>200)=200;
tMatrix(tMatrix<-200)=-200;

tMatrixImg=zeros(size(tMatrix,1)*10, size(tMatrix,2)*10);
for x=1:size(tMatrix,1)
    for y=1:size(tMatrix,1)
        tMatrixImg((x-1)*10+1:x*10, (y-1)*10+1:y*10)=tMatrix(x,y);
    end
end
tMatrixImg(isnan(tMatrixImg))=0;
figure; subplot(1,2,2);
imagesc(tMatrixImg);


%Histograms of duration and period preferences
%subjectOrder=["dataS8", "dataS9", "dataS10", "dataS11", "dataS12", "dataS13", "dataHarv", "dataNell"];
subjectOrder={'dataAllSub'};
mapNames={'TLO', 'TTOP', 'TTOA', 'TPO', 'TLS', 'TPCI', 'TPCM', 'TPCS', 'TFI' 'TFS'};%Revision
hemispheres={'Left', 'Right'};

%Get data locations
veThresh=0.2;
whichSub=1;
figure; hold on;
for whichMap=1:length(mapNames)
    for whichHemi=1:2;
            data=eval(strcat(subjectOrder{whichSub}, '.', mapNames{whichMap}, '.', hemispheres{whichHemi}, '.DurationPeriodHRF.All'));
            subplot(length(hemispheres),length(mapNames),(whichHemi-1)*length(mapNames)+whichMap);
            histogram(data.y0s(data.y0s>0.05 & data.y0s<1 & data.ves>veThresh), 0.05:0.05:1, 'FaceColor', 'none');
            axis square;
            title(strcat(hemispheres(whichHemi), ' ', mapNames(whichMap)))
    end
end


%Timing vs distance etc (Revision process)
subData={'dataS8', 'dataS9', 'dataS10', 'dataS11', 'dataS12', 'dataS13', 'dataHarv', 'dataNell', 'dataAllSub'};
mapNames={'TLO', 'TTOP', 'TTOA', 'TPO', 'TLS', 'TPCI', 'TPCM', 'TPCS', 'TFI' 'TFS'};%Revision
%mapNames={'TPCI'}
hemispheres={'Left', 'Right'};

whichSubs=[7];%[7 8 2 3 4 6 ];
mapList=[1:length(mapNames)];
DTnames={'All', 'Odd', 'Even'};


for thisSub=whichSubs
    clear plotdata;
    eval(['data=', char(subData(thisSub)), ';'])
    for whichMap=mapList %1:length(mapNames)
        for whichHemi=1:2
            count=(whichMap-1)*2+whichHemi;
            dataPresent=char(strcat('isfield(data, ''', mapNames{whichMap}, ''') && isfield(data.',  mapNames{whichMap}, ',''',hemispheres{whichHemi}, ''')' ));
            if eval(dataPresent)
                for whichDT=1:3;
                    dataName=char(strcat('data.', mapNames{whichMap}, '.', hemispheres{whichHemi},'.DurationPeriod.', DTnames(whichDT)));
                    eval(['dataTmp=', dataName, ';'])
                    dataTmp.ROItitle=char(strcat(subData{thisSub},mapNames{whichMap}, hemispheres{whichHemi}));
                    
                    if thisSub==9
                        dataTmp.meanDist=25;
                    end
                    
                    plotdata{count, whichDT}=dataTmp(1,1);
                end
            end
        end
    end
    %try
    stats{thisSub}=PlotRoiDistanceAllRoiDurationPeriod(plotdata, 0.1, [0 1 0 0 0 0 0]); %Plot everything [1 1 1 0 1 1 1] %Tuning widths only [0 0 0 0 0 1 1]
    %end
end

subsTmp=whichSubs;
for whichSubs=subsTmp
    for map=mapList
        index=(map-1)*2+[1:2];
            stats{whichSubs}.pXprog(index,:)=statsTmp{whichSubs}.pXprog(index,:);
            stats{whichSubs}.pYprog(index,:)=statsTmp{whichSubs}.pYprog(index,:);
            stats{whichSubs}.rs(index)=statsTmp{whichSubs}.rs(index);
            stats{whichSubs}.ns(index)=statsTmp{whichSubs}.ns(index);
            stats{whichSubs}.rsxOddEven(index)=statsTmp{whichSubs}.rsxOddEven(index);
            stats{whichSubs}.rsyOddEven(index)=statsTmp{whichSubs}.rsyOddEven(index);
            stats{whichSubs}.nsOddEven(index)=statsTmp{whichSubs}.nsOddEven(index);
        try
            stats{whichSubs}.pSMajProg(index,:)=statsTmp{whichSubs}.pSMajProg(index,:);
            stats{whichSubs}.pSMinProg(index,:)=statsTmp{whichSubs}.pSMinProg(index,:);
        catch
            stats{whichSubs}=rmfield(stats{whichSubs}, 'pSMajProg');
            stats{whichSubs}=rmfield(stats{whichSubs}, 'pSMinProg');
            stats{whichSubs}.pSMajProg(index,:)=statsTmp{whichSubs}.pSMajProg(index,:);
            stats{whichSubs}.pSMinProg(index,:)=statsTmp{whichSubs}.pSMinProg(index,:);
        end
        stats{whichSubs}.tuningSlopesMajor(index,:)=statsTmp{whichSubs}.tuningSlopesMajor(index,:);
        stats{whichSubs}.tuningSlopesMinor(index,:)=statsTmp{whichSubs}.tuningSlopesMinor(index,:);
        stats{whichSubs}.tuningXs(index,:)=statsTmp{whichSubs}.tuningXs(index,:);
        stats{whichSubs}.tuningYmajor(index,:)=statsTmp{whichSubs}.tuningYmajor(index,:);
        stats{whichSubs}.tuningYminor(index,:)=statsTmp{whichSubs}.tuningYminor(index,:);
    end
end
whichSubs=subsTmp;

whichmaps=[3 9 10 11];
clear plotdatatmp
for n=1:length(whichmaps)
    for m=1:3
        plotdatatmp{n,m}=plotdata{whichmaps(n), m};
    end
end

%ANOVAs
subOrder={'dataS8', 'dataS9', 'dataS10', 'dataS11', 'dataS12', 'dataS13', 'dataHarv', 'dataNell'};
mapNames={'TLO', 'TTOP', 'TTOA', 'TPO', 'TLS', 'TPCI', 'TPCM', 'TPCS', 'TFI' 'TFS'};
hemispheres={'Left', 'Right'};

VEs=nan([length(subjectOrder) length(mapNames) length(hemispheres)]);
Exps=VEs;
SigmaMajor=VEs;
SigmaMinor=VEs;
SigmaRatio=VEs;
SigmaTheta=VEs;
Q1D=VEs;
Q2D=VEs;
Q3D=VEs;
IQRD=VEs;
Q1P=VEs;
Q2P=VEs;
Q3P=VEs;
IQRP=VEs;
linearSlopeMajor=VEs;
vSlopeMajor=VEs;
linearSlopeMinor=VEs;
vSlopeMinor=VEs;
nVoxels=VEs;
%figure; hold on;
for n=1:length(subjectOrder)
    for m=1:length(mapNames)
        for whichHemi=1:length(hemispheres)
            if eval(char(strcat('isfield(', subjectOrder{n}, ', ''', mapNames{m}, ''') && isfield(', subjectOrder{n}, '.', mapNames{m}, ',''',hemispheres{whichHemi}, ''')' ))) 
                
                eval(char(strcat('data=', subjectOrder{n},'.', mapNames{m},'.', hemispheres{whichHemi},'.DurationPeriodHRF.All',';')));
                
                veIndices=data.ves>0.1 & data.x0s>0.05 & data.x0s<1 ;
                if any(veIndices)
                    nVoxels(n,m,whichHemi)=sum(veIndices);
                end
                VEs(n,m,whichHemi)=mean(data.ves(veIndices));
                %Exps(n,m,whichHemi)=mean(data.exp(veIndices));
                SigmaMajor(n,m,whichHemi)=mean(data.sigmas(veIndices));
                SigmaMinor(n,m,whichHemi)=mean(data.sigmaMinor(veIndices));
                SigmaRatio(n,m,whichHemi)=mean(data.sigmas(veIndices)./data.sigmaMinor(veIndices));
                SigmaTheta(n,m,whichHemi)=mean(data.sigmaTheta(veIndices));
                Q1D(n,m,whichHemi)=prctile(data.x0s(veIndices), 25);
                Q2D(n,m,whichHemi)=mean(data.x0s(veIndices));
                Q3D(n,m,whichHemi)=prctile(data.x0s(veIndices), 75);
                IQRD(n,m,whichHemi)=prctile(data.x0s(veIndices), 75)-prctile(data.x0s(veIndices), 25);
                
                %Plot distribution
%                 if m==2
%                     plot(n.*ones(size(data.x0s(veIndices)))+rand(size(data.x0s(veIndices))), data.x0s(veIndices), '.');
%                 end
                
                Q1P(n,m,whichHemi)=prctile(data.y0s(veIndices), 25);
                Q2P(n,m,whichHemi)=mean(data.y0s(veIndices));
                Q3P(n,m,whichHemi)=prctile(data.y0s(veIndices), 75);
                IQRP(n,m,whichHemi)=prctile(data.y0s(veIndices), 75)-prctile(data.y0s(veIndices), 25);
                
%                 if sum(~isnan(stats{n}.tuningYmajor((m-1)*2+whichHemi,1:10)))>3 && sum(~isnan(stats{n}.tuningYmajor((m-1)*2+whichHemi,12:end)))>3
%                     linearSlopeMajor(n,m,whichHemi)=stats{n}.tuningSlopesMajor((m-1)*2+whichHemi,1);
%                     vSlopeMajor(n,m,whichHemi)=stats{n}.tuningSlopesMajor((m-1)*2+whichHemi,2);
%                     linearSlopeMinor(n,m,whichHemi)=stats{n}.tuningSlopesMinor((m-1)*2+whichHemi,1);
%                     vSlopeMinor(n,m,whichHemi)=stats{n}.tuningSlopesMinor((m-1)*2+whichHemi,2);
%                 else
%                     stats{n}.pSMajProg((m-1)*2+(3-whichHemi) , :)=nan;
%                     stats{n}.pSMinProg((m-1)*2+(3-whichHemi) , :)=nan;
%                 end%Making tables of stats
                
            end
        end
    end
end

%Setting up ANOVA structures
hemisphereGroups=cat(3, ones(size(VEs(:,:,1))), ones(size(VEs(:,:,1)))*2);
tmp=1:length(subjectOrder);
subjectLabels=cat(3, repmat(tmp(:), [1,length(mapNames)]), repmat(tmp(:), [1,length(mapNames)]));
mapLabels=cat(3, repmat(1:length(mapNames), [length(subjectOrder),1]), repmat(1:length(mapNames), [length(subjectOrder),1]));
mapLabels=mapLabels(:);
subjectLabels=subjectLabels(:);
hemisphereLabels=hemisphereGroups(:);
subjectLabels=subjectOrder(subjectLabels);
mapLabels=mapNames(mapLabels);
hemiTmp={'L', 'R'};
hemisphereLabels=hemiTmp(hemisphereLabels);

%Stats and plot of exponents (general code for 3 factors, use when
%hemisphere difference)
[p, tmp, statsOut] = anovan(Exps(:),{hemisphereLabels subjectLabels mapLabels}, 'varnames', {'hemisphere', 'subject', 'map'});
%Multiple comparison tests on ANOVA output
figure; results=multcompare(statsOut, 'Dimension', [1 3])
axis([0 10 0.5 18.5])
axis square

[p, tmp, statsOut] = anovan(SigmaTheta(:),{subjectLabels mapLabels}, 'varnames', {'subject', 'map'});
%Multiple comparison tests on ANOVA output
figure; results=multcompare(statsOut, 'Dimension', [2])
axis([0 pi/2 0.5 18.5])
axis square

%Plot of VEs
[p, tmp, statsOut] = anovan(VEs(:),{hemisphereLabels subjectLabels mapLabels}, 'varnames', {'hemisphere', 'subject', 'map'});
%Multiple comparison tests on ANOVA output
figure; results=multcompare(statsOut, 'Dimension', [1 3])
axis([0 0.6 0.5 20.5])
axis square

%Plot of ROI voxel count
[p, tmp, statsOut] = anovan(nVoxels(:),{hemisphereLabels subjectLabels mapLabels}, 'varnames', {'hemisphere', 'subject', 'map'});
%Multiple comparison tests on ANOVA output
figure; results=multcompare(statsOut, 'Dimension', [1 3])
axis([0 800 0.5 20.5])
axis square

%Plot of sigma majors
[p, tmp, statsOut] = anovan(SigmaMajor(:),{hemisphereLabels subjectLabels mapLabels}, 'varnames', {'hemisphere', 'subject', 'map'});
%Multiple comparison tests on ANOVA output
figure; results=multcompare(statsOut, 'Dimension', [1 3])
axis([0 1.2 0.5 20.5])
axis square
%Stats, because there is no hemisphere effect
[p, tmp, statsOut] = anovan(SigmaMajor(:),{subjectLabels mapLabels}, 'varnames', {'subject', 'map'});

%Plot of sigma minors
[p, tmp, statsOut] = anovan(SigmaMinor(:),{hemisphereLabels subjectLabels mapLabels}, 'varnames', {'hemisphere', 'subject', 'map'});
%Multiple comparison tests on ANOVA output
figure; results=multcompare(statsOut, 'Dimension', [1 3])
axis([0 1.2 0.5 20.5])
axis square

%Plot of sigma ratios
[p, tmp, statsOut] = anovan(SigmaRatio(:),{hemisphereLabels subjectLabels mapLabels}, 'varnames', {'hemisphere', 'subject', 'map'});
%Multiple comparison tests on ANOVA output
figure; results=multcompare(statsOut, 'Dimension', [1 3])
axis([1 7 0.5 20.5])
axis square
%Stats, because there is no hemisphere effect
[p, tmp, statsOut] = anovan(SigmaRatio(:),{subjectLabels mapLabels}, 'varnames', {'subject', 'map'});

%Plot of sigma thetas
[p, tmp, statsOut] = anovan(0.5*pi-SigmaTheta(:),{hemisphereLabels subjectLabels mapLabels}, 'varnames', {'hemisphere', 'subject', 'map'});
%Multiple comparison tests on ANOVA output
figure; results=multcompare(statsOut, 'Dimension', [1 3])
axis([1/8*pi 3/8*pi 0.5 20.5])
axis square
set(gca,'XTick',[pi/8 2*pi/8 3*pi/8])
[p, tmp, statsOut] = anovan(SigmaTheta(:),{subjectLabels mapLabels}, 'varnames', {'subject', 'map'});

%Plot of duration mean
[p, tmp, statsOut] = anovan(Q2D(:),{hemisphereLabels subjectLabels mapLabels}, 'varnames', {'hemisphere', 'subject', 'map'});
%Multiple comparison tests on ANOVA output
figure; results=multcompare(statsOut, 'Dimension', [1 3])
axis([0 1 0.5 20.5])
axis square

%Plot of duration inter quartile range
[p, tmp, statsOut] = anovan(IQRD(:),{hemisphereLabels subjectLabels mapLabels}, 'varnames', {'hemisphere', 'subject', 'map'});
%Multiple comparison tests on ANOVA output
figure; results=multcompare(statsOut, 'Dimension', [1 3])
axis([0 0.5 0.5 20.5])
axis square
[p, tmp, statsOut] = anovan(IQRD(:),{subjectLabels mapLabels}, 'varnames', {'subject', 'map'});

%Plot of period inter quartile range
[p, tmp, statsOut] = anovan(IQRP(:),{hemisphereLabels subjectLabels mapLabels}, 'varnames', {'hemisphere', 'subject', 'map'});
%Multiple comparison tests on ANOVA output
figure; results=multcompare(statsOut, 'Dimension', [1 3])
axis([0 0.5 0.5 20.5])
axis square
[p, tmp, statsOut] = anovan(IQRP(:),{subjectLabels mapLabels}, 'varnames', {'subject', 'map'});

%Cartoon plot of exponent scaling
Yvals=0:0.01:20;
exps=0.1:0.1:1;
figure; hold on;
for count=1:length(exps)
    plot(Yvals, ((Yvals.^exps(count))));
end
axis square

%Getting segmentations and maps as nifti
%Segmentation ROI first. Start by copying an existing folder, installing a segmentation and keeping all
%gray nodes. Create a gray ROI and restrict to layer 1.
mrVista inplane
installSegmentation(1, 1);
mrVista 3;
VOLUME{1}=makeGrayROI(VOLUME{1}); VOLUME{1}=refreshScreen(VOLUME{1},0);
VOLUME{1} = roiRestricttoLayer1(VOLUME{1},VOLUME{1}.selectedROI); 
VOLUME{1} = refreshScreen(VOLUME{1},0);

segName='T1_1mm_class_20.nii';
T1=readFileNifti('t1_1mm.nii');
seg=readFileNifti(segName);
segTmp=seg.data;
seg=T1;
seg.fname=segName;
seg.data=segTmp;
seg
figure; imagesc(seg.data(:,:,100))
writeFileNifti(seg);

segTmp(segTmp==6)=16;
segTmp(segTmp==5)=15;
for n=1:size(VOLUME{1}.ROIs(1).coords, 2)
    segTmp(VOLUME{1}.ROIs(1).coords(3,n), seg.dim(2)+1-VOLUME{1}.ROIs(1).coords(2,n), seg.dim(1)+1-VOLUME{1}.ROIs(1).coords(1,n))=segTmp(VOLUME{1}.ROIs(1).coords(3,n), seg.dim(2)+1-VOLUME{1}.ROIs(1).coords(2,n), seg.dim(1)+1-VOLUME{1}.ROIs(1).coords(1,n))-10;
end
%figure; imagesc(segTmp(:,:,100))
segTmp(segTmp>14)=1;
%figure; imagesc(segTmp(:,:,100))
seg=T1;
seg.data=segTmp;
seg
figure; imagesc(seg.data(:,:,100))
seg.fname='SegmentationLayer1.nii.gz';
writeFileNifti(seg);
%While you have this window open. prepare a box for putting model data in, and save it.
allTmp=0-ones(size(seg.data));
for n=1:size(VOLUME{1}.ROIs(1).coords, 2)
    allTmp(VOLUME{1}.ROIs(1).coords(3,n), seg.dim(2)+1-VOLUME{1}.ROIs(1).coords(2,n), seg.dim(1)+1-VOLUME{1}.ROIs(1).coords(1,n))=0;
end
figure; imagesc(allTmp(:,:,100))
save('allTmp.mat', 'allTmp')

%Open a new mrVista window and load the model. 
%Load the ROI in which the model was run.
%Now, to pull the data out
load('allTmp.mat')
T1=readFileNifti('t1_1mm.nii');
ves=rmGet(VOLUME{1}.rm.retinotopyModels{1}, 've');
x=rmGet(VOLUME{1}.rm.retinotopyModels{1}, 'x');
y=rmGet(VOLUME{1}.rm.retinotopyModels{1}, 'y');

smoothing=5;
if smoothing>0
    [tmp, ii]=intersectCols(VOLUME{1}.coords, VOLUME{1}.ROIs(end).coords);
    mask=zeros(size(VOLUME{1}.coords));
    mask(ii)=1;
    [ves conMat]=dhkGraySmooth(VOLUME{1},ves,[smoothing 0.6], [], mask);
    [x conMat]=dhkGraySmooth(VOLUME{1},x,[smoothing 0.6], conMat, mask);
    [y conMat]=dhkGraySmooth(VOLUME{1},y,[smoothing 0.6], conMat, mask);
end
%Checking it worked, put the result back in the VOLUME{1} window
% VOLUME{1}.rm.retinotopyModels{1}=rmSet(VOLUME{1}.rm.retinotopyModels{1}, 'x', x);
% VOLUME{1}.rm.retinotopyModels{1}=rmSet(VOLUME{1}.rm.retinotopyModels{1}, 'y', y);
% VOLUME{1}.rm.retinotopyModels{1}=rmSet(VOLUME{1}.rm.retinotopyModels{1}, 've', ves);
% VOLUME{1} = rmLoad(VOLUME{1}, 1, 'x0', 'map');
% VOLUME{1} = rmLoad(VOLUME{1}, 1, 'y0', 'amp');
% VOLUME{1} = rmLoad(VOLUME{1}, 1, 'varexplained', 'co');
% VOLUME{1} = refreshScreen(VOLUME{1}, 1);

ves(x>1)=0;
ves(y>1)=0;
ves(x<0.05)=0;
ves(y<0.05)=0;
xtmp=x;
ytmp=y;
xtmp(x>1)=0;
xtmp(y>1)=0;
xtmp(x<0.05)=0;
xtmp(y<0.05)=0;
ytmp(x>1)=0;
ytmp(y>1)=0;
ytmp(x<0.05)=0;
ytmp(y<0.05)=0;
x=xtmp;
y=ytmp;
%Determine which indices to set
[tmp, ii]=intersectCols(VOLUME{1}.coords, VOLUME{1}.ROIs(1).coords);
veTmp=allTmp;
xTmp=allTmp;
yTmp=allTmp;
for n=1:length(ii)
    veTmp(VOLUME{1}.coords(3,ii(n)), T1.dim(2)+1-VOLUME{1}.coords(2,ii(n)), T1.dim(1)+1-VOLUME{1}.coords(1,ii(n)))=ves(ii(n));
    xTmp(VOLUME{1}.coords(3,ii(n)), T1.dim(2)+1-VOLUME{1}.coords(2,ii(n)), T1.dim(1)+1-VOLUME{1}.coords(1,ii(n)))=x(ii(n));
    yTmp(VOLUME{1}.coords(3,ii(n)), T1.dim(2)+1-VOLUME{1}.coords(2,ii(n)), T1.dim(1)+1-VOLUME{1}.coords(1,ii(n)))=y(ii(n));
end
venii=T1;
xii=T1;
yii=T1;
venii.data=veTmp;
xii.data=xTmp;
yii.data=yTmp;
if smoothing>0
    venii.fname='VarExpSmooth.nii.gz';
    xii.fname='PrefDurationSmooth.nii.gz';
    yii.fname='PrefPeriodSmooth.nii.gz';   
else
    venii.fname='VarExp.nii.gz';
    xii.fname='PrefDuration.nii.gz';
    yii.fname='PrefPeriod.nii.gz';
end
writeFileNifti(venii);
writeFileNifti(xii);
writeFileNifti(yii);

%AFNI clustering process
paths{1}='/media/Storage3/WholeBrainNewSeg/ForAFNIclustering/FromAlessioMatlab/S08_Data';
paths{2}='/media/Storage3/WholeBrainNewSeg/ForAFNIclustering/FromAlessioMatlab/S09_Data';
paths{3}='/media/Storage3/WholeBrainNewSeg/ForAFNIclustering/FromAlessioMatlab/S10_Data';
paths{4}='/media/Storage3/WholeBrainNewSeg/ForAFNIclustering/FromAlessioMatlab/S11ForAFNI';
paths{5}='/media/Storage3/WholeBrainNewSeg/ForAFNIclustering/FromAlessioMatlab/S12_Data';
paths{6}='/media/Storage3/WholeBrainNewSeg/ForAFNIclustering/FromAlessioMatlab/S13_Data';
paths{7}='/media/Storage3/WholeBrainNewSeg/ForAFNIclustering/FromAlessioMatlab/S01ForAFNI';
paths{8}='/media/Storage3/WholeBrainNewSeg/ForAFNIclustering/FromAlessioMatlab/S03ForAFNI';

%selSizes=26:2:102; %Original 
selSizes=86;%;50:2:102;
%selSizes=76;  
selVarExps= 0.76:0.002:0.998;  %For percentiles
selVarExps=0.01:0.01:0.40; %For actual VE
actualVE=1;
% #selSizes <- c(65, 75)
% #selVarExps <- c( 0.895, 0.897, 0.899 )

outputClusteringSubjLeft=zeros([length(selSizes),length(selVarExps)]);
outputClusteringSubjRight=zeros([length(selSizes),length(selVarExps)]);
mainDir='/media/Storage3/WholeBrainNewSeg/ForAFNIclustering/FromAlessioMatlab';
%mainDir <- 'C:/Users/alessiof/Dropbox/miscellaneousCode/BenCurrBio_revision'

subjDirectories=paths;
whichSubs=[5:8];
 parpool(4)
parfor nSubj=1:length(whichSubs);%(subjDirectories)
  
  subjDir=subjDirectories{whichSubs(nSubj)}
  cd(subjDir)
  
  !rm *ttt.nii.gz;
  !rm -R surfaces_folder_left;
  !rm -R surfaces_folder_right;
  !rm -R modelResults_ttt_interp_surfaces_folder_left;
  !rm -R modelResults_ttt_interp_surfaces_folder_right;
  
  instr= ['!3dcalc -a VarExp.nii.gz -expr "a*step(a)" -prefix VarExp_corr_ttt.nii.gz'];
  commitEval(instr)    
  %combine modelling results
  !3dTcat -prefix modelResults_ttt.nii.gz PrefDuration.nii.gz PrefPeriod.nii.gz VarExp_corr_ttt.nii.gz; 
  %modeling in orig space
  !3drefit -space ORIG -view orig modelResults_ttt.nii.gz; 
  %anat in orig space
  !3dcopy t1_1mm.nii.gz t1_1mm_copy_ttt.nii.gz; 
  !3drefit -space ORIG -view orig t1_1mm_copy_ttt.nii.gz;
  %segmentation in orig space
  !3dcopy SegmentationLayer1.nii.gz SegmentationLayer1_ttt.nii.gz;
  !3drefit -space ORIG -view orig SegmentationLayer1_ttt.nii.gz;
  %Left segmentation
  instr=['!3dcalc -a SegmentationLayer1_ttt.nii.gz -expr "within(a,2.9,3.1)" -prefix leftSeg_ttt.nii.gz'];
  commitEval(instr)
  %Right segmentation
  instr=['!3dcalc -a SegmentationLayer1_ttt.nii.gz -expr "within(a,3.9,4.1)" -prefix rightSeg_ttt.nii.gz'];
  commitEval(instr)
  
  %generate surfaces
  instr=sprintf('!Rscript %s/generateSurfacesFromBoundaries.R leftSeg_ttt.nii.gz 10 0 800 1', mainDir);
  commitEval(instr)
  !mv surfaces_folder/ surfaces_folder_left;
  instr=sprintf('!Rscript %s/generateSurfacesFromBoundaries.R rightSeg_ttt.nii.gz 10 0 800 1', mainDir);
  commitEval(instr)
  !mv surfaces_folder/ surfaces_folder_right/;
  
  %copy maps on surface_left
  instr=sprintf('!Rscript %s/interpolateSurface.R modelResults_ttt.nii.gz surfaces_folder_left/', mainDir); 
  commitEval(instr)
  %copy maps on surface_right
  instr=sprintf('!Rscript %s/interpolateSurface.R modelResults_ttt.nii.gz surfaces_folder_right/', mainDir);
  commitEval(instr)
  
%   %load the data, to test
%   sizeThr=25;
%   varThr=0.1;
%   %cluster the surfaces based on variance explained
%   instr=sprintf('!SurfClust -spec surfaces_folder_left/spec.surfaces.smoothed -surf_A surfaces_folder_left/boundary00_sm.1D.coord  -input modelResults_ttt_interp_surfaces_folder_left/boundary00_sm_modelResults_ttt_surf.1D.dset 8 -rmm -1 -amm2 %s -thresh_col 8 -thresh %s > leftTable_ttt.1D.txt', num2str(sizeThr), num2str(varThr) );
%   commitEval(instr)
%   instr=sprintf('!SurfClust -spec surfaces_folder_right/spec.surfaces.smoothed -surf_A surfaces_folder_right/boundary00_sm.1D.coord  -input modelResults_ttt_interp_surfaces_folder_right/boundary00_sm_modelResults_ttt_surf.1D.dset 8 -rmm -1 -amm2 %s -thresh_col 8 -thresh %s > rightTable_ttt.1D.txt', num2str(sizeThr), num2str(varThr) );
%   commitEval(instr)
%   %cluster the surfaces based on variance explained, show table
% %   instr=sprintf('!SurfClust -spec surfaces_folder_left/spec.surfaces.smoothed -surf_A surfaces_folder_left/boundary00_sm.1D.coord  -input modelResults_ttt_interp_surfaces_folder_left/boundary00_sm_modelResults_ttt_surf.1D.dset 8 -rmm -1 -amm2 %s -thresh_col 8 -thresh %s', num2str(sizeThr), num2str(varThr) ); 
% %   commitEval(instr)
% %   instr=sprintf('!SurfClust -spec surfaces_folder_right/spec.surfaces.smoothed -surf_A surfaces_folder_right/boundary00_sm.1D.coord  -input modelResults_ttt_interp_surfaces_folder_right/boundary00_sm_modelResults_ttt_surf.1D.dset 8 -rmm -1 -amm2 %s -thresh_col 8 -thresh %s', num2str(sizeThr), num2str(varThr) );
% %   commitEval(instr)
% %   leftTable <- read.table('leftTable_ttt.1D', as.is=TRUE)
% %   rightTable <- read.table('rightTable_ttt.1D', as.is=TRUE)
%   !cp leftTable_ttt.1D leftTable_ttt.1D.txt
%   tmp=readtable('leftTable_ttt.1D.txt');%, 'Delimiter', '\t', 'HeaderLines', 25, 'ReadVariableNames', true);
%   tmp=tmp(26:end, 1);
%   leftTable=(table2array(tmp));
%   !cp rightTable_ttt.1D rightTable_ttt.1D.txt
%   tmp=readtable('rightTable_ttt.1D.txt');
%   tmp=tmp(26:end, 1);
%   rightTable=(table2array(tmp));  %Split removed here
  
  
%   # # load the data, threshold and clusterize according to pars
%   # sizeThrLoop <- seq( 2, 80, 4 )
%   # varThrLoop <- seq( 0.06, 0.9, 0.04 )
%   # outputClusteringLeft <- array(999,c(length(sizeThrLoop), length(varThrLoop)))
%   # outputClusteringRight <- array(999,c(length(sizeThrLoop), length(varThrLoop)))
%   # for ( sizeCount in 1:length(sizeThrLoop) ) {
%   #   for ( varCount in 1:length(varThrLoop) ) {
%   #     selSize <- sizeThrLoop[sizeCount]
%   #     selVar <- varThrLoop[varCount]
%   #     instr <- sprintf( 'SurfClust -spec surfaces_folder_left/spec.surfaces.smoothed -surf_A surfaces_folder_left/boundary00_sm.1D.coord  -input modelResults_ttt_interp_surfaces_folder_left/boundary00_sm_modelResults_ttt_surf.1D.dset 8 -rmm -1 -amm2 %s -thresh_col 8 -thresh %s > leftTableLoop_ttt.1D', as.character(selSize), as.character(selVar) ); system( instr ); # cluster the surface left based on variance explained
%   #     instr <- sprintf( 'SurfClust -spec surfaces_folder_right/spec.surfaces.smoothed -surf_A surfaces_folder_right/boundary00_sm.1D.coord  -input modelResults_ttt_interp_surfaces_folder_right/boundary00_sm_modelResults_ttt_surf.1D.dset 8 -rmm -1 -amm2 %s -thresh_col 8 -thresh %s > rightTableLoop_ttt.1D', as.character(selSize), as.character(selVar) ); system( instr ); # cluster the surface right based on variance explained
%   #     leftTableLoop <- read.table('leftTableLoop_ttt.1D', as.is=TRUE)
%   #     rightTableLoop <- read.table('rightTableLoop_ttt.1D', as.is=TRUE)
%   #     
%   #     if ( sum( dim(leftTableLoop)==c(1,3) )==2 ) {
%   #       outputClusteringLeft[sizeCount, varCount] <- 0
%   #     } else {
%   #       outputClusteringLeft[sizeCount, varCount] <- dim( leftTableLoop )[1]
%   #     }
%   #     
%   #     if ( sum( dim(rightTableLoop)==c(1,3) )==2 ) {
%   #       outputClusteringRight[sizeCount, varCount] <- 0
%   #     } else {
%   #       outputClusteringRight[sizeCount, varCount] <- dim( rightTableLoop )[1]
%   #     }
%   #     
%   #   }
%   # }
  
  % load model data for each hemisphere and get varExp distribution (quantiles)
  !cp modelResults_ttt_interp_surfaces_folder_left/boundary00_sm_modelResults_ttt_surf.1D.dset modelResults_ttt_interp_surfaces_folder_left/boundary00_sm_modelResults_ttt_surf.1D.dset.txt
  modelLeft=readtable('modelResults_ttt_interp_surfaces_folder_left/boundary00_sm_modelResults_ttt_surf.1D.dset.txt');
  modelLeft=table2array(modelLeft);
  !cp modelResults_ttt_interp_surfaces_folder_right/boundary00_sm_modelResults_ttt_surf.1D.dset modelResults_ttt_interp_surfaces_folder_right/boundary00_sm_modelResults_ttt_surf.1D.dset.txt
  modelRight=readtable('modelResults_ttt_interp_surfaces_folder_right/boundary00_sm_modelResults_ttt_surf.1D.dset.txt');
  modelRight=table2array(modelRight);
  
  varExp_nodeLeft=modelLeft((find(modelLeft(:,9)>0.001)), 9);
  varExp_nodeRight=modelRight((find(modelRight(:,9)>0.001)) ,9);
%   #sizeThrLoop <- seq( 15, 85, 15 )
%   #varQuantThrLoopLeft <- quantile( varExp_nodeLeft, c( seq( 0.40, 0.95, 0.05 ), seq( 0.96, 0.998, 0.002 ) ) )
%   #varQuantThrLoopRight <- quantile( varExp_nodeRight,  c( seq( 0.40, 0.95, 0.05 ), seq( 0.96, 0.998, 0.002 ) )  )
  sizeThrLoop=selSizes;
  
  %if using percentile
  if actualVE==0
      varQuantThrLoopLeft=quantile(varExp_nodeLeft, selVarExps);
      varQuantThrLoopRight=quantile(varExp_nodeRight,  selVarExps);
  else
      varQuantThrLoopLeft=selVarExps;
      varQuantThrLoopRight=selVarExps;
  end
  
  outputClusteringLeftQuant=zeros(length(sizeThrLoop), length(varQuantThrLoopLeft))+999;
  outputClusteringRightQuant=zeros(length(sizeThrLoop), length(varQuantThrLoopRight))+999;
  for sizeCount=1:length(sizeThrLoop)
    for varCount=1:length(varQuantThrLoopLeft)
      
      thisSelSize=sizeThrLoop(sizeCount);
      selVarLeft=varQuantThrLoopLeft(varCount);
      selVarRight=varQuantThrLoopRight(varCount);
      
      disp( sprintf('size: %s, varExp (left): %s, varExp (right): %s', num2str(thisSelSize ), num2str(round(selVarLeft, 3 )), num2str( round( selVarRight, 3 ) ) ) )
      % cluster the surfaces based on variance explained
      instr=sprintf( '!SurfClust -spec surfaces_folder_left/spec.surfaces.smoothed -surf_A surfaces_folder_left/boundary00_sm.1D.coord  -input modelResults_ttt_interp_surfaces_folder_left/boundary00_sm_modelResults_ttt_surf.1D.dset 8 -rmm -1 -amm2 %s -thresh_col 8 -thresh %s > leftTableLoop_ttt.1D', num2str(thisSelSize), num2str(selVarLeft) );
      commitEval(instr);
      instr=sprintf( '!SurfClust -spec surfaces_folder_right/spec.surfaces.smoothed -surf_A surfaces_folder_right/boundary00_sm.1D.coord  -input modelResults_ttt_interp_surfaces_folder_right/boundary00_sm_modelResults_ttt_surf.1D.dset 8 -rmm -1 -amm2 %s -thresh_col 8 -thresh %s > rightTableLoop_ttt.1D', num2str(thisSelSize), num2str(selVarRight) );
      commitEval(instr);
      
      !cp leftTableLoop_ttt.1D leftTableLoop_ttt.1D.txt
      !cp rightTableLoop_ttt.1D rightTableLoop_ttt.1D.txt
      
      fid=fopen('leftTableLoop_ttt.1D.txt');
      leftTableLoop=textscan(fid, repmat('%f', 1, 23), 'HeaderLines', 26);
      fid=fopen('rightTableLoop_ttt.1D.txt');
      rightTableLoop=textscan(fid, repmat('%f', 1, 23), 'HeaderLines', 26);
      %leftTableLoop=leftTableLoop(26:end, 1);
      %leftTableLoop=split(table2array(leftTableLoop));
%       rightTableLoop=readtable('rightTableLoop_ttt.1D.txt');
%       rightTableLoop=rightTableLoop(26:end, 1);
      %rightTableLoop=split(table2array(rightTableLoop));
      
        outputClusteringLeftQuant(sizeCount, varCount)=size(leftTableLoop{1},1);
      
        outputClusteringRightQuant(sizeCount, varCount)=size(rightTableLoop{1},1);
    end
  end
  cd(mainDir)
  filename=sprintf('dataClustering_Subject_%s_Date_%s.mat', num2str(whichSubs(nSubj)), datestr(now));
  commitSave(filename, outputClusteringLeftQuant, outputClusteringRightQuant);
end


    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clustering trends and plot #
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(mainDir);
allSubjLeftQuant=zeros([length(selSizes), length(selVarExps), length(paths)]);
allSubjRightQuant=zeros([length(selSizes), length(selVarExps), length(paths)]);
for nSubj=5;%1:length(paths)
    fname=dir(sprintf('dataClustering_Subject_%s_Date*', num2str(nSubj)));
    load(fname.name)
    allSubjLeftQuant(:,:,nSubj)=outputClusteringLeftQuant;
    allSubjRightQuant(:,:,nSubj)=outputClusteringRightQuant;
end

% takes the average clustering result between hemispheres
clustOutput=(allSubjLeftQuant + allSubjRightQuant) / 2;
% takes the average clustering result between participants (or optionally
% choose one participant)
clustOutputMean =clustOutput(:,:,5);  %mean(clustOutput, 3);
% takes the standard error of the mean clustering result
clustOutputSe =std(clustOutput, 0, 3) ./ sqrt( length( paths ) );

% #clustOutput <- abind( outputClusteringSubjLeft, outputClusteringSubjRight, along=3 ) 
% #clustOutputMean <- apply( clustOutput, c(1,2), mean )
% #clustOutputSe <- apply( clustOutput, c(1,2), sd ) / sqrt( dim( clustOutput )[3] )

% performs the linear and piecewise regression 
clusterSize=86;
% set.seed(123)
% graphics.off()
% x11(width = 5, height = 4.5)
% library(segmented)
%this is the selected clustering size corresponding to selSizes[31] = 86 squared millimeters
lineCode =find(selSizes==clusterSize); 
% y axis (cluster number) limit
maxY =16; 
%X axis candidates to plot
plotVals=1:(size(clustOutputMean,2));
figure; plot(selVarExps(plotVals), clustOutputMean(lineCode,plotVals), 'k', 'LineWidth', 2);
axis([selVarExps(plotVals(1)) selVarExps(plotVals(end)) 0 maxY]);
axis square;
ylabel('Number of clusters');
xlabel('Variance explained');
% axis(1, selVarExps[plotVals[seq(1,length(plotVals),15)]], round(selVarExps[plotVals[seq(1,length(plotVals),15)]],3) )
% axis(2, seq(0,maxY,10), seq(0,maxY,10), las=1  )
hold on; plot(selVarExps(plotVals), clustOutputMean(lineCode,plotVals)-clustOutputSe(lineCode,plotVals), 'Color', [0.5 0.5 0.5]);
hold on; plot(selVarExps(plotVals), clustOutputMean(lineCode,plotVals)+clustOutputSe(lineCode,plotVals), 'Color', [0.5 0.5 0.5]);
%Linear fit
afit=linreg(selVarExps(plotVals), clustOutputMean(lineCode,plotVals));
xfit=[selVarExps(plotVals(1)) selVarExps(plotVals(end))];
hold on; plot(xfit, xfit.*afit(2)+afit(1), 'b');

%Piecewise broken line fit
xs=selVarExps(plotVals);
ys=clustOutputMean(lineCode,plotVals);
ptdif=diff(ys')./diff(xs');
L1=[true; ptdif<2];
L2=~L1;
B1=polyfit(xs(L1), ys(L1), 1);
F1=polyval(B1, xs(L1));
B2=polyfit(xs(L2), ys(L2), 1);
F2=polyval(B2, xs(L2));
hold on; plot(xs(L1), F1);

dati <- data.frame( y = clustOutputMean[lineCode,plotVals], x = selVarExps[plotVals] ) 
out.lm <- lm(y ~ x, data = dati)
o <- segmented( out.lm, seg.Z = ~x, psi = list(x = c(0.85,0.92)),
                control = seg.control(display = FALSE), n.boot=0, it.max = 1000, random=TRUE, quant=FALSE )
dat2 = data.frame(x = selVarExps[plotVals], y = broken.line(o)$fit)
abline( out.lm, col='blue', lwd=2  )
lines( dat2, col='darkorange', lwd=2 )
o
selSizes[lineCode]
storedClustVals <- clustOutputMean[lineCode,plotVals]
storedVars <- selVarExps[plotVals]
mean( storedClustVals[ o$id.group==1 ] )
sd( storedClustVals[ o$id.group==1 ] )
startPoint <- which( o$id.group==1 )[1]
endPoint <- which( o$id.group==2 )[1] - 1
abline( v=selVarExps[ plotVals[startPoint] ], lty=2, lwd=1, col='gray50' )
abline( v=selVarExps[ plotVals[endPoint] ], lty=2, lwd=1, col='gray50' )
abline( h=mean( storedClustVals[ o$id.group==1 ] ), lty=2, lwd=1, col='gray50'   )
legend( 0.88, 20, legend=c('linear fit','piecewise fit'), fill=c('blue','darkorange'), bty='n' )
filename <- sprintf('clust_%1.0f_meanFlat_%1.2f_sdFlat_%1.2f.pdf', selSizes[lineCode], round( mean( storedClustVals[ o$id.group==1 ] ),2), round(sd( storedClustVals[ o$id.group==1 ] ),2) )
lmValsClust <- storedClustVals[ o$id.group==1 ]; lmValsClust <- lmValsClust[ seq(1,length(lmValsClust),length.out = 5) ]
lmValsVar <- storedVars[ o$id.group==1 ]; lmValsVar <- storedVars[ seq(1,length(lmValsVar),length.out = 5) ]
summary( lm( lmValsClust ~ lmValsVar ) )
lmValsClust <- storedClustVals[ o$id.group==0 ]; lmValsClust <- lmValsClust[ seq(1,length(lmValsClust),length.out = 5) ]
lmValsVar <- storedVars[ o$id.group==0 ]; lmValsVar <- storedVars[ seq(1,length(lmValsVar),length.out = 5) ]
summary( lm( lmValsClust ~ lmValsVar ) )
lmValsClust <- storedClustVals[ o$id.group==2 ]; lmValsClust <- lmValsClust[ seq(1,length(lmValsClust),length.out = 5) ]
lmValsVar <- storedVars[ o$id.group==2 ]; lmValsVar <- storedVars[ seq(1,length(lmValsVar),length.out = 5) ]
summary( lm( lmValsClust ~ lmValsVar ) )
anova( out.lm, o )
dev.copy2pdf( file = filename )
dev.off()

outputDir <- 'roisClustered';
source( sprintf('%s/AFNIio.R',mainDir) )
source( sprintf('%s/coordinateFromLinearIndex.r',mainDir) )
source( sprintf('%s/linearIndexFromCoordinate.r',mainDir) )
library('RANN')

%based on selSizes[31] = 86 squared millimeters
selectedSize=86;
%based on the mean( storedVars[ o$id.group==1 ] ), average of the variance explained quantile for the 'flat' portion of the piecewise regression
selectedQuantile=0.90;
actualVE=1;
selectedVarExp=0.12
for nSubj=5%1:length(subjDirectories)
  subjDir=subjDirectories{nSubj}
  cd( subjDir )
  !rm *0.1D *.niml.dset roisVol_left_zzz.nii.gz roisVol_right_zzz.nii.gz *_zzz.nii.gz *_zzz.nii.gz
  %load model data for each hemisphere and get varExp distribution (quantiles)
  modelLeft=readtable('modelResults_ttt_interp_surfaces_folder_left/boundary00_sm_modelResults_ttt_surf.1D.dset.txt');
  modelLeft=table2array(modelLeft);
  modelRight=readtable('modelResults_ttt_interp_surfaces_folder_right/boundary00_sm_modelResults_ttt_surf.1D.dset.txt');
  modelRight=table2array(modelRight);
  varExp_nodeLeft=modelLeft((find(modelLeft(:,9)>0.001)), 9);
  varExp_nodeRight=modelRight((find(modelRight(:,9)>0.001)) ,9);
  sizeThrLoop=selectedSize;
  varQuantThrLoopLeft=quantile(varExp_nodeLeft, selectedQuantile);
  varQuantThrLoopRight=quantile(varExp_nodeRight,  selectedQuantile);
  if exist('actualVE', 'var') && actualVE==1
    varQuantThrLoopLeft=selectedVarExp;
    varQuantThrLoopRight=selectedVarExp;
  end
  %cluster the surface left based on variance explained
  disp( sprintf('size: %s, varExp (left): %s, varExp (right): %s', num2str(selectedSize ), num2str(round(varQuantThrLoopLeft, -3 )), num2str( round( varQuantThrLoopRight, -3 ) ) ) )
  instr=sprintf( '!SurfClust -spec surfaces_folder_left/spec.surfaces.smoothed -surf_A surfaces_folder_left/boundary00_sm.1D.coord  -input modelResults_ttt_interp_surfaces_folder_left/boundary00_sm_modelResults_ttt_surf.1D.dset 8 -rmm -1 -amm2 %s -thresh_col 8 -thresh %s -out_roidset -prefix left_clus_roi_', num2str(sizeThrLoop), num2str(varQuantThrLoopLeft) );   
  commitEval(instr);
  instr=sprintf( '!SurfClust -spec surfaces_folder_right/spec.surfaces.smoothed -surf_A surfaces_folder_right/boundary00_sm.1D.coord  -input modelResults_ttt_interp_surfaces_folder_right/boundary00_sm_modelResults_ttt_surf.1D.dset 8 -rmm -1 -amm2 %s -thresh_col 8 -thresh %s -out_roidset -prefix right_clus_roi_', num2str(sizeThrLoop), num2str(varQuantThrLoopRight) ); 
  commitEval(instr);
  
  !3dSurf2Vol -spec surfaces_folder_left/spec.surfaces.smoothed -surf_A surfaces_folder_left/boundary00_sm.1D.coord -sdata left_clus_roi__ClstMsk_e1_a86.0.niml.dset -grid_parent SegmentationLayer1_ttt.nii.gz -sv SegmentationLayer1_ttt.nii.gz -map_func ave -prefix roisVol_left_zzz.nii.gz;
  !3dSurf2Vol -spec surfaces_folder_right/spec.surfaces.smoothed -surf_A surfaces_folder_right/boundary00_sm.1D.coord -sdata right_clus_roi__ClstMsk_e1_a86.0.niml.dset -grid_parent SegmentationLayer1_ttt.nii.gz -sv SegmentationLayer1_ttt.nii.gz -map_func ave -prefix roisVol_right_zzz.nii.gz;
  roisLeftFile=readFileNifti('roisVol_left_zzz.nii.gz');
  roisRightFile=readFileNifti('roisVol_right_zzz.nii.gz');
  segmentationFile=readFileNifti('SegmentationLayer1_ttt.nii.gz');
  
  roisLeftVolume=roisLeftFile.data;
  roisRightVolume=roisRightFile.data;
  segAppVol=segmentationFile.data;
  indexLeft=find(segAppVol==5);
  indexRight=find(segAppVol==6);
  [x y z]=ind2sub(size(segAppVol), indexLeft);
  coordsLeft=[x,y,z];
  [x y z]=ind2sub(size(segAppVol), indexRight);
  coordsRight=[x,y,z];
  [x y z]=ind2sub(size(roisLeftVolume), find(roisLeftVolume>0));
  coordsRoiLeft=[x,y,z];
  [x y z]=ind2sub(size(roisRightVolume), find(roisRightVolume>0));
  coordsRoiRight=[x,y,z];
  
  distOutLeft=nearpoints(coordsRoiLeft',coordsLeft');
  distOutRight=nearpoints(coordsRoiRight',coordsRight');
  leftVolOut=zeros(size(segAppVol)); 
  rightVolOut=zeros(size(segAppVol));
  
  leftVolOut(indexLeft(distOutLeft))=roisLeftVolume(roisLeftVolume>0);
  rightVolOut(indexRight(distOutRight))=roisRightVolume(roisRightVolume>0);
  
  segmentationFile.fname=sprintf('roiLeftLayer1_Sub%s_zzz.nii.gz',num2str(nSubj))
  segmentationFile.data=leftVolOut;
  writeFileNifti(segmentationFile);
  segmentationFile.fname=sprintf('roiRightLayer1_Sub%s_zzz.nii.gz',num2str(nSubj))
  segmentationFile.data=rightVolOut;
  writeFileNifti(segmentationFile);


  instr=sprintf('!cp roiLeftLayer1_Sub%s_zzz.nii.gz %s/roiLeftLayer1_Sub%s_zzz.nii.gz',num2str(nSubj), mainDir,num2str(nSubj));
  commitEval(instr);
  instr=sprintf('!cp roiRightLayer1_Sub%s_zzz.nii.gz %s/roiRightLayer1_Sub%s_zzz.nii.gz',num2str(nSubj), mainDir, num2str(nSubj));
  commitEval(instr);
  
end

%project maps on standard template 
paths{1}='/media/Storage3/WholeBrainNewSeg/ForAFNIclustering/FromAlessioMatlab/S08_Data';
paths{2}='/media/Storage3/WholeBrainNewSeg/ForAFNIclustering/FromAlessioMatlab/S09_Data';
paths{3}='/media/Storage3/WholeBrainNewSeg/ForAFNIclustering/FromAlessioMatlab/S10_Data';
paths{4}='/media/Storage3/WholeBrainNewSeg/ForAFNIclustering/FromAlessioMatlab/S11ForAFNI';
paths{5}='/media/Storage3/WholeBrainNewSeg/ForAFNIclustering/FromAlessioMatlab/S12_Data';
paths{6}='/media/Storage3/WholeBrainNewSeg/ForAFNIclustering/FromAlessioMatlab/S13_Data';
paths{7}='/media/Storage3/WholeBrainNewSeg/ForAFNIclustering/FromAlessioMatlab/S01ForAFNI';
paths{8}='/media/Storage3/WholeBrainNewSeg/ForAFNIclustering/FromAlessioMatlab/S03ForAFNI';

mainDir='/media/Storage3/WholeBrainNewSeg/ForAFNIclustering/FromAlessioMatlab';
templateDir='/media/Storage3/WholeBrainNewSeg/ForAFNIclustering/FromAlessioMatlab/suma_TT_N27';

subjDirectories=paths;
%prepares files for each individual
for nSubj=1:length(subjDirectories)
  subjDir=subjDirectories{nSubj};
  cd(subjDir); 
  !rm *_nnn.nii.gz 
  !rm Qwarp*; 
  !rm *_nnn_*';
  !3dcopy Segmentation.nii.gz Segmentation_nnn.nii.gz;
  !3drefit -space ORIG -view orig Segmentation_nnn.nii.gz;
  !3dcalc -a Segmentation_nnn.nii.gz -b t1_1mm_copy_ttt.nii.gz -expr "b*step(a-1)" -prefix t1_skullStripped_nnn.nii.gz;
end
%prepares files in the template
cd(templateDir)
!rm *_nnn.nii;
!rm *_nnn.nii.gz;
!3dcopy aparc+aseg.nii aparc+aseg_nnn.nii;
!3dcalc -a aparc+aseg_nnn.nii -b aparc+aseg_nnn.nii -expr " not( or( within(a,6.5,8.5), within(a,45.5,47.5), within(a,14.5,16.5) ) )*b " -prefix aparc+aseg_nocer_nnn.nii;
!3dcalc -a TT_N27_SurfVol.nii -b aparc+aseg_nocer_nnn.nii -expr " a*step(b) " -prefix TT_N27_SurfVol_target_nnn.nii;

%non-lin coregistration to the template for each individual 
for nSubj=1:length(subjDirectories)
    subjDir=subjDirectories{nSubj};
    cd(subjDir);
    !rm t1_skullStripped_nnn_at.nii t1_skullStripped_nnn_at.Xat.1D t1_skullStripped_nnn_at.nii_WarpDrive.log t1_skullStripped_nnn_at.nii.Xaff12.1D;
  
    instr=sprintf('!@auto_tlrc -base %s/TT_N27_SurfVol_target_nnn.nii -input t1_skullStripped_nnn.nii.gz -no_ss -init_xform AUTO_CENTER', templateDir);
    commitEval(instr);
    instr=sprintf('!3dQwarp -iwarp -blur 0 0 -base %s/TT_N27_SurfVol_target_nnn.nii -source t1_skullStripped_nnn_at.nii', templateDir );
    commitEval(instr);
end

%Specify whichROIs
mapNames=["TLO", "TTOP", "TTOA", "TPO", "TLS", "TPCI", "TPCM", "TPCS", "TFI" "TFS"];
hemispheres=["Left", "Right"];
whichRoi=["Map", "Center"];
subOrder={'S8', 'S9', 'S10', 'S11', 'S12', 'S13', 'Harv', 'Nell'};

pathsRoi{1}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/S08_Combined';
pathsRoi{2}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/S09_Combined';
pathsRoi{3}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/S10_Combined';
pathsRoi{4}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/S11_Combined';
pathsRoi{5}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/S12_Combined';
pathsRoi{6}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/S13_Combined_Seg20';
pathsRoi{7}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/HarvCombined';
pathsRoi{8}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/NellCombined';

%projects model results and clustering results on the standard template for each participant
cd(templateDir) 
!rm *_fff.nii.gz;
for nSubj=1:length(subjDirectories)
  subjDir=subjDirectories{nSubj}; 
  cd(subjDir);
  !rm *_fff.nii.gz;  
  !rm Roi*.nii.gz; 
  stringSplit=split(subjDir, '/');
  subjName=char(stringSplit(end));
  
  segmentationAfni=readFileNifti('SegmentationLayer1.nii.gz');
  segAppVol=segmentationAfni.data;
  indexLeft=find(segAppVol==5);
  indexRight=find(segAppVol==6);
  indexAll=find(segAppVol==6 | segAppVol==5);
  [x y z]=ind2sub(size(segAppVol), indexLeft);
  coordsLeft=[x,y,z];
  [x y z]=ind2sub(size(segAppVol), indexRight);
  coordsRight=[x,y,z];
  [x y z]=ind2sub(size(segAppVol), indexAll);
  coordsAll=[x,y,z];
  
  
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%   %Pass gray-Layer1 ROI through warp to determine correspondance. Useful
%   %for testing, not used in main analysis
%   load(char(strcat(pathsRoi{nSubj}, '/Gray/ROIs/gray-Layer1.mat')));
%   t1=readFileNifti('Segmentation.nii.gz');
%   t1.fname=char(strcat('RoiGray-Layer1.nii.gz'));
%   t1.data=zeros(size(t1.data));
% 
%   whichAfniVertex=nearpoints(double([ROI.coords(3,:); (size(t1.data,2)-ROI.coords(2,:)); (size(t1.data,3)-ROI.coords(1,:))]), coordsAll');
%   t1.data(indexAll(whichAfniVertex))=1;
%   whichAfniVertex=nearpoints(coordsAll', double([ROI.coords(3,:); (size(t1.data,2)-ROI.coords(2,:)); (size(t1.data,3)-ROI.coords(1,:))]));
%   t1.data(indexAllVista(whichAfniVertex))=1;    
%   
%   writeFileNifti(t1);
%   %This part warps the result to template
%   instr=sprintf('!3dAllineate -1Dmatrix_apply t1_skullStripped_nnn_at.Xat.1D -master %s/TT_N27_SurfVol_target_nnn.nii -interp linear -prefix %sInAnatomy_fff.nii.gz -input %s', templateDir, char(strcat('RoiGray-Layer1')), char(strcat('RoiGray-Layer1.nii.gz')));
%   commitEval(instr);
%   instr=sprintf('!3dNwarpApply -nwarp Qwarp_WARP+tlrc -source %sInAnatomy_fff.nii.gz -master %s/TT_N27_SurfVol_target_nnn.nii -prefix %sInAnatomy_warp_fff.nii.gz -ainterp linear', char(strcat('RoiGray-Layer1')), templateDir, char(strcat('RoiGray-Layer1')));
%   commitEval(instr);
%   instr=sprintf('!cp %sInAnatomy_warp_fff.nii.gz %s/%s/%sInAnatomy_mask_fff.nii.gz', char(strcat('RoiGray-Layer1')), templateDir, 'Map', char(strcat(subOrder{nSubj}, 'Gray-Layer1')));
%   commitEval(instr);
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





    %Pass All Map ROI through warp to determine correspondance
  load(char(strcat(pathsRoi{nSubj}, '/Gray/ROIs/LeftAllMap.mat')));
  t1=readFileNifti('Segmentation.nii.gz');
  t1.fname=char(strcat('RoiLeftAllMap.nii.gz'));
  t1.data=zeros(size(t1.data));

  whichAfniVertex=nearpoints(double([ROI.coords(3,:); (size(t1.data,2)-ROI.coords(2,:)); (size(t1.data,3)-ROI.coords(1,:))]), coordsAll');
  t1.data(indexAll(whichAfniVertex))=1;
  whichAfniVertex=nearpoints(coordsAll', double([ROI.coords(3,:); (size(t1.data,2)-ROI.coords(2,:)); (size(t1.data,3)-ROI.coords(1,:))]));
  t1.data(indexAllVista(whichAfniVertex))=1;    
  
  writeFileNifti(t1);
  %This part warps the result to template
  instr=sprintf('!3dAllineate -1Dmatrix_apply t1_skullStripped_nnn_at.Xat.1D -master %s/TT_N27_SurfVol_target_nnn.nii -interp linear -prefix %sInAnatomy_fff.nii.gz -input %s', templateDir, char(strcat('RoiLeftAllMap')), char(strcat('RoiLeftAllMap.nii.gz')));
  commitEval(instr);
  instr=sprintf('!3dNwarpApply -nwarp Qwarp_WARP+tlrc -source %sInAnatomy_fff.nii.gz -master %s/TT_N27_SurfVol_target_nnn.nii -prefix %sInAnatomy_warp_fff.nii.gz -ainterp linear', char(strcat('RoiLeftAllMap')), templateDir, char(strcat('RoiLeftAllMap')));
  commitEval(instr);
  instr=sprintf('!cp %sInAnatomy_warp_fff.nii.gz %s/%s/%sInAnatomy_mask_fff.nii.gz', char(strcat('RoiLeftAllMap')), templateDir, 'Map', char(strcat(subOrder{nSubj}, 'LeftAllMap')));
  commitEval(instr); 
  
  load(char(strcat(pathsRoi{nSubj}, '/Gray/ROIs/RightAllMap.mat')));
  t1=readFileNifti('Segmentation.nii.gz');
  t1.fname=char(strcat('RoiRightAllMap.nii.gz'));
  t1.data=zeros(size(t1.data));

  whichAfniVertex=nearpoints(double([ROI.coords(3,:); (size(t1.data,2)-ROI.coords(2,:)); (size(t1.data,3)-ROI.coords(1,:))]), coordsAll');
  t1.data(indexAll(whichAfniVertex))=1;
  whichAfniVertex=nearpoints(coordsAll', double([ROI.coords(3,:); (size(t1.data,2)-ROI.coords(2,:)); (size(t1.data,3)-ROI.coords(1,:))]));
  t1.data(indexAllVista(whichAfniVertex))=1;    
  
  writeFileNifti(t1);
  %This part warps the result to template
  instr=sprintf('!3dAllineate -1Dmatrix_apply t1_skullStripped_nnn_at.Xat.1D -master %s/TT_N27_SurfVol_target_nnn.nii -interp linear -prefix %sInAnatomy_fff.nii.gz -input %s', templateDir, char(strcat('RoiRightAllMap')), char(strcat('RoiRightAllMap.nii.gz')));
  commitEval(instr);
  instr=sprintf('!3dNwarpApply -nwarp Qwarp_WARP+tlrc -source %sInAnatomy_fff.nii.gz -master %s/TT_N27_SurfVol_target_nnn.nii -prefix %sInAnatomy_warp_fff.nii.gz -ainterp linear', char(strcat('RoiRightAllMap')), templateDir, char(strcat('RoiRightAllMap')));
  commitEval(instr);
  instr=sprintf('!cp %sInAnatomy_warp_fff.nii.gz %s/%s/%sInAnatomy_mask_fff.nii.gz', char(strcat('RoiRightAllMap')), templateDir, 'Map', char(strcat(subOrder{nSubj}, 'RightAllMap')));
  commitEval(instr); 
  
  
  okROIs=ones(length(hemispheres), length(mapNames), length(whichRoi));
  for nHemi=1:length(hemispheres)

      
      for nMap=1:length(mapNames)
          for nRoiType=2;%1:length(whichRoi)
              %This part makes ROIs into niftis, and can be commented if
              %already done
             try
                load(char(strcat(pathsRoi{nSubj}, '/Gray/ROIs/', hemispheres(nHemi), mapNames(nMap), whichRoi(nRoiType), '.mat')));
             catch
                 okROIs(nHemi, nMap,nRoiType)=0;
             end
             t1=readFileNifti('Segmentation.nii.gz');
             t1.fname=char(strcat('Roi', hemispheres(nHemi), mapNames(nMap), whichRoi(nRoiType), '.nii.gz'));
             t1.data=zeros(size(t1.data));
             if okROIs(nHemi, nMap,nRoiType)==1
                 if nHemi==1
                     %for m=1:size(ROI.coords,2)
                         whichAfniVertex=nearpoints(double([ROI.coords(3,:); (size(t1.data,2)-ROI.coords(2,:)); (size(t1.data,3)-ROI.coords(1,:))]), coordsLeft');
                         t1.data(indexLeft(whichAfniVertex))=1;
                         if nRoiType==1
                             whichAfniVertex=nearpoints(coordsLeft', double([ROI.coords(3,:); (size(t1.data,2)-ROI.coords(2,:)); (size(t1.data,3)-ROI.coords(1,:))]));
                             t1.data(indexLeft(whichAfniVertex))=1;
                         end
                     %end
                 else
                         whichAfniVertex=nearpoints(double([ROI.coords(3,:); (size(t1.data,2)-ROI.coords(2,:)); (size(t1.data,3)-ROI.coords(1,:))]), coordsRight');
                         t1.data(indexRight(whichAfniVertex))=1;
                         if nRoiType==1
                             whichAfniVertex=nearpoints(coordsRight', double([ROI.coords(3,:); (size(t1.data,2)-ROI.coords(2,:)); (size(t1.data,3)-ROI.coords(1,:))]));
                             t1.data(indexRight(whichAfniVertex))=1;
                         end
                 end
             end
             writeFileNifti(t1);
             
             %This part warps the result to template
             if nRoiType==1
                 instr=sprintf('!3dAllineate -1Dmatrix_apply t1_skullStripped_nnn_at.Xat.1D -master %s/TT_N27_SurfVol_target_nnn.nii -interp linear -prefix %sInAnatomy_fff.nii.gz -input %s', templateDir, char(strcat('Roi', hemispheres(nHemi), mapNames(nMap), whichRoi(nRoiType))), char(strcat('Roi', hemispheres(nHemi), mapNames(nMap), whichRoi(nRoiType), '.nii.gz')));
                 commitEval(instr);
                 instr=sprintf('!3dNwarpApply -nwarp Qwarp_WARP+tlrc -source %sInAnatomy_fff.nii.gz -master %s/TT_N27_SurfVol_target_nnn.nii -prefix %sInAnatomy_warp_fff.nii.gz -ainterp linear', char(strcat('Roi', hemispheres(nHemi), mapNames(nMap), whichRoi(nRoiType))), templateDir, char(strcat('Roi', hemispheres(nHemi), mapNames(nMap), whichRoi(nRoiType))));
                 commitEval(instr);
                 instr=sprintf('!cp %sInAnatomy_warp_fff.nii.gz %s/%s/%sInAnatomy_mask_fff.nii.gz', char(strcat('Roi', hemispheres(nHemi), mapNames(nMap), whichRoi(nRoiType))), templateDir, whichRoi(nRoiType), char(strcat(subOrder{nSubj}, hemispheres(nHemi), mapNames(nMap), whichRoi(nRoiType))));
                 commitEval(instr);
             else
                 instr=sprintf('!3dAllineate -1Dmatrix_apply t1_skullStripped_nnn_at.Xat.1D -master %s/TT_N27_SurfVol_target_nnn.nii -interp NN -prefix %sInAnatomy_fff.nii.gz -input %s', templateDir, char(strcat('Roi', hemispheres(nHemi), mapNames(nMap), whichRoi(nRoiType))), char(strcat('Roi', hemispheres(nHemi), mapNames(nMap), whichRoi(nRoiType), '.nii.gz')));
                 commitEval(instr);
                 instr=sprintf('!3dNwarpApply -nwarp Qwarp_WARP+tlrc -source %sInAnatomy_fff.nii.gz -master %s/TT_N27_SurfVol_target_nnn.nii -prefix %sInAnatomy_warp_fff.nii.gz -ainterp linear', char(strcat('Roi', hemispheres(nHemi), mapNames(nMap), whichRoi(nRoiType))), templateDir, char(strcat('Roi', hemispheres(nHemi), mapNames(nMap), whichRoi(nRoiType))));
                 commitEval(instr);
                 instr=sprintf('!cp %sInAnatomy_warp_fff.nii.gz %s/%s/%sInAnatomy_mask_fff.nii.gz', char(strcat('Roi', hemispheres(nHemi), mapNames(nMap), whichRoi(nRoiType))), templateDir, whichRoi(nRoiType), char(strcat(subOrder{nSubj}, hemispheres(nHemi), mapNames(nMap), whichRoi(nRoiType))));
                 commitEval(instr);
             end
          end
      end
  end
  



%   %!3dcalc -a roisVol_left_zzz.nii.gz -b roisVol_right_zzz.nii.gz -expr "step(a+b)" -prefix clusterVolume_fff.nii.gz;
%   instr=sprintf('!3dAllineate -1Dmatrix_apply t1_skullStripped_nnn_at.Xat.1D -master %s/TT_N27_SurfVol_target_nnn.nii -interp NN -prefix clusterInAnatomy_fff.nii.gz -input clusterVolume_fff.nii.gz', templateDir);
%   commitEval(instr);
%   instr=sprintf('!3dNwarpApply -nwarp Qwarp_WARP+tlrc -source clusterInAnatomy_fff.nii.gz -master %s/TT_N27_SurfVol_target_nnn.nii -prefix clusterInAnatomy_warp_fff.nii.gz -ainterp NN', templateDir);
%   commitEval(instr);
%   %instr <- sprintf('3dmerge -1blur_fwhm 1.0 -doall -prefix clusterInAnatomy_warp_blur_fff.nii.gz clusterInAnatomy_warp_fff.nii.gz'); system(instr);
%   %instr <- '3dcalc -a clusterInAnatomy_warp_fff.nii.gz -expr "step(a)" -prefix clusterInAnatomy_mask_fff.nii.gz'; system(instr);
%   %stringSplit <- strsplit( subjDir, '/')
%   %subjName <- stringSplit[[1]][ length( stringSplit[[1]] ) ]
%   instr=sprintf('!cp clusterInAnatomy_warp_fff.nii.gz %s/clusterInAnatomy_mask_%s_fff.nii.gz', templateDir, subjName); 
%   commitEval(instr);
end 

%Now we need to make a mrVista structure for the template.
%First, copy any existing session.
%Open the inplane window.
%Install the template anatomy
%Install the template segmentation. Keep all gray nodes.
mrVista inplane
installSegmentation(1, 1);
mrVista 3;
VOLUME{1}=makeGrayROI(VOLUME{1}); VOLUME{1}=refreshScreen(VOLUME{1},0);
VOLUME{1} = roiRestricttoLayer1(VOLUME{1},VOLUME{1}.selectedROI); 
VOLUME{1} = refreshScreen(VOLUME{1},0);

segName='segBen3.nii.gz';
T1=readFileNifti('T1.nii.gz');
seg=readFileNifti(segName);
segTmp=seg.data;
seg=T1;
seg.fname=segName;
seg.data=segTmp;
seg
figure; imagesc(seg.data(:,:,100))
writeFileNifti(seg);

segTmp(segTmp==6)=16;
segTmp(segTmp==5)=15;
for n=1:size(VOLUME{1}.ROIs(1).coords, 2)
    segTmp(VOLUME{1}.ROIs(1).coords(3,n), seg.dim(2)+1-VOLUME{1}.ROIs(1).coords(2,n), seg.dim(1)+1-VOLUME{1}.ROIs(1).coords(1,n))=segTmp(VOLUME{1}.ROIs(1).coords(3,n), seg.dim(2)+1-VOLUME{1}.ROIs(1).coords(2,n), seg.dim(1)+1-VOLUME{1}.ROIs(1).coords(1,n))-10;
end
%figure; imagesc(segTmp(:,:,100))
segTmp(segTmp>14)=1;
%figure; imagesc(segTmp(:,:,100))
seg=T1;
seg.data=segTmp;
seg
figure; imagesc(seg.data(:,:,100))
seg.fname='SegmentationLayer1.nii.gz';
writeFileNifti(seg);
%While you have this window open. prepare a box for putting model data in, and save it.
allTmp=0-ones(size(seg.data));
for n=1:size(VOLUME{1}.ROIs(1).coords, 2)
    allTmp(VOLUME{1}.ROIs(1).coords(3,n), seg.dim(2)+1-VOLUME{1}.ROIs(1).coords(2,n), seg.dim(1)+1-VOLUME{1}.ROIs(1).coords(1,n))=0;
end
figure; imagesc(allTmp(:,:,100))
save('allTmp.mat', 'allTmp')

%Now transform each ROI onto that template surface
segmentationVista=readFileNifti('SegmentationLayer1.nii.gz');

segVista=segmentationVista.data;
indexAllVista=find(segVista==6 | segVista==5);
[x y z]=ind2sub(size(segVista), indexLeftVista);
coordsLeftVista=[x,y,z];
[x y z]=ind2sub(size(segVista), indexRightVista);
coordsRightVista=[x,y,z];
[x y z]=ind2sub(size(segVista), indexAllVista);
coordsAllVista=[x,y,z];

%Get the niftis out of the resulting directories and re-warp to surface
InRoiThreshold=0.2; %A threshold of the proportion overlap beween the ROI and vertices in the warped surface
HowMuchExtraDistance=2; %Dilation of ROI to fill gaps, in mm
cd([templateDir, '/Map']);
!rm sumMaps.nii.gz
files=dir('*AllMap*.nii.gz');
for nFile=1:length(files)
    segmentationAfni=readFileNifti(files(nFile).name);
    segAppVol=segmentationAfni.data;
    indexAll=find(segAppVol>InRoiThreshold);
    [x y z]=ind2sub(size(segAppVol), indexAll);
    coordsAll=[x,y,z];
    t1=segmentationVista;
    t1.fname=char(strcat('GoodWarp_', files(nFile).name));
    
    t1.data=zeros(size(t1.data));
    [bestVistaVertex, bestDistances]=nearpoints(coordsAll', coordsAllVista');
    for nVertex=1:length(bestVistaVertex)
        allDists=sqrt(sum((coordsAllVista-coordsAll(nVertex,:)).^2, 2));
        passedVistaVertices=find(allDists<(sqrt(bestDistances(nVertex))+HowMuchExtraDistance));
        t1.data(indexAllVista(passedVistaVertices))=1;
    end
    writeFileNifti(t1);
end

%Now sum resulting niftis
cd([templateDir, '/Map']);
!rm sumMaps.nii.gz
%files=dir('GoodWarp_*AllMap*.nii.gz');
files=[dir('GoodWarp_*LeftT*.nii.gz'); dir('GoodWarp_*RightT*.nii.gz')];
sumMaps=readFileNifti(files(1).name);
sumMaps.fname='sumMaps.nii.gz';
for nFile=2:length(files)
    mapNii=readFileNifti(files(nFile).name);
    max(mapNii.data(:))
    sumMaps.data=sumMaps.data+mapNii.data;
end
writeFileNifti(sumMaps);

%Read the result into mrVista
sumMap=readFileNifti('sumMaps.nii.gz');

%Load an example model
load('/media/Storage3/WholeBrainNewSeg/AFNI_preproc/NellCombined/Gray/TimingSweeps (L1)/HrfFitFreeExponent/retModel-20191013-215916-Lin-2dOvalGaussian-DurationPeriod-DT0.5-maxValue-2-expIntensity-free-fFit-fFit-fFit-fFit.mat')

%Convert nifti data into coordinates list
sumData=zeros(1, size(VOLUME{1}.coords, 2));
for whichCoord=1:size(VOLUME{1}.coords, 2)
    sumData(whichCoord)=sumMap.data(VOLUME{1}.coords(3,whichCoord),  sumMap.dim(2)+1-VOLUME{1}.coords(2,whichCoord), sumMap.dim(1)+1-VOLUME{1}.coords(1,whichCoord));
end
sumData(sumData>16)=0;
model{1}.x0=sumData;
model{1}.y0=zeros(size(sumData));
model{1}.sigma.major=sumData;
model{1}.sigma.minor=sumData;
model{1}.sigma.theta=zeros(size(sumData));
model{1}.rawrss=ones(size(sumData)).*10;
model{1}.rss=10-sumData;

save('/media/Storage3/WholeBrainNewSeg/ForAFNIclustering/FromAlessioMatlab/TemplateMrVista/Gray/Original/nSubsAtlas.mat', 'model', 'params')

%To set the view parameters
VOLUME{1} = rmLoad(VOLUME{1}, 1, 'x0', 'map');
VOLUME{1} = setClipMode(VOLUME{1}, 'map', [0 8]);
VOLUME{1} = refreshScreen(VOLUME{1}, 1);

%Now for the map centers
%colors=["r", "g", "b", "c", "m", "y", "k", "w"]
segmentationVista=readFileNifti('SegmentationLayer1.nii.gz');

%segmentationVista=readFileNifti('Seg_2.nii.gz');
segVista=segmentationVista.data;
indexAllVista=find(segVista==6 | segVista==5);
[x y z]=ind2sub(size(segVista), indexLeftVista);
coordsLeftVista=[x,y,z];
[x y z]=ind2sub(size(segVista), indexRightVista);
coordsRightVista=[x,y,z];
[x y z]=ind2sub(size(segVista), indexAllVista);
coordsAllVista=[x,y,z];

colours={'r', 'g', 'b', 'c', 'm', 'y', 'k', 'w', 'r', 'g', 'b', 'c', 'm', 'y', 'k', 'w',};
VOLUME{1}=deleteAllROIs(VOLUME{1}); VOLUME{1}=refreshScreen(VOLUME{1},0);
for whichMap=1:length(mapNames)
    cd([templateDir, '/Center']);
    files=dir(['*Right', char(mapNames(whichMap)),'Center*.nii.gz']);
%     cd('/media/Storage3/WholeBrainNewSeg/ForAFNIclustering/FromAlessioMatlab/TemplateMrVista')
%     [VOLUME{1}, ok] = loadROI(VOLUME{1}, 'gray-Layer1', 1,colours{whichMap},0, 1);
%     VOLUME{1}.ROIs(whichMap).coords=[];
%     %VOLUME{1}.ROIs(whichMap).color=colours{whichMap}
%     VOLUME{1}.ROIs(whichMap).name=char(mapNames(whichMap))
    cd([templateDir, '/Center']);
    eval(sprintf('CenterCoordsRight%s=[];', mapNames(whichMap)))
    
    for nFile=1:length(files)
        segmentationAfni=readFileNifti(files(nFile).name);
        segAppVol=segmentationAfni.data;
        [~, indexAll]=max(segAppVol(:));
        [x y z]=ind2sub(size(segAppVol), indexAll);
        coordsAll=[x,y,z]
        t1=segmentationVista;
        t1.fname=char(strcat('GoodWarp_', files(nFile).name));
        t1.data=zeros(size(t1.data));
        
        if ~isempty(coordsAll);
            [bestVistaVertex, bestDistances]=nearpoints(coordsAll', coordsAllVista');
            coords=coordsAllVista((bestVistaVertex), :);
            eval(sprintf('CenterCoordsRight%s=[CenterCoordsRight%s; coords];', mapNames(whichMap), mapNames(whichMap)))
            t1.data(indexAllVista(bestVistaVertex))=whichMap;
            
%             VOLUME{1}.ROIs(whichMap).coords=[VOLUME{1}.ROIs(whichMap).coords, [seg.dim(1)+1-coords(3) seg.dim(2)+1-coords(2) coords(1)]'];
        end
        writeFileNifti(t1);
    end
end


cd([templateDir, '/Center']);
!rm sumMaps.nii.gz
files=dir('GoodWarp_*.nii.gz');
sumMaps=readFileNifti(files(1).name);
sumMaps.fname='sumMaps.nii.gz';
%sumMaps.data=zeros(size(sumMaps.data));
for nFile=2:length(files)
    mapNii=readFileNifti(files(nFile).name);
    sumMaps.data=sumMaps.data+mapNii.data;
end
writeFileNifti(sumMaps);

%Read the result into mrVista
sumMap=readFileNifti('sumMaps.nii.gz');

%Load an example model
load('/media/Storage3/WholeBrainNewSeg/AFNI_preproc/NellCombined/Gray/TimingSweeps (L1)/HrfFitFreeExponent/retModel-20191013-215916-Lin-2dOvalGaussian-DurationPeriod-DT0.5-maxValue-2-expIntensity-free-fFit-fFit-fFit-fFit.mat')

%Convert nifti data into coordinates list
sumData=zeros(1, size(VOLUME{1}.coords, 2));
for whichCoord=1:size(VOLUME{1}.coords, 2)
    sumData(whichCoord)=sumMap.data(VOLUME{1}.coords(3,whichCoord),  sumMap.dim(2)+1-VOLUME{1}.coords(2,whichCoord), sumMap.dim(1)+1-VOLUME{1}.coords(1,whichCoord));
end
model{1}.x0=sumData;
model{1}.sigma.major=sumData;
model{1}.sigma.minor=sumData;
model{1}.sigma.theta=zeros(size(sumData));
model{1}.rawrss=ones(size(sumData)).*10;
model{1}.rss=1-(sumData>0);
save('Gray/Original/CentersAtlas.mat', 'model', 'params')



segmentationAfni=readFileNifti('HarvGray-Layer1InAnatomy_mask_fff.nii.gz');
segAppVol=segmentationAfni.data;
indexAll=find(segAppVol==1);
[x y z]=ind2sub(size(segAppVol), indexAll);
coordsAll=[x,y,z];



 t1=segmentationVista;
 t1.fname=char(strcat('RoiGray-Layer1_GoodWarp.nii.gz'));
 t1.data=zeros(size(t1.data));
 whichVistaVertex=nearpoints(coordsAll', coordsAllVista');
 t1.data(indexAllVista(whichVistaVertex))=1;
 
 
 %Try this again later
%  t1b=t1;
%  t1b.data=zeros(size(t1b.data));
%  whichAfniVertex=nearpointsBackwards(coordsAll', coordsAll', coordsAllVista');
%  t1b.data(indexAllVista(whichAfniVertex))=1;
 


%Open a new mrVista window and load the model. 
%Load the ROI in which the model was run.
%Now, to pull the data out
load('allTmp.mat')
T1=readFileNifti('t1_1mm.nii');
ves=rmGet(VOLUME{1}.rm.retinotopyModels{1}, 've');
x=rmGet(VOLUME{1}.rm.retinotopyModels{1}, 'x');
y=rmGet(VOLUME{1}.rm.retinotopyModels{1}, 'y');

smoothing=5;
if smoothing>0
    [tmp, ii]=intersectCols(VOLUME{1}.coords, VOLUME{1}.ROIs(end).coords);
    mask=zeros(size(VOLUME{1}.coords));
    mask(ii)=1;
    [ves conMat]=dhkGraySmooth(VOLUME{1},ves,[smoothing 0.6], [], mask);
    [x conMat]=dhkGraySmooth(VOLUME{1},x,[smoothing 0.6], conMat, mask);
    [y conMat]=dhkGraySmooth(VOLUME{1},y,[smoothing 0.6], conMat, mask);
end











%Get the niftis out of the resulting directories and sum
cd([templateDir, '/Map']);
!rm sumMaps.nii.gz
files=dir('*.nii.gz');
sumMaps=readFileNifti(files(1).name);
sumMaps.fname='sumMaps.nii.gz';
for nFile=2:length(files)
    mapNii=readFileNifti(files(nFile).name);
    sumMaps.data=sumMaps.data+mapNii.data;
end
writeFileNifti(sumMaps);
    


cd(templateDir) 
for nSubj=1:length(subjDirectories)
  subjDir=subjDirectories{nSubj};
  stringSplit=split(subjDir, '/');
  subjName=char(stringSplit(end));
  if nSubj==1
    instr=sprintf('!3dMean -prefix averageData_fff.nii.gz modelTemplate_%s_fff.nii.gz', subjName);
  else
    instr=sprintf('!%s modelTemplate_%s_fff.nii.gz', instr, subjName)
  end
end
commitEval(instr);
%this creates the volume superimposing single participant clustering
%results (Ben will probably do this another way)
!3dmask_tool -input clusterInAnatomy* -dilate_input 2 -1 -prefix mask_counts_fff.nii.gz -count'; system( instr ); 


%Once clusters are determined, load the cluster nifti, then:
leftClust=readFileNifti('roiLeftLayer1_Sub5_zzz.nii.gz');
rightClust=readFileNifti('roiRightLayer1_Sub5_zzz.nii.gz');
colours={'r', 'g', 'b', 'c', 'm', 'y', 'k', 'w', 'r', 'g', 'b', 'c', 'm', 'y', 'k', 'w',};
for whichRoi=1:max(leftClust.data(:))
    [VOLUME{1}, ok] = loadROI(VOLUME{1}, 'gray-Layer1', 1,colours{whichRoi},0, 1);
    [x,y,z]=ind2sub(size(leftClust.data), find(leftClust.data==whichRoi));
    VOLUME{1}.ROIs(end).coords=[leftClust.dim(3)+1-z'; leftClust.dim(2)+1-y'; x'];
    VOLUME{1}.ROIs(end).name=strcat('cluster', num2str(whichRoi), 'Left');
end
for whichRoi=1:max(rightClust.data(:))
    [VOLUME{1}, ok] = loadROI(VOLUME{1}, 'gray-Layer1', 1,colours{whichRoi},0, 1);
    [x,y,z]=ind2sub(size(rightClust.data), find(rightClust.data==whichRoi));
    VOLUME{1}.ROIs(end).coords=[rightClust.dim(3)+1-z'; rightClust.dim(2)+1-y'; x'];
    VOLUME{1}.ROIs(end).name=strcat('cluster', num2str(whichRoi), 'Right');
end
for whichROI=1:size(VOLUME{1}.ROIs, 2)
    VOLUME{1}.selectedROI=whichROI;
    VOLUME{1}=mrv_dilateCurrentROI(VOLUME{1});   
end

VOLUME{1}.ROIs=VOLUME{1}.ROIs((whichROI+1):size(VOLUME{1}.ROIs, 2));
VOLUME{1}=refreshScreen(VOLUME{1},0);

for whichROI=1:size(VOLUME{1}.ROIs, 2)
    VOLUME{1}.selectedROI=whichROI;
    VOLUME{1} = roiRestricttoLayer1(VOLUME{1});
end
VOLUME{1}=refreshScreen(VOLUME{1},0);


%VF mapping with css
paths{1}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/S08_VFM_Num';
paths{2}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/S09_VFM_Num';
paths{3}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/S10_VFM_Num';
paths{4}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/S11_VFM_Num';
paths{5}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/S12_VFM_Num';
paths{6}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/S13_VFM_Num';
paths{7}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/Harv_VFM_Num';
paths{8}='/media/Storage3/WholeBrainNewSeg/AFNI_preproc/Nell_VFM_Num';

whichSubs=[6];

for thisSub=1:length(whichSubs)
    cd(paths{whichSubs(thisSub)})
    startup vo;
    mrVista 3;
    startup vm;
    matfilename = 'retModel-FreeHRF_NoDC_OneG_CSS';
    tmpv=rmMain([1 length(dataTYPES)],'gray-Layer1',4,...
        'prf model', {'css'},...
        'matfilename',matfilename,...
        'coarsetofine',false,...
        'coarseDecimate',0,...
        'datadrivendc',false);
end

%For negative betas, change line 65 of rmGridFit_oneGaussianNonlinear.m and
%line 49 of rmModelSearchFit_oneGaussianNonlinear and line 105 of rmSearchFit_oneGaussianNonlinear
for thisSub=1:length(whichSubs)
    cd(paths{whichSubs(thisSub)})
    startup vo;
    mrVista 3;
    startup vm;
    matfilename = 'retModel-FreeHRF_NoDC_OneG_CSS_NegativeBetas';
    tmpv=rmMain([1 length(dataTYPES)],'gray-Layer1',4,...
        'prf model', {'css'},...
        'matfilename',matfilename,...
        'coarsetofine',false,...
        'coarseDecimate',0,...
        'datadrivendc',false);
end

%When something goes wrong with that:
folderName=[pwd '/Gray/' dataTYPES(14).name];
modelFile='retModel-FreeHRF_NoDC_OneG_CSS_NegativeBetas-gFit.mat'
rmMainPostGrid([1 14],'gray-Layer1',4, [folderName '/' modelFile]);


