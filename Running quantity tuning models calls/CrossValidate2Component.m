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
whichSubs=6:8 % 1:8; TODO change back
%Fitting the improved model
for thisSub=1:length(whichSubs)
    cd(paths{whichSubs(thisSub)})
    mrVista 3;
    
    %The code to run for each
    if whichSubs(thisSub)<=6
        load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsMix.mat')
        for whichModel=[15]
            setAllRetParams(paramsDurationPeriod, combinedDT);
            rmRunDurFreq2d(VOLUME{1},combinedDT, 'gray-Layer1',whichModel,{'1g'},[],[],[],0);
        end
    else
        load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsSessions.mat')
        for whichModel=[15]
            setAllRetParams(paramsDurationPeriod, combinedDT);
            rmRunDurFreq2d(VOLUME{1},combinedDT, 'gray-Layer1',whichModel,{'1g'},[],[],[],1);
        end
    end
    close(1); mrvCleanWorkspace;
    
end

%Fitting with compressive duration
for thisSub=1:length(whichSubs)
    cd(paths{whichSubs(thisSub)})
    mrVista 3;
    
    %The code to run for each
    if whichSubs(thisSub)<=6
        load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsMix.mat')
        for whichModel=[13]
            setAllRetParams(paramsDurationPeriod, combinedDT);
            rmRunDurFreq2d(VOLUME{1},combinedDT, 'gray-Layer1',whichModel,{'1g'},[],[],[],0);
        end
    else
        load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsSessions.mat')
        for whichModel=[13]
            setAllRetParams(paramsDurationPeriod, combinedDT);
            rmRunDurFreq2d(VOLUME{1},combinedDT, 'gray-Layer1',whichModel,{'1g'},[],[],[],1);
        end
    end
    close(1); mrvCleanWorkspace;
    
end

%Cross validating the resulting model. Move the models around first
for thisSub=whichSubs
   cd(paths{thisSub})
   mrVista 3;
   combinedDT=5:7;

   %The code to run for each
   if thisSub<=6
       load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsMix.mat')
   else
       load('/media/Storage3/WholeBrainNewSeg/TimingModelParamsSessions.mat')
   end
   allXvalDTs=combinedDT(2:3);
   
   
   
   %Monotonic component models
    for n=1:length(allXvalDTs)
       files=dir(['Gray/' dataTYPES(allXvalDTs(n)).name, '/Monotonic/', '*XcompressiveYNoNormOccupancy*.mat']);
       thisPath=['Gray/' dataTYPES(allXvalDTs(n)).name, '/Monotonic/'];
       otherPath=['Gray/' dataTYPES(allXvalDTs(3-n)).name, '/Monotonic/'];
%        eval(['!mkdir ',  '"',otherPath, 'xval"']);
%        eval(['!mkdir ',  '"',otherPath, 'xvalRefit"']);
       
       for whichFile=1:length(files)
           eval(['!cp ', '"', thisPath, files(whichFile).name, '" "', otherPath, 'xval/xval-', files(whichFile).name, '"']);
       end
   end
   
   
%  
   for whichDT=1:length(allXvalDTs)
       folderName=[pwd '/Gray/' dataTYPES(allXvalDTs(whichDT)).name, '/Monotonic/xval'];
       modelFiles=dir([folderName '/*XcompressiveYNoNormOccupancy*']);
       for whichModel=1:length(modelFiles)
           rmMainPostSearch([1 allXvalDTs(whichDT)],'gray-Layer1',4, [folderName '/' modelFiles(whichModel).name]);
       end
   end
   delete(gcp('nocreate'))
   
    close(1); mrvCleanWorkspace;
   
end