function CrossValidateCandidateTimingModels(paths, whichSubs, combinedDT, FinalModelFolder)%, ParamsPathMix, ParamsPathSessions)
%Cross validates a set of timing selectivity pRF models in a specified
%folder

%BMH, 01/2020

for thisSub=whichSubs
   cd(paths{thisSub})
   mrVista 3;
   
   %Unnecessary
%    if thisSub<=6
%        load(ParamsPathMix)
%    else
%        load(ParamsPathSessions)
%    end
    
   allXvalDTs=combinedDT(2:3);
   
   
   for n=1:length(allXvalDTs)
       files=dir(['Gray/' dataTYPES(allXvalDTs(n)).name, '/', FinalModelFolder, '/', '*.mat']);
       thisPath=['Gray/' dataTYPES(allXvalDTs(n)).name, '/', FinalModelFolder, '/'];
       otherPath=['Gray/' dataTYPES(allXvalDTs(3-n)).name, '/', FinalModelFolder, '/'];
       eval(['!mkdir ',  '"',otherPath, 'xval"']);
       eval(['!mkdir ',  '"',otherPath, 'xvalRefit"']);
       
       for whichFile=1:length(files)
           eval(['!cp ', '"', thisPath, files(whichFile).name, '" "', otherPath, 'xval/xval-', files(whichFile).name, '"']);
       end
   end
   parpool(2)
   parfor whichDT=1:length(allXvalDTs)
       folderName=[pwd '/Gray/' dataTYPES(allXvalDTs(whichDT)).name, '/', FinalModelFolder, '/xval'];
       modelFiles=dir([folderName '/*.mat']);
       for whichModel=1:length(modelFiles)
           rmMainPostSearch([1 allXvalDTs(whichDT)],'gray-Layer1',4, [folderName '/' modelFiles(whichModel).name]);
       end
   end
   delete(gcp('nocreate'))

end

