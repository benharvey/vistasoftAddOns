function RunCandidateTimingModels(paths, whichSubs, combinedDT, maxCores, ParamsPathMix, ParamsPathSessions, FullModelSet)
%Fits the main candidate timing selectivity models from Harvey et al 2020
%If FullModelSet input is false, only runs the best fitting model
%Further candidate models can also be run (and were), but each takes a long
%time and this is the most informative set.

%BMH, 01/2020

for thisSub=whichSubs
   cd(paths{thisSub})
   mrVista 3;
   
   %The code to run for each
   if thisSub<=6
       load(ParamsPathMix)
       setAllRetParams(paramsDurationPeriod, combinedDT);
       rmRunDuration2dOval(VOLUME{1},combinedDT, 'gray-Layer1',5,{'1g'},[],[],[],0,'free'); %This is the best fitting model. Setting the 4th input arguement to '5' fits the HRF parameters too.
       
       if FullModelSet
           setAllRetParams(paramsDurationPeriod, combinedDT);
           rmRunDuration2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],maxCores,0);
           setAllRetParams(paramsDurationPeriod, combinedDT);
           rmRunLogDuration2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0,'free');
           for whichModel=[12 14 15 16]
               setAllRetParams(paramsDurationPeriod, combinedDT);
               rmRunDurFreq2d(VOLUME{1},combinedDT, 'gray-Layer1',whichModel,{'1g'},[],[],[],0);
           end
           
           setAllRetParams(paramsOccupancyPeriod, combinedDT);
           rmRunOccupancyPeriod2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0,'free');
           for whichModel=[12 14]
               setAllRetParams(paramsOccupancyPeriod, combinedDT);
               rmRunOccupancyFreq2d(VOLUME{1},combinedDT, 'gray-Layer1',whichModel,{'1g'},[],[],[],0);
           end
           
           setAllRetParams(paramsOnTimeOffTime, combinedDT);
           rmRunOnTimeOffTime2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],0,'free');
           
           setAllRetParams(paramsTemporalFrequency, combinedDT);
           rmRunTemporalFreq1d(VOLUME{1},combinedDT, 'gray-Layer1',14,{'1g'},[],[],[],0);
           
           %For temporal frequency models run per TR (rather than per event)
           load(ParamsPathMix)
           for n=1:4
               paramsTemporalFrequency(n).paramsFile=[paramsTemporalFrequency(n).paramsFile(1:(end-9)), 'TF', paramsTemporalFrequency(n).paramsFile((end-5):end)];
           end
           setAllRetParams(paramsTemporalFrequency, combinedDT);
           rmRunTemporalFreq1d(VOLUME{1},combinedDT, 'gray-Layer1',14,{'1g'},[],sprintf('retModel-%s-Lin-1dGaussianXnoY-TemporalFreq-%s',datestr(now,'yyyymmdd-HHMMSS'), num2str(20,2)),[],0);
           load(ParamsPathMix)
       end
       
   else
       load(ParamsPathSessions)
       setAllRetParams(paramsDurationPeriod, combinedDT);
       rmRunDuration2dOval(VOLUME{1},combinedDT, 'gray-Layer1',5,{'1g'},[],[],[],1,'free'); %This is the best fitting model. Setting the 4th input arguement to '5' fits the HRF parameters too.
       
       if FullModelSet
           setAllRetParams(paramsDurationPeriod, combinedDT);
           rmRunDuration2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],maxCores,1);
           
           setAllRetParams(paramsDurationPeriod, combinedDT);
           rmRunLogDuration2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],1,'free');
           for whichModel=[12 14 15 16]
               setAllRetParams(paramsDurationPeriod, combinedDT);
               rmRunDurFreq2d(VOLUME{1},combinedDT, 'gray-Layer1',whichModel,{'1g'},[],[],[],1);
           end
           
           setAllRetParams(paramsOccupancyPeriod, combinedDT);
           rmRunOccupancyPeriod2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],1,'free');
           for whichModel=[12 14]
               setAllRetParams(paramsOccupancyPeriod, combinedDT);
               rmRunOccupancyFreq2d(VOLUME{1},combinedDT, 'gray-Layer1',whichModel,{'1g'},[],[],[],1);
           end
           
           setAllRetParams(paramsOnTimeOffTime, combinedDT);
           rmRunOnTimeOffTime2dOval(VOLUME{1},combinedDT, 'gray-Layer1',4,{'1g'},[],[],[],1,'free');
           
           setAllRetParams(paramsTemporalFrequency, combinedDT);
           rmRunTemporalFreq1d(VOLUME{1},combinedDT, 'gray-Layer1',14,{'1g'},[],[],[],1);
           
           %For temporal frequency models run per TR (rather than per event)
           load(ParamsPathSessions)
           for n=1:4
               paramsTemporalFrequency(n).paramsFile=[paramsTemporalFrequency(n).paramsFile(1:(end-9)), 'TF', paramsTemporalFrequency(n).paramsFile((end-5):end)];
           end
           setAllRetParams(paramsTemporalFrequency, combinedDT);
           rmRunTemporalFreq1d(VOLUME{1},combinedDT, 'gray-Layer1',14,{'1g'},[],sprintf('retModel-%s-Lin-1dGaussianXnoY-TemporalFreq-%s',datestr(now,'yyyymmdd-HHMMSS'), num2str(20,2)),[],1);
           load(ParamsPathSessions)
       end
   end
   close(1); mrvCleanWorkspace; 
end
end

