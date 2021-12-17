function view = makeCombinedDTWholeBrainNumber(view)
mrGlobals
% viewSet(view, 'curdt', firstdt);
% saveSession;
% savePrefs(view);
duplicateDataType(view,'NumberConditionsL1');
curDataType = length(dataTYPES);%viewGet(view,'curDataType');
dataTYPES(curDataType).scanParams(2:4)=dataTYPES(curDataType).scanParams(1);
dataTYPES(curDataType).blockedAnalysisParams(2:4)=dataTYPES(curDataType).blockedAnalysisParams(1);
dataTYPES(curDataType).eventAnalysisParams(2:4)=dataTYPES(curDataType).eventAnalysisParams(1);
%dataTYPES(curDataType).retinotopyModelParams(2:4)=dataTYPES(curDataType).retinotopyModelParams(1);
dataTYPES(curDataType).scanParams(1).annotation='area';
dataTYPES(curDataType).scanParams(2).annotation='size';
dataTYPES(curDataType).scanParams(3).annotation='circ';
dataTYPES(curDataType).scanParams(4).annotation='dense';
saveSession;
savePrefs(view);

!cp -r Gray/AveragesSize/TSeries/Scan1 Gray/NumberConditions/TSeries/Scan2
!cp -r Gray/AveragesCirc/TSeries/Scan1 Gray/NumberConditions/TSeries/Scan3
!cp -r Gray/AveragesDense/TSeries/Scan1 Gray/NumberConditions/TSeries/Scan4
end

