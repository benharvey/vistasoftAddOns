function [view, connectionMatrix] = collapseOverLayersTseries(view, connectionMatrix, dosave, ROI, dataType)
% collapseOverLayersTseries - Function to collapse data from all layers 
% onto layer 1 of the gray structure, weighted by varience explained.
% Reduces noise and makes a more 2-dimensional structure.
%
% If dosave=true (default), the new data is saved in a new datatype.
%
% BMH & SOD wrote it

if exist('dataType','var') && ~isempty(dataType)
    view=selectDataType(view,dataType);
end
if notDefined('connectionMatrix') || isempty(connectionMatrix)
    % for all layers
    connectionMatrix=mrmCollapseToLayer1(view.nodes, view.mmPerVox, view.edges);
end
if notDefined('dosave')
    dosave = true;
end

% if dosave duplicate datatype instead of overwriting
if dosave
    % get name
    dtnum  = viewGet(view,'curdt');
    dtname = viewGet(view,'dtname');
    % add extension to identify it as collapsed over layers
    dtname = [dtname ' (L1)'];
    % duplicate
    duplicateDataType(view,dtname);
    mrGlobals;
    view=viewSet(view,'curdt',length(dataTYPES));
end

% actual routine
nScans = viewGet(view,'numScans',view.curDataType);
layer1Nodes=find(connectionMatrix(1,:)>0);
for scans=1:nScans
    
    % load tSeries
    view.tSeries = loadtSeries(view, scans, 1);
    withData = isfinite(sum(view.tSeries));
    
    % fill with NaNs
    tSeriesOut=NaN(size(view.tSeries));
    
    % Determine which nodes are in ROI (if ROI is passed in)
    if exist('ROI', 'var') && ~isempty(ROI)
        allCrds  = viewGet(view,'coords');
        [tmp, RoiCrds] = intersectCols(allCrds,ROI.coords);
        inRoi=zeros(size(withData));
        inRoi(RoiCrds)=1;
        inRoi=logical(inRoi);
    end
    
    for node=layer1Nodes
        % add columns
        column=connectionMatrix(connectionMatrix(:,node)>0, node);
        column=column(withData(column));
        if exist('RoiCrds', 'var')
            column=column(inRoi(column));
        end

        % to prevent devisions by zero due to empty column
        if ~isempty(column)
            tSeriesOut(:,node) = mean(view.tSeries(:,column),2);
        end
    end
    view.tSeries=tSeriesOut;
    
    % save
    if dosave
        savetSeries(view.tSeries,view,scans,1)
    end
end