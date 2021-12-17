function VOLUME = NiftiClustersToVistaROIs(VOLUME, ClusterNifti, LeftOrRight)
colours = {'r', 'g', 'b', 'c', 'm', 'y', 'k', 'w', 'r', 'g', 'b', 'c', 'm', 'y', 'k', 'w',};

VOLUME{1}=deleteAllROIs(VOLUME{1}); VOLUME{1}=refreshScreen(VOLUME{1},0);

%Pull ROIs out of nifti
for whichRoi = 1:max(ClusterNifti.data(:))
    [VOLUME{1}, ok] = loadROI(VOLUME{1}, 'gray-Layer1', 1,colours{whichRoi},0, 1);
    [x,y,z] = ind2sub(size(ClusterNifti.data), find(ClusterNifti.data==whichRoi));
    VOLUME{1}.ROIs(end).coords = [ClusterNifti.dim(3)+1-z'; ClusterNifti.dim(2)+1-y'; x'];
    VOLUME{1}.ROIs(end).name = strcat('cluster', num2str(whichRoi), LeftOrRight);
end
%Dilate these ROIs to fill gaps
for whichROI = 1:size(VOLUME{1}.ROIs, 2)
    VOLUME{1}.selectedROI = whichROI;
    VOLUME{1} = mrv_dilateCurrentROI(VOLUME{1});   
end

%Delete undilated ROIs
VOLUME{1}.ROIs = VOLUME{1}.ROIs((whichROI+1):size(VOLUME{1}.ROIs, 2));
VOLUME{1} = refreshScreen(VOLUME{1},0);

%Restrict to layer 1
for whichROI = 1:size(VOLUME{1}.ROIs, 2)
    VOLUME{1}.selectedROI = whichROI;
    VOLUME{1} = roiRestricttoLayer1(VOLUME{1});
end
VOLUME{1} = refreshScreen(VOLUME{1},0);
end