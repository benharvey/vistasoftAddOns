function VOLUME = VistaSurfaceToNifti(T1name, segName, seg_path, numGrayLayers)

%BEWARE: PART OF THIS PROCESS WILL OVERWRITE THE DATA AND MODELS IN YOUR
%CURRENT FOLDER. DO THIS ON A COPY OF THE MRVISTA SESSION YOU ARE USING.

%Saves layer 1 mrVista surfaces as niftis for clustering
%01/2020 BMH

mrVista inplane
installSegmentation(1, 1, seg_path, numGrayLayers); % No, this really is broken.
VOLUME{1} = mrVista('3');
VOLUME{1} = makeGrayROI(VOLUME{1});
VOLUME{1} = refreshScreen(VOLUME{1},0);
VOLUME{1} = roiRestricttoLayer1(VOLUME{1},VOLUME{1}.selectedROI); 
VOLUME{1} = refreshScreen(VOLUME{1},0);

T1 = readFileNifti(T1name);
seg = readFileNifti(segName);
segTmp = seg.data;
seg = T1;
seg.fname = segName;
seg.data = segTmp;
seg %#ok<NOPRT>
figure; imagesc(seg.data(:,:,100))
writeFileNifti(seg);

segTmp(segTmp==6) = 16;
segTmp(segTmp==5) = 15;
for n = 1:size(VOLUME{1}.ROIs(1).coords, 2)
    segTmp(VOLUME{1}.ROIs(1).coords(3,n), seg.dim(2)+1-VOLUME{1}.ROIs(1).coords(2,n), seg.dim(1)+1-VOLUME{1}.ROIs(1).coords(1,n)) = segTmp(VOLUME{1}.ROIs(1).coords(3,n), seg.dim(2)+1-VOLUME{1}.ROIs(1).coords(2,n), seg.dim(1)+1-VOLUME{1}.ROIs(1).coords(1,n))-10;
end
%figure; imagesc(segTmp(:,:,100))
segTmp(segTmp>14) = 1;
%figure; imagesc(segTmp(:,:,100))
seg = T1;
seg.data = segTmp;
seg %#ok<NOPRT>
figure; imagesc(seg.data(:,:,100))
seg.fname = 'SegmentationLayer1.nii.gz';
writeFileNifti(seg);
%While you have this window open. prepare a box for putting model data in, and save it.
allTmp = 0-ones(size(seg.data));
for n = 1:size(VOLUME{1}.ROIs(1).coords, 2)
    allTmp(VOLUME{1}.ROIs(1).coords(3,n), seg.dim(2)+1-VOLUME{1}.ROIs(1).coords(2,n), seg.dim(1)+1-VOLUME{1}.ROIs(1).coords(1,n)) = 0;
end
figure; imagesc(allTmp(:,:,100))
save('allTmp.mat', 'allTmp')

end