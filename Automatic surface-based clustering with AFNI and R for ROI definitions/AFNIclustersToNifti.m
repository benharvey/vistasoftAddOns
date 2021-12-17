function AFNIclustersToNifti(selSizes, selectedQuantile, selectedVarExp, actualVE)

!rm *0.1D *.niml.dset roisVol_left_zzz.nii.gz roisVol_right_zzz.nii.gz *_zzz.nii.gz *_zzz.nii.gz
  %load model data for each hemisphere and get varExp distribution (quantiles)
  modelLeft = readtable('modelResults_ttt_interp_surfaces_folder_left/boundary00_sm_modelResults_ttt_surf.1D.dset.txt');
  modelLeft = table2array(modelLeft);
  modelRight = readtable('modelResults_ttt_interp_surfaces_folder_right/boundary00_sm_modelResults_ttt_surf.1D.dset.txt');
  modelRight = table2array(modelRight);
  varExp_nodeLeft = modelLeft((find(modelLeft(:,9)>0.001)), 9);
  varExp_nodeRight = modelRight((find(modelRight(:,9)>0.001)) ,9);
  selectedSize = selSizes;
  sizeThrLoop = selectedSize;
  varQuantThrLoopLeft = quantile(varExp_nodeLeft, selectedQuantile);
  varQuantThrLoopRight = quantile(varExp_nodeRight,  selectedQuantile);
  if exist('actualVE', 'var') && actualVE == 1
    varQuantThrLoopLeft = selectedVarExp;
    varQuantThrLoopRight = selectedVarExp;
  end
  %cluster the surface based on variance explained
  disp(sprintf('size: %s, varExp (left): %s, varExp (right): %s', num2str(selectedSize), num2str(round(varQuantThrLoopLeft, -3)), num2str(round(varQuantThrLoopRight, -3)))) %#ok<DSPS>
  instr = sprintf( '!SurfClust -spec surfaces_folder_left/spec.surfaces.smoothed -surf_A surfaces_folder_left/boundary00_sm.1D.coord  -input modelResults_ttt_interp_surfaces_folder_left/boundary00_sm_modelResults_ttt_surf.1D.dset 8 -rmm -1 -amm2 %s -thresh_col 8 -thresh %s -out_roidset -prefix left_clus_roi_', num2str(sizeThrLoop), num2str(varQuantThrLoopLeft) );   
  eval(instr);
  instr = sprintf( '!SurfClust -spec surfaces_folder_right/spec.surfaces.smoothed -surf_A surfaces_folder_right/boundary00_sm.1D.coord  -input modelResults_ttt_interp_surfaces_folder_right/boundary00_sm_modelResults_ttt_surf.1D.dset 8 -rmm -1 -amm2 %s -thresh_col 8 -thresh %s -out_roidset -prefix right_clus_roi_', num2str(sizeThrLoop), num2str(varQuantThrLoopRight) ); 
  eval(instr);
  
  
  instr = sprintf('!3dSurf2Vol -spec surfaces_folder_left/spec.surfaces.smoothed -surf_A surfaces_folder_left/boundary00_sm.1D.coord -sdata left_clus_roi__ClstMsk_e1_a%s.0.niml.dset -grid_parent SegmentationLayer1_ttt.nii.gz -sv SegmentationLayer1_ttt.nii.gz -map_func ave -prefix roisVol_left_zzz.nii.gz', num2str(selSizes));
  eval(instr)
%   !3dSurf2Vol -spec surfaces_folder_left/spec.surfaces.smoothed -surf_A surfaces_folder_left/boundary00_sm.1D.coord -sdata left_clus_roi__ClstMsk_e1_a86.0.niml.dset -grid_parent SegmentationLayer1_ttt.nii.gz -sv SegmentationLayer1_ttt.nii.gz -map_func ave -prefix roisVol_left_zzz.nii.gz;
  instr = sprintf('!3dSurf2Vol -spec surfaces_folder_right/spec.surfaces.smoothed -surf_A surfaces_folder_right/boundary00_sm.1D.coord -sdata right_clus_roi__ClstMsk_e1_a%s.0.niml.dset -grid_parent SegmentationLayer1_ttt.nii.gz -sv SegmentationLayer1_ttt.nii.gz -map_func ave -prefix roisVol_right_zzz.nii.gz', num2str(selSizes));
  eval(instr)
%   !3dSurf2Vol -spec surfaces_folder_right/spec.surfaces.smoothed -surf_A surfaces_folder_right/boundary00_sm.1D.coord -sdata right_clus_roi__ClstMsk_e1_a86.0.niml.dset -grid_parent SegmentationLayer1_ttt.nii.gz -sv SegmentationLayer1_ttt.nii.gz -map_func ave -prefix roisVol_right_zzz.nii.gz;
  roisLeftFile = readFileNifti('roisVol_left_zzz.nii.gz');
  roisRightFile = readFileNifti('roisVol_right_zzz.nii.gz');
  segmentationFile = readFileNifti('SegmentationLayer1_ttt.nii.gz');
  
  roisLeftVolume = roisLeftFile.data;
  roisRightVolume = roisRightFile.data;
  segAppVol = segmentationFile.data;
  indexLeft = find(segAppVol==5);
  indexRight = find(segAppVol==6);
  [x, y, z] = ind2sub(size(segAppVol), indexLeft);
  coordsLeft = [x,y,z];
  [x, y, z] = ind2sub(size(segAppVol), indexRight);
  coordsRight = [x,y,z];
  [x, y, z] = ind2sub(size(roisLeftVolume), find(roisLeftVolume>0));
  coordsRoiLeft = [x,y,z];
  [x, y, z] = ind2sub(size(roisRightVolume), find(roisRightVolume>0));
  coordsRoiRight = [x,y,z];
  
  %addpath(genpath('/mnt/WholeBrainNewSeg-afni/matlab/vistasoft_nearpoints'));
  
  distOutLeft = nearpoints(coordsRoiLeft',coordsLeft');
  distOutRight = nearpoints(coordsRoiRight',coordsRight');
  leftVolOut = zeros(size(segAppVol)); 
  rightVolOut = zeros(size(segAppVol));
  
  leftVolOut(indexLeft(distOutLeft)) = roisLeftVolume(roisLeftVolume>0);
  rightVolOut(indexRight(distOutRight)) = roisRightVolume(roisRightVolume>0);
  
  segmentationFile.fname = sprintf('ClustersLeftLayer1_zzz.nii.gz') %#ok<NOPTS>
  segmentationFile.data = leftVolOut;
  writeFileNifti(segmentationFile);
  segmentationFile.fname = sprintf('ClustersRightLayer1_zzz.nii.gz') %#ok<NOPTS>
  segmentationFile.data = rightVolOut;
  writeFileNifti(segmentationFile);
  
end