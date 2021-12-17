function AFNIclusterMultiThreshold(paths, whichSubs, task, smoothLevel, selSizes, selVarExps, actualVE, mainDir)
%Surface-based clusterig at a range of variance explained thresholds
%This is pretty heavy, so runs on multiple cores for multiple subjects
%simultaneously

%01/2020, BMH & AF

% parpool(4)
% par
for nSubj = whichSubs
  
  subjDir = paths{nSubj}
  cd(subjDir)
  
  !rm *ttt.nii.gz;
  !rm -R surfaces_folder_left;
  !rm -R surfaces_folder_right;
  !rm -R modelResults_ttt_interp_surfaces_folder_left;
  !rm -R modelResults_ttt_interp_surfaces_folder_right;
  if smoothLevel > 0
      instr = '!3dcalc -a VarExpSmooth.nii.gz -expr "a*step(a)" -prefix VarExpSmooth_corr_ttt.nii.gz';
      commitEval(instr)
  else
      instr = '!3dcalc -a VarExp.nii.gz -expr "a*step(a)" -prefix VarExp_corr_ttt.nii.gz';
      commitEval(instr)
  end
  %combine modelling results
  if task == 1 % timing
      if smoothLevel > 0
          !3dTcat -prefix modelResults_ttt.nii.gz PrefDurationSmooth.nii.gz PrefPeriodSmooth.nii.gz VarExpSmooth_corr_ttt.nii.gz;
      else
          !3dTcat -prefix modelResults_ttt.nii.gz PrefDuration.nii.gz PrefPeriod.nii.gz VarExp_corr_ttt.nii.gz;
      end
  elseif task == 2 % numerosity
      if smoothLevel > 0
          !3dTcat -prefix modelResults_ttt.nii.gz PrefNumerositySmooth.nii.gz PrefYSmooth.nii.gz VarExpSmooth_corr_ttt.nii.gz;
      else
          !3dTcat -prefix modelResults_ttt.nii.gz PrefNumerosity.nii.gz PrefY.nii.gz VarExp_corr_ttt.nii.gz;
      end
  end
  %modeling in orig space
  !3drefit -space ORIG -view orig modelResults_ttt.nii.gz; 
  %anat in orig space
  !3dcopy t1_1mm.nii.gz t1_1mm_copy_ttt.nii.gz; 
  !3drefit -space ORIG -view orig t1_1mm_copy_ttt.nii.gz;
  %segmentation in orig space
  !3dcopy SegmentationLayer1.nii.gz SegmentationLayer1_ttt.nii.gz;
  !3drefit -space ORIG -view orig SegmentationLayer1_ttt.nii.gz;
  %Left segmentation
  instr = '!3dcalc -a SegmentationLayer1_ttt.nii.gz -expr "within(a,2.9,3.1)" -prefix leftSeg_ttt.nii.gz';
  commitEval(instr)
  %Right segmentation
  instr = '!3dcalc -a SegmentationLayer1_ttt.nii.gz -expr "within(a,3.9,4.1)" -prefix rightSeg_ttt.nii.gz';
  commitEval(instr)
  
  %generate surfaces
  instr = sprintf('!Rscript %s/generateSurfacesFromBoundaries.R leftSeg_ttt.nii.gz 10 0 800 1', mainDir);
  commitEval(instr)
  !mv surfaces_folder/ surfaces_folder_left;
  instr = sprintf('!Rscript %s/generateSurfacesFromBoundaries.R rightSeg_ttt.nii.gz 10 0 800 1', mainDir);
  commitEval(instr)
  !mv surfaces_folder/ surfaces_folder_right/;
  
  %copy maps on surface_left
  instr = sprintf('!Rscript %s/interpolateSurface.R modelResults_ttt.nii.gz surfaces_folder_left/', mainDir); 
  commitEval(instr)
  %copy maps on surface_right
  instr = sprintf('!Rscript %s/interpolateSurface.R modelResults_ttt.nii.gz surfaces_folder_right/', mainDir);
  commitEval(instr)
  
  % load model data for each hemisphere and get varExp distribution (quantiles)
  !cp modelResults_ttt_interp_surfaces_folder_left/boundary00_sm_modelResults_ttt_surf.1D.dset modelResults_ttt_interp_surfaces_folder_left/boundary00_sm_modelResults_ttt_surf.1D.dset.txt
  modelLeft = readtable('modelResults_ttt_interp_surfaces_folder_left/boundary00_sm_modelResults_ttt_surf.1D.dset.txt');
  modelLeft = table2array(modelLeft);
  !cp modelResults_ttt_interp_surfaces_folder_right/boundary00_sm_modelResults_ttt_surf.1D.dset modelResults_ttt_interp_surfaces_folder_right/boundary00_sm_modelResults_ttt_surf.1D.dset.txt
  modelRight = readtable('modelResults_ttt_interp_surfaces_folder_right/boundary00_sm_modelResults_ttt_surf.1D.dset.txt');
  modelRight = table2array(modelRight);
  
  varExp_nodeLeft = modelLeft((find(modelLeft(:,9)>0.001)), 9);
  varExp_nodeRight = modelRight((find(modelRight(:,9)>0.001)) ,9);
%   #sizeThrLoop <- seq( 15, 85, 15 )
%   #varQuantThrLoopLeft <- quantile( varExp_nodeLeft, c( seq( 0.40, 0.95, 0.05 ), seq( 0.96, 0.998, 0.002 ) ) )
%   #varQuantThrLoopRight <- quantile( varExp_nodeRight,  c( seq( 0.40, 0.95, 0.05 ), seq( 0.96, 0.998, 0.002 ) )  )
  sizeThrLoop = selSizes;
  
  %if using percentile
  if actualVE == 0
      varQuantThrLoopLeft = quantile(varExp_nodeLeft, selVarExps);
      varQuantThrLoopRight = quantile(varExp_nodeRight,  selVarExps);
  else
      varQuantThrLoopLeft = selVarExps;
      varQuantThrLoopRight = selVarExps;
  end
  
  outputClusteringLeftQuant = zeros(length(sizeThrLoop), length(varQuantThrLoopLeft))+999;
  outputClusteringRightQuant = zeros(length(sizeThrLoop), length(varQuantThrLoopRight))+999;
  for sizeCount = 1:length(sizeThrLoop)
    for varCount = 1:length(varQuantThrLoopLeft)
      
      thisSelSize = sizeThrLoop(sizeCount);
      selVarLeft = varQuantThrLoopLeft(varCount);
      selVarRight = varQuantThrLoopRight(varCount);
      
      disp( sprintf('size: %s, varExp (left): %s, varExp (right): %s', num2str(thisSelSize ), num2str(round(selVarLeft, 3 )), num2str( round( selVarRight, 3 ) ) ) ) %#ok<DSPS>
      % cluster the surfaces based on variance explained
      instr = sprintf( '!SurfClust -spec surfaces_folder_left/spec.surfaces.smoothed -surf_A surfaces_folder_left/boundary00_sm.1D.coord  -input modelResults_ttt_interp_surfaces_folder_left/boundary00_sm_modelResults_ttt_surf.1D.dset 8 -rmm -1 -amm2 %s -thresh_col 8 -thresh %s > leftTableLoop_ttt.1D', num2str(thisSelSize), num2str(selVarLeft) );
      commitEval(instr);
      instr = sprintf( '!SurfClust -spec surfaces_folder_right/spec.surfaces.smoothed -surf_A surfaces_folder_right/boundary00_sm.1D.coord  -input modelResults_ttt_interp_surfaces_folder_right/boundary00_sm_modelResults_ttt_surf.1D.dset 8 -rmm -1 -amm2 %s -thresh_col 8 -thresh %s > rightTableLoop_ttt.1D', num2str(thisSelSize), num2str(selVarRight) );
      commitEval(instr);
      
      !cp leftTableLoop_ttt.1D leftTableLoop_ttt.1D.txt
      !cp rightTableLoop_ttt.1D rightTableLoop_ttt.1D.txt
      
      fid = fopen('leftTableLoop_ttt.1D.txt');
      leftTableLoop = textscan(fid, repmat('%f', 1, 23), 'HeaderLines', 26);
      fid = fopen('rightTableLoop_ttt.1D.txt');
      rightTableLoop = textscan(fid, repmat('%f', 1, 23), 'HeaderLines', 26);
      %leftTableLoop = leftTableLoop(26:end, 1);
      %leftTableLoop = split(table2array(leftTableLoop));
%       rightTableLoop = readtable('rightTableLoop_ttt.1D.txt');
%       rightTableLoop = rightTableLoop(26:end, 1);
      %rightTableLoop = split(table2array(rightTableLoop));
      
        outputClusteringLeftQuant(sizeCount, varCount) = size(leftTableLoop{1},1);
      
        outputClusteringRightQuant(sizeCount, varCount) = size(rightTableLoop{1},1);
    end
  end
  cd(mainDir)
  filename = sprintf('dataClustering_Date_%s.mat', datestr(now));
  commitSave(filename, outputClusteringLeftQuant, outputClusteringRightQuant);
end
end