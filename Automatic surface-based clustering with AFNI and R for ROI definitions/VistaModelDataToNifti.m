function VistaModelDataToNifti(VOLUME,T1name,expression,smoothLevel)

%Saves timing model data as niftis for clustering
%01/2020 BMH

load('allTmp.mat')
T1 = readFileNifti(T1name);
ves = rmGet(VOLUME{1}.rm.retinotopyModels{1}, 've');
%if task == 1 % timing
    x = rmGet(VOLUME{1}.rm.retinotopyModels{1}, 'x'); % duration
    y = rmGet(VOLUME{1}.rm.retinotopyModels{1}, 'y'); % period
    
    smoothing = smoothLevel;
    if smoothing > 0
        [~, ii] = intersectCols(VOLUME{1}.coords, VOLUME{1}.ROIs(end).coords);
        mask = zeros(size(VOLUME{1}.coords));
        mask(ii) = 1;
        [ves, conMat] = dhkGraySmooth(VOLUME{1},ves,[smoothing 0.6], [], mask);
        [x, conMat] = dhkGraySmooth(VOLUME{1},x,[smoothing 0.6], conMat, mask);
        [y, conMat] = dhkGraySmooth(VOLUME{1},y,[smoothing 0.6], conMat, mask); %#ok<ASGLU>
    end
    
    %Threshold data for a specific range of duration and period preferences
    %that demonstrate tuned selectivity.
    
    mask=eval(expression);
    x(~mask)=0;
    y(~mask)=0;
    ves(~mask)=0;
%     ves(x>1) = 0;
%     ves(y>1) = 0;
%     ves(x<0.05) = 0;
%     ves(y<0.05) = 0;
%     xtmp = x;
%     ytmp = y;
%     xtmp(x>1) = 0;
%     xtmp(y>1) = 0;
%     xtmp(x<0.05) = 0;
%     xtmp(y<0.05) = 0;
%     ytmp(x>1) = 0;
%     ytmp(y>1) = 0;
%     ytmp(x<0.05) = 0;
%     ytmp(y<0.05) = 0;
%     x = xtmp;
%     y = ytmp;
    %Determine which indices to set
    [~, ii] = intersectCols(VOLUME{1}.coords, VOLUME{1}.ROIs(1).coords);
    veTmp = allTmp;
    xTmp = allTmp;
    yTmp = allTmp;
    for n = 1:length(ii)
        veTmp(VOLUME{1}.coords(3,ii(n)), T1.dim(2)+1-VOLUME{1}.coords(2,ii(n)), T1.dim(1)+1-VOLUME{1}.coords(1,ii(n))) = ves(ii(n));
        xTmp(VOLUME{1}.coords(3,ii(n)), T1.dim(2)+1-VOLUME{1}.coords(2,ii(n)), T1.dim(1)+1-VOLUME{1}.coords(1,ii(n))) = x(ii(n));
        yTmp(VOLUME{1}.coords(3,ii(n)), T1.dim(2)+1-VOLUME{1}.coords(2,ii(n)), T1.dim(1)+1-VOLUME{1}.coords(1,ii(n))) = y(ii(n));
    end
    venii = T1;
    xii = T1;
    yii = T1;
    venii.data = veTmp;
    xii.data = xTmp;
    yii.data = yTmp;
    if smoothing > 0
        venii.fname = 'VarExpSmooth.nii.gz';
        xii.fname = 'PrefDurationSmooth.nii.gz';
        yii.fname = 'PrefPeriodSmooth.nii.gz';
    else
        venii.fname = 'VarExp.nii.gz';
        xii.fname = 'PrefDuration.nii.gz';
        yii.fname = 'PrefPeriod.nii.gz';
    end
    writeFileNifti(venii);
    writeFileNifti(xii);
    writeFileNifti(yii);
    
% elseif task == 2 % numerosity
%     x = rmGet(VOLUME{1}.rm.retinotopyModels{1}, 'x'); % numerosity
%     y = rmGet(VOLUME{1}.rm.retinotopyModels{1}, 'y'); % empty
%     
%     smoothing = smoothLevel;
%     if smoothing > 0
%         [~, ii] = intersectCols(VOLUME{1}.coords, VOLUME{1}.ROIs(end).coords);
%         mask = zeros(size(VOLUME{1}.coords));
%         mask(ii) = 1;
%         [ves, conMat] = dhkGraySmooth(VOLUME{1},ves,[smoothing 0.6], [], mask);
%         [x, conMat] = dhkGraySmooth(VOLUME{1},x,[smoothing 0.6], conMat, mask); %#ok<ASGLU>
%     end
%     
%     %Threshold data for a specific range of duration and period preferences
%     %that demonstrate tuned selectivity.
%     ves(x>log(6.99)) = 0;
%     ves(y>0) = 0;
%     ves(x<log(1.01)) = 0;
%     ves(y<0) = 0;
%     xtmp = x;
%     ytmp = y;
%     xtmp(x>log(6.99)) = 0;
%     xtmp(y>0) = 0;
%     xtmp(x<log(1.01)) = 0;
%     xtmp(y<0) = 0;
%     ytmp(x>log(6.99)) = 0;
%     ytmp(y>1) = 0;
%     ytmp(x<log(1.01)) = 0;
%     ytmp(y<0) = 0;
%     x = xtmp;
%     y = ytmp;
%     %Determine which indices to set
%     [~, ii] = intersectCols(VOLUME{1}.coords, VOLUME{1}.ROIs(1).coords);
%     veTmp = allTmp;
%     xTmp = allTmp;
%     yTmp = allTmp;
%     for n = 1:length(ii)
%         veTmp(VOLUME{1}.coords(3,ii(n)), T1.dim(2)+1-VOLUME{1}.coords(2,ii(n)), T1.dim(1)+1-VOLUME{1}.coords(1,ii(n))) = ves(ii(n));
%         xTmp(VOLUME{1}.coords(3,ii(n)), T1.dim(2)+1-VOLUME{1}.coords(2,ii(n)), T1.dim(1)+1-VOLUME{1}.coords(1,ii(n))) = x(ii(n));
%         yTmp(VOLUME{1}.coords(3,ii(n)), T1.dim(2)+1-VOLUME{1}.coords(2,ii(n)), T1.dim(1)+1-VOLUME{1}.coords(1,ii(n))) = y(ii(n));
%     end
%     venii = T1;
%     xii = T1;
%     yii = T1;
%     venii.data = veTmp;
%     xii.data = xTmp;
%     yii.data = yTmp;
%     if smoothing > 0
%         venii.fname = 'VarExpSmooth.nii.gz';
%         xii.fname = 'PrefNumerositySmooth.nii.gz';
%         yii.fname = 'PrefYSmooth.nii.gz';
%     else
%         venii.fname = 'VarExp.nii.gz';
%         xii.fname = 'PrefNumerosity.nii.gz';
%         yii.fname = 'PrefY.nii.gz';
%     end
%     writeFileNifti(venii);
%     writeFileNifti(xii);
%     writeFileNifti(yii);
% end
end