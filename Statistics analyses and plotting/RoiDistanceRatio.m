function [data] = RoiDistanceRatio(view, targetRoiA, targetRoiB, sourceRoi, meansDT, modelROI)

if ~exist('meansDT','var') || isempty(meansDT)
    meansDT=viewGet(view, 'curdt');
end
 
if ~exist('targetRoiA','var') || isempty(targetRoiA)
else
    distancesA = RoiToRoiDist(targetRoiA, sourceRoi, view);
    distancesB = RoiToRoiDist(targetRoiB, sourceRoi, view);
    data.meanDist=mean([RoiToRoiDist(targetRoiA,targetRoiB, view) RoiToRoiDist(targetRoiB,targetRoiA, view)]);
    data.ratio=distancesA./(distancesB+distancesA);
end

[tmp, roiIndices]=intersectCols(view.coords, view.ROIs(sourceRoi).coords);
x0s=view.rm.retinotopyModels{1}.x0(roiIndices);
y0s=view.rm.retinotopyModels{1}.y0(roiIndices);
sigmas=view.rm.retinotopyModels{1}.sigma.major(roiIndices);
sigmaMinor=view.rm.retinotopyModels{1}.sigma.minor(roiIndices);
sigmaTheta=view.rm.retinotopyModels{1}.sigma.theta(roiIndices);
if isfield(view.rm.retinotopyModels{1}, 'exp')
    exp=view.rm.retinotopyModels{1}.exp(roiIndices);
end

ves=1-(view.rm.retinotopyModels{1}.rss(roiIndices)./view.rm.retinotopyModels{1}.rawrss(roiIndices));
ves(~isfinite(ves)) = 0;
ves = max(ves, 0);
ves = min(ves, 1);

if isfield(view.rm.retinotopyModels{1}, 'sigma2')
    sigmas2=view.rm.retinotopyModels{1}.sigma2.major(roiIndices);
    if size(view.rm.retinotopyModels{1}.x0,2)>size(view.rm.retinotopyModels{1}.beta,2)
        roiCoords = view.ROIs(modelROI).coords;
        allCrds  = viewGet(view,'coords');
        [tmp, roiCoords] = intersectCols(allCrds,roiCoords);
        sparseBetas=zeros([1 size(view.rm.retinotopyModels{1}.x0,2) size(view.rm.retinotopyModels{1}.beta,3)]);
        sparseBetas(1 ,roiCoords, :)=view.rm.retinotopyModels{1}.beta;
        betas=squeeze(sparseBetas(1, roiIndices,:));
    else
        betas=squeeze(view.rm.retinotopyModels{1}.beta(1, roiIndices,:));
    end
    data.sigmas2=sigmas2;
    data.betas=betas;
end

view=viewSet(view, 'curdt', meansDT);
% view=computeMeanMap(view,1,-1);
% meanSignal=view.map{1}(roiIndices);
% data.meanSignal=meanSignal;
data.meanSignal=zeros(size(x0s));


data.x0s=x0s;
data.y0s=y0s;
data.sigmas=sigmas;
data.sigmaMinor=sigmaMinor;
data.sigmaTheta=sigmaTheta;
data.ves=ves;

data.roiIndices=roiIndices';
if isfield(view.rm.retinotopyModels{1}, 'exp')
    data.exp=exp;
end

end

