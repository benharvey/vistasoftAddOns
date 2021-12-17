for n=1:1000
tmp=paramsOrig.dotOrder(7:end);
tmp=tmp(randperm(length(tmp)));
params.dotOrder=[tmp(end-5:end) tmp];
savefile=['params_zperm_log_numbers_' num2str(n, '%04d') '.mat'];
save(savefile, 'params')
end

for n=1:1000;
for m=1:5;
    dataTYPES(dts(m)).retinotopyModelParams.paramsFile=['/media/Storage1/Ben/NumbersStudy/WholeBrain/Nell/NellCombined/Stimuli/params_zperm_log_numbers_' num2str(n, '%04d') '.mat'];
end
saveSession;
VOLUME{1}= rmClearSettings(VOLUME{1});
rmRunNumbersScriptLog(VOLUME{1},dts, 'BothNumROIs',1, [], hrfParams, ['nullPermutationsNoCycle-logNumber-', num2str(n, '%04d')], 5);
end

[tmp ROIindices]=intersectCols(VOLUME{1}.coords,  VOLUME{1}.ROIs(3).coords);

VE=zeros(1000, length(ROIindices));
for n=1:1000
    load(['nullPermutationsNoCycle-logNumber-', num2str(n, '%04d'), '-gFit.mat'])
    VE(n,:)=1-(model{1}.rss(ROIindices) ./ model{1}.rawrss(ROIindices));
end
VE(~isfinite(VE)) = 0;
VE = max(VE, 0);
VE = min(VE, 1);