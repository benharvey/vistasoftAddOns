params = viewGet(VOLUME{1}, 'rmParams');
params = rmRecomputeParams(VOLUME{1}, params);
roi = VOLUME{1}.selectedROI;
roi = tc_roiStruct(VOLUME{1}, roi);

M = rmPlotGUI_getModel(VOLUME{1}, roi, 0);
M.pred=zeros(size(M.tSeries));
for ii=1:length(M.coords)
    coords = M.coords(ii);
    M.ui.voxel.sliderHandle=ii;
    [pred RF rfParams variance_explained blanks] = rmPlotGUI_makePrediction(M, coords,ii, params.seperateRunBetas);
    M.pred(:,ii)=pred;
    for i=1:4
        M.rawrss(i,ii)=sum(M.tSeries(((i-1)*44+1):(i*44),ii).^2);
        M.rss(i,ii)=sum((M.tSeries(((i-1)*44+1):(i*44),ii)-M.pred((((i-1)*44)+1):(i*44), ii)).^2);
        M.ve(i,ii)=1 - (M.rss(i,ii) ./ M.rawrss(i,ii));
    end
    M.rawrssAll(ii)=sum(M.tSeries(:,ii).^2);
    M.rssAll(ii)=sum((M.tSeries(:,ii)-M.pred(:, ii)).^2);
    M.veAll(ii)=1 - (M.rssAll(ii) ./ M.rawrssAll(ii));
end
M.ve(~isfinite(M.ve)) = 0;
M.ve = max(M.ve, 0);
M.ve = min(M.ve, 1);
M.veAll(~isfinite(M.veAll)) = 0;
M.veAll = max(M.veAll, 0);
M.veAll = min(M.veAll, 1);


% 
% 
% rss=(groupParams.TotalCirc(1).allConditionsRSSlin.NumRightIPS1 + groupParams.TotalCirc(2).allConditionsRSSlin.NumRightIPS1 + groupParams.TotalCirc(3).allConditionsRSSlin.NumRightIPS1 + groupParams.TotalCirc(4).allConditionsRSSlin.NumRightIPS1);
% rawrss=(groupParams.TotalCirc(1).allConditionsRawRSSlin.NumRightIPS1 + groupParams.TotalCirc(2).allConditionsRawRSSlin.NumRightIPS1 + groupParams.TotalCirc(3).allConditionsRawRSSlin.NumRightIPS1 + groupParams.TotalCirc(4).allConditionsRawRSSlin.NumRightIPS1);
% 
% rss=(groupParams.TotalCirc(3).allConditionsRSSlin.NumRightIPS1);
% rawrss=(groupParams.TotalCirc(3).allConditionsRawRSSlin.NumRightIPS1);
% 
% 
% ves=1-(rss./rawrss);
% mean(ves)
% groupParams
% groupParams.TotalCirc
% ves(ves>0)=0
% ves=1-(rss./rawrss);
% ves(ves<0)=0;
% ves(ves>1)=1;
% mean(ves)