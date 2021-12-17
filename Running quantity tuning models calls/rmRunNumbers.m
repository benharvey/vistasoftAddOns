function vw = rmRunNumbers(vw,roi,wSearch,models)

if ~exist('vw','var') || isempty(vw)
    vw = getCurView;
end
if ~exist('roi','var')
    roi = [];
end
if ~exist('wSearch','var') || isempty(wSearch)
    wSearch = 1; % Grid fit only: make sure we turn of coarse to fine
end
if ~exist('models','var') || isempty(models)
    models = {'1g'};%,'dog'};
end

% restrict roi to non NaN data or create one that covers all gray matter
%vw = roiRestrictToNonNan(vw);


% If you sample pRF size at 0.1
%'minrf',.25,'maxrf',14,'numbersigmas',139,...
% 0.2
%'minrf',.25,'maxrf',14,'numbersigmas',69,...

% If you sample pRF position at 0.1
% 'relativeGridStep',.4,...
% 'minrf',.25,...

%models = {'dog'}
for n=1:numel(models)
    switch lower(models{n})
        case '1g'
            % actual call with modified parameters - 1D 1 Gaussian model
            matFileName = sprintf('retModel-%s-1DGaussian-ZeroBlanks',datestr(now,'yyyymmdd-HHMMSS'));
            vw  = rmMain(vw,roi,wSearch,...
                'prf model',{'1dgaussian'},...
                'coarsetofine',false,...
                'coarseDecimate',0,... % alldata no blurring
                'minFieldSize',1,...
                'sampleRate',1,... % samplerate to make stimulus on
                'fieldsize',7,...
                'relativeGridStep',.8,... % sample x at relativeGridStep*minRF
                'minrf',.25,'maxrf',14,'numbersigmas',69,... % sample prfs at <0.1 (whatever the unit is)
                'spacesigmas','linear',...
                'matfilename',matFileName);
            
        case 'dog'
            % actual call with modified parameters - 1D DoG model
            matFileName = sprintf('retModel-%s-1DDoGaussian-ZeroBlanks',datestr(now,'yyyymmdd-HHMMSS'));
            sigmaRatio  = [1.1:.2:3 5 10];
            sigmaRatioMaxVal = 14.*sigmaRatio(1);
            
            vw  = rmMain(vw,roi,wSearch,...
                'prf model',{'1ddog'},...
                'coarsetofine',false,...
                'coarseDecimate',0,... % alldata no blurring
                'minFieldSize',1,...
                'sampleRate',1,... % samplerate to make stimulus on
                'fieldsize',7,...
                'relativeGridStep',.8,... % sample x at relativeGridStep*minRF
                'minrf',.25,'maxrf',14,'numbersigmas',69,... % sample prfs at <0.1 (whatever the unit is)
                'spacesigmas','linear',...
                'sigmaRatio',sigmaRatio,...
                'sigmaRatioMaxVal',sigmaRatioMaxVal,...
                'matfilename',matFileName);
    end
end
   
return
