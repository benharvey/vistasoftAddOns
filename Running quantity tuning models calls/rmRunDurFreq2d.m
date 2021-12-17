function vw = rmRunDurFreq2d(vw,dataTypes, roi,wSearch,models, hrfParams, matFileName, maxCores, separateBetas)

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
    models = {'1g'};
end

if ~exist('hrfParams','var') || isempty(hrfParams)
    hrfParams = {'two gammas (SPM style)', [5.4000 5.2000 10.8000 7.3500 0.3500]};
end

if ~exist('matFileName','var') || isempty(matFileName)
    matFileName = sprintf('retModel-%s-2dDGaussian-lin-FreqInvDur',datestr(now,'yyyymmdd-HHMMSS')); 
end
if ~exist('maxCores','var') || isempty(maxCores)
    maxCores=5; % Grid fit only: make sure we turn of coarse to fine
end

if ~exist('separateBetas','var') || isempty(separateBetas)
    separateBetas=1;
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

%For Utrecht 7T
highestValue=20;
sampleRate=highestValue*(0.025);
gridStep=0.1/sampleRate; %Maybe smaller, also numberSigmas
blankValue=20;
nSigmas=20;

if wSearch==12
    matFileName = sprintf('retModel-%s-Lin-1dGaussianXcompressiveY-DurFreq-%s',datestr(now,'yyyymmdd-HHMMSS'), num2str(highestValue,2)); 
    modelType='one gaussian x compressive y';
    wSearch=4;
elseif wSearch==13
    matFileName = sprintf('retModel-%s-Lin-compressiveXcompressiveYNoNormOccupancy-DurFreq-%s',datestr(now,'yyyymmdd-HHMMSS'), num2str(highestValue,2)); 
    modelType='compressive x plus compressive y';
    wSearch=1;
elseif wSearch==14
    matFileName = sprintf('retModel-%s-Lin-1dGaussianXlinearY-DurFreq-%s',datestr(now,'yyyymmdd-HHMMSS'), num2str(highestValue,2)); 
    modelType='one gaussian x linear y';
    %To avoid search fit, particularly avoid searching compressive exponent y
    wSearch=1;    
    sampleRate=highestValue*(0.0025);
    gridStep=0.01/sampleRate;
    nSigmas=100;
elseif wSearch==15 %Zhou model
    matFileName = sprintf('retModel-%s-Lin-linearXcompressiveYNoNormOccupancy-DurFreq-%s',datestr(now,'yyyymmdd-HHMMSS'), num2str(highestValue,2)); 
    modelType='linear x plus compressive y';
    wSearch=1;
elseif wSearch==16 %Stigliani model
    matFileName = sprintf('retModel-%s-Lin-linearXlinearYNoNormOccupancy-DurFreq-%s',datestr(now,'yyyymmdd-HHMMSS'), num2str(highestValue,2)); 
    modelType='linear x plus linear y';
    wSearch=1;
else
    matFileName = sprintf('retModel-%s-Lin-2dOvalGaussian-DurFreq-%s',datestr(now,'yyyymmdd-HHMMSS'), num2str(highestValue,2)); 
    modelType='one gaussian';
end

%models = {'dog'}
ncores=length(dataTypes);
if ncores>maxCores
    ncores=maxCores;
end
% if ncores>1
%     parpool(ncores)
% end
%par
for dt=1:length(dataTypes)
    for n=1:numel(models)
        switch lower(models{n})
            case '1g'
                % actual call with modified parameters - 1D 1 Gaussian model

                tmpv=rmMain([1 dataTypes(dt)],roi,wSearch,...
                    'prf model',{modelType},...
                    'coarsetofine',false,...
                    'coarseDecimate',0,... % alldata no blurring
                    'minFieldSize',0,...
                    'sampleRate',sampleRate,... % samplerate to make stimulus on
                    'fieldsize',highestValue,...
                    'relativeGridStep',gridStep,... % sample x at relativeGridStep*minRF
                    'minrf',sampleRate,'maxrf',blankValue,'numbersigmas',nSigmas,... % sample prfs at <0.1 (whatever the unit is)
                    'spacesigmas','linear',...
                    'hrffitthreshve', 0.2,...
                    'hrf', hrfParams,...
                    'matfilename',matFileName,...
                    'separate betas', separateBetas,...
                    'stimx', [0.05:0.05:1 2],...
                    'stimy', 1./[0.05:0.05:1 2.1]);
        end
%             hrftmp=viewGet(tmpv, 'rmhrf')
%             hrfOut{dt}=hrftmp{2};
    end
end
delete(gcp('nocreate'))

return