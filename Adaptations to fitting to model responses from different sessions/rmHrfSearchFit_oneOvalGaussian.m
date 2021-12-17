function prediction = rmHrfSearchFit_oneOvalGaussian(model, params, loopSlices, wProcess,varexp,allstimimages)
% rmHrfSearchFit_oneGaussian - make predictions without HRF
%
% prediction = rmHrfSearchFit_oneGaussian(model, params, loopSlices, wProcess)
%
% 10/2009: SD & WZ wrote it.


% intiate data stucture for our pRF values
tmp.x    = [];
tmp.y    = [];
tmp.s    = [];
tmp.s_major    = [];
tmp.s_minor    = [];
tmp.s_theta    = [];
tmp.ve   = [];


% Get the data required to build our pRF
for slice=loopSlices,
    % now we extract only the data from that slice
    s = rmSliceGet(model,slice);
    
    % store
    tmp.x = [tmp.x s{1}.x0(wProcess)];
    tmp.y = [tmp.y s{1}.y0(wProcess)];
    tmp.s = [tmp.s s{1}.s(wProcess)];
    tmp.s_major = [tmp.s_major s{1}.s_major(wProcess)];
    tmp.s_minor = [tmp.s_minor s{1}.s_minor(wProcess)];
    tmp.s_theta = [tmp.s_theta s{1}.s_theta(wProcess)];
    tmp.ve = [tmp.ve varexp(wProcess)];
end


if isfield(s{1}, 'exp')
    s{1}.exponent=s{1}.exp;
end

if isfield(s{1}, 'exponent')
    tmp.exponent =[];
    for slice=loopSlices,
        s = rmSliceGet(model,slice);
        tmp.exponent = [tmp.exponent s{1}.exponent(wProcess)];
    end
    periods=[0.05:0.05:1 2.1 0.1 0.6 0.9 1 0.15 0.65 0.85 1 0.2 0.7 0.8 1 0.25 0.75 1 0.3 0.7 0.8 1 0.35 0.65 0.85 1 0.4 0.6 0.9 1 0.45 0.55 0.95 1 0.5 1 0.55 1 0.6 1 0.65 1 0.7 1 0.75 1 0.8 1 0.85 1 0.9  1 0.95 1 1 2.1];
    freq=5./periods;
end

% build pRF and make our predictions
n = numel(tmp.x);
s = [[1:ceil(n./1000):n-2] n+1]; %#ok<NBRAK>
prediction = zeros(size(allstimimages,1),n);
fprintf(1,'[%s]:Making predictions without HRF for each voxel (%d):',mfilename,n);
drawnow;tic;
for n=1:numel(s)-1,
    % make rfs
    rf   = rfGaussian2d(params.analysis.X, params.analysis.Y,...
        tmp.s_major(s(n):s(n+1)-1), ...
        tmp.s_minor(s(n):s(n+1)-1),...
        tmp.s_theta(s(n):s(n+1)-1), ...
        tmp.x(s(n):s(n+1)-1), ...
        tmp.y(s(n):s(n+1)-1));
        
    if isfield(tmp, 'exponent')
        count=1;
        for sn=s(n):(s(n+1)-1)
            scale=1./(freq'.^tmp.exponent(sn));
            scale=scale.*freq';
            rf(:,count)=rf(:,count)./scale;
            count=count+1;
        end
    end
    
    % convolve with stimulus
    pred = allstimimages*rf;
    
    % store
    prediction(:,s(n):s(n+1)-1) = pred;
    fprintf(1,'.');drawnow;
end;
clear n s rf pred;
fprintf(1, 'Done[%d min].\t(%s)\n', round(toc/60), datestr(now));
drawnow;
