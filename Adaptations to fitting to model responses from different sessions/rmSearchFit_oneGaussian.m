function model = rmSearchFit_oneGaussian(model, data, params, wProcess, t)
% rmSearchFit_oneGaussian - wrapper for 'fine' one Gaussian fit
%
% model = rmSearchFit_oneGaussian(model, data, params, wProcess, t);
%
% 2008/01 SOD: split of from rmSearchFit.
% 2010/02 SOD: cleanup.

% fminsearch options
searchOptions = params.analysis.fmins.options;
expandRange   = params.analysis.fmins.expandRange;
 
% convert to double just in case
params.analysis.X = double(params.analysis.X);
params.analysis.Y = double(params.analysis.Y);
params.analysis.allstimimages = double(params.analysis.allstimimages);
data = double(data);

% get starting upper and lower range and reset TolFun 
% (raw rss computation (similar to norm) and TolFun adjustments)
[range TolFun] = rmSearchFit_range(params,model,data);

% amount of negative fits
nNegFit  = 0;
vethresh = params.analysis.fmins.vethresh;
trends   = t.trends;
if strcmp(model.desc, 'compressive monotonic x plus compressive monotonic y(xcomp,ycomp, positive only)') || strcmp(model.desc, '1D pRF fit plus compressive monotonic Y(x,ycomp,sigma, positive only)')
    if params.seperateRunBetas
        t_id   = t.dcid+(2*length(params.stim));
    else
        t_id   = t.dcid+2;
    end
elseif length(params.stim)>1 && params.seperateRunBetas;
    t_id           = t.dcid+length(params.stim);
else
    t_id           = t.dcid+1;
end

%-----------------------------------
% Go for each voxel
%-----------------------------------
progress = 0;tic;
for ii = 1:numel(wProcess),
    % progress monitor (10 dots)
    if floor(ii./numel(wProcess)*10)>progress,
        % print out estimated time left
        if progress==0,
            esttime = toc.*10;
            if floor(esttime./3600)>0,
                fprintf(1,'[%s]:Estimated processing time: %d voxels: %d hours.\n',...
                    mfilename,numel(wProcess),ceil(esttime./3600));
            else
                fprintf(1,'[%s]:Estimated processing time: %d voxels: %d minutes.\n',...
                    mfilename,numel(wProcess),ceil(esttime./60));
            end;
            fprintf(1,'[%s]:Nonlinear optimization (x,y,sigma):',mfilename);
        end;
        fprintf(1,'.');drawnow;
        progress = progress + 1;
    end;

    % volume index
    vi = wProcess(ii);
    vData = data(:,ii);
    
    % reset tolFun: Precision of evaluation function. 
    % We define RMS improvement relative to the initial raw 'no-fit' data
    % RMS. So, 1 means stop if there is less than 1% improvement on the fit:
    % searchOptions = optimset(searchOptions,'tolFun',optimget(params.analysis.fmins.options,'tolFun')./100.*rawrss);
    % optimset and optimget are a little slow so:
    searchOptions.TolFun = TolFun(ii);
    
    % actual fitting routine
    if searchOptions.MaxIter>0
        try
            if isfield(params, 'seperateRunBetas') && params.seperateRunBetas && isfield(params.analysis, 'exp')
                outParams = ...
                    fmincon(@(x) rmModelSearchFit_oneGaussian_exp_seperateBetas(x,vData,...
                    params.analysis.X,...
                    params.analysis.Y,...
                    params.analysis.allstimimages,trends,...
                    length(params.stim)),...
                    range.start(:,vi),[],[],[],[],range.lower(:,vi),range.upper(:,vi),...
                    @(x) distanceCon(x,range.start(:,vi),range.step(:,vi).*expandRange),searchOptions);
            elseif isfield(params, 'seperateRunBetas') && params.seperateRunBetas && strcmp(model.desc, '1D pRF fit with compressive monotonic y(x,ycomp,sigma, positive only)')
                outParams = ...
                    fmincon(@(x) rmModelSearchFit_xGaussianYCompressive_seperateBetas(x,vData,...
                    params.analysis.X,...
                    params.analysis.Y,...
                    params.analysis.allstimimages,trends,...
                    length(params.stim)),...
                    range.start(:,vi),[],[],[],[],range.lower(:,vi),range.upper(:,vi),...
                    @(x) distanceCon(x,range.start(:,vi),range.step(:,vi).*expandRange),searchOptions);    
            elseif isfield(params, 'seperateRunBetas') && params.seperateRunBetas
                outParams = ...
                    fmincon(@(x) rmModelSearchFit_oneGaussian_seperateBetas(x,vData,...
                    params.analysis.X,...
                    params.analysis.Y,...
                    params.analysis.allstimimages,trends,...
                    length(params.stim)),...
                    range.start(:,vi),[],[],[],[],range.lower(:,vi),range.upper(:,vi),...
                    @(x) distanceCon(x,range.start(:,vi),range.step(:,vi).*expandRange),searchOptions);
            elseif isfield(params.analysis, 'exp')
                outParams = ...
                    fmincon(@(x) rmModelSearchFit_oneGaussian_exp(x,vData,...
                    params.analysis.X,...
                    params.analysis.Y,...
                    params.analysis.allstimimages,trends),...
                    range.start(:,vi),[],[],[],[],range.lower(:,vi),range.upper(:,vi),...
                    @(x) distanceCon(x,range.start(:,vi),range.step(:,vi).*expandRange),searchOptions);
            elseif strcmp(model.desc, '1D pRF fit with compressive monotonic y(x,ycomp,sigma, positive only)')
                outParams = ...
                    fmincon(@(x) rmModelSearchFit_xGaussianYCompressive(x,vData,...
                    params.analysis.X,...
                    params.analysis.Y,...
                    params.analysis.allstimimages,trends),...
                    range.start(:,vi),[],[],[],[],range.lower(:,vi),range.upper(:,vi),...
                    @(x) distanceCon(x,range.start(:,vi),range.step(:,vi).*expandRange),searchOptions);
            else
                outParams = ...
                    fmincon(@(x) rmModelSearchFit_oneGaussian(x,vData,...
                    params.analysis.X,...
                    params.analysis.Y,...
                    params.analysis.allstimimages,trends),...
                    range.start(:,vi),[],[],[],[],range.lower(:,vi),range.upper(:,vi),...
                    @(x) distanceCon(x,range.start(:,vi),range.step(:,vi).*expandRange),searchOptions);
            end
        catch
            outParams = range.start(:,vi);
        end
    else
        outParams = range.start(:,vi);
    end

    %[ outParams range.lower(:,vi) range.start(:,vi) range.upper(:,vi)]

    % make RF, prediction and get rss,b
    if strcmp(model.desc, '1D pRF fit with compressive monotonic y(x,ycomp,sigma, positive only)')
        Xv = params.analysis.X-outParams(1);
        rf = exp( -.5 *((Xv.*Xv ./ outParams(3)).^2)).*((params.analysis.Y.^outParams(2))./params.analysis.Y);
        %For flat prediction
        %rf=ones(size(rf));
        X  = [params.analysis.allstimimages * rf trends];
        if isfield(params, 'seperateRunBetas') && params.seperateRunBetas
            framesPerStim=length(vData)/length(params.stim);
            b=zeros(size(X,2), length(params.stim));
            for run=1:length(params.stim)
                indices=((run-1)*framesPerStim+1):(run*framesPerStim);
                b(:,run) = pinv(X(indices,:))*vData(indices);
                residual(indices)=vData(indices)-X(indices,:)*b(:,run);
            end
            btmp=b(1,:);
            for run=1:length(params.stim)
                btmp=[btmp b((1+run),run)];
            end
            b=btmp;
            rss=norm(residual).^2;
        else
            b    = pinv(X)*vData;
            rss  = norm(vData-X*b).^2;
        end
%     elseif strcmp(model.desc, '1D pRF fit with compressive monotonic y(x,ycomp,sigma, positive only)')
%         Xv = params.analysis.X-outParams(1);
%         rf = exp( -.5 *((Xv.*Xv ./ outParams(3)).^2)).*((params.analysis.Y.^outParams(2))./params.analysis.Y);
%         X  = [params.analysis.allstimimages * rf trends];
%         b    = pinv(X)*vData;
%         rss  = norm(vData-X*b).^2;
        
    elseif strcmp(model.desc, 'compressive monotonic x plus compressive monotonic y(xcomp,ycomp, positive only)')
        if searchOptions.MaxIter>0
            rf=params.analysis.X.^outParams(1);
            prediction = params.analysis.allstimimages*rf;
            rf2   = (params.analysis.Y.^outParams(2))./params.analysis.Y;
            prediction2 = params.analysis.allstimimages*rf2;
            X  = [prediction prediction2 trends];
            if isfield(params, 'seperateRunBetas') && params.seperateRunBetas
                framesPerStim=length(vData)/length(params.stim);
                b=zeros(size(X,2), length(params.stim));
                for run=1:length(params.stim)
                    indices=((run-1)*framesPerStim+1):(run*framesPerStim);
                    b(:,run) = pinv(X(indices,:))*vData(indices);
                    residual(indices)=vData(indices)-X(indices,:)*b(:,run);
                end
                btmp=b(1,:);
                btmp=[btmp b(2,:)];
                for run=1:length(params.stim)
                    btmp=[btmp b((2+run),run)];
                end
                b=btmp;
                rss=norm(residual).^2;
            else
                b    = pinv(X)*vData;
                rss  = norm(vData-X*b).^2;
            end
            if b(2)<0
                b(1)=b(2);
            end
        else
            rf=params.analysis.X.^outParams(1);
            prediction1 = params.analysis.allstimimages*rf;
            rf2   = (params.analysis.Y.^outParams(2))./params.analysis.Y;
            prediction2 = params.analysis.allstimimages*rf2;
            
           

            if isfield(params, 'seperateRunBetas') && params.seperateRunBetas
                framesPerStim=length(vData)/length(params.stim);
                for run=1:length(params.stim)
                    %Because the model has two componenets with independnt
                    %amplitudes, fix the ratio between these when evaluation the
                    %model fit. Otherwise, we provide an extra degree of freedom to
                    %this model.
                    indices=((run-1)*framesPerStim+1):(run*framesPerStim);
                    ratio=model.b(run, vi)/model.b(run+length(params.stim),vi);
                    prediction=prediction1(indices).*model.b(run, vi)+prediction2(indices).*model.b(run+length(params.stim),vi);
                    X  = [prediction trends(indices,:)];  
                    b = pinv(X)*vData(indices);
                    
                    if b(1)<0
                        b(1)=0;
                        b(2:end)= pinv(X(:,2:end))*vData(indices);
                    end
                    rsstotal(run)  = norm(vData(indices)-X*b).^2;
                    
                    %To account for ratio of betas, add an extra line with
                    %first beta / ratio
                    b2=zeros(2, size(b,2));
                    if model.b(run,vi)>0 && model.b(run+length(params.stim),vi)>0
                        b2(1)=b(1);
                        b2(2)=b(1)./ratio;
                    elseif model.b(run,vi)==0
                        b2(1)=zeros(1, size(b,2));
                        b2(2)=b(1);
                    elseif model.b(run+length(params.stim),vi)==0
                        b2(1)=b(1);
                        b2(2)=zeros(1, size(b,2));
                    end
                    bout([run run+length(params.stim) run+length(params.stim)*2])=[b2; b(run+1)];
                end
                b=bout;
                rss=sum(rsstotal);

            else
                ratio=model.b(1, vi)/model.b(2,vi);
                prediction=prediction1.*model.b(1, vi)+prediction2.*model.b(2,vi);
                X  = [prediction trends];
                b    = pinv(X)*vData;
                if b(1)<0
                    b(1)=0;
                    b(2:end)= pinv(X(:,2:end))*vData;
                end
                rss  = norm(vData-X*b).^2;
                
                %To account for ratio of betas, add an extra line with
                %first beta / ratio
                b2=zeros(2, size(b,2));
                if model.b(1,vi)>0 && model.b(2,vi)>0
                    b2(1,:)=b(1,:);
                    b2(2,:)=b(1,:)./ratio;
                elseif model.b(1,vi)==0
                    b2(1,:)=zeros(1, size(b,2));
                    b2(2,:)=b(1,:);
                elseif model.b(2,vi)==0
                    b2(1,:)=b(1,:);
                    b2(2,:)=zeros(1, size(b,2));
                end
                b=[b2; b(2:end)]; 
            end
        end
    elseif strcmp(model.desc, '1D pRF fit plus compressive monotonic Y(x,ycomp,sigma, positive only)')
        Xv = params.analysis.X-outParams(1);
        Yv = params.analysis.Y.*0;
        rf = exp( (Yv.*Yv + Xv.*Xv) ./ (-2.*(outParams(3).^2)) );
        prediction = params.analysis.allstimimages*rf;
        rf2   = (params.analysis.Y.^outParams(2))./params.analysis.Y;
        prediction2 = params.analysis.allstimimages*rf2;
        X  = [prediction prediction2 trends];
        if isfield(params, 'seperateRunBetas') && params.seperateRunBetas
            framesPerStim=length(vData)/length(params.stim);
            b=zeros(size(X,2), length(params.stim));
            for run=1:length(params.stim)
                indices=((run-1)*framesPerStim+1):(run*framesPerStim);
                b(:,run) = pinv(X(indices,:))*vData(indices);
                residual(indices)=vData(indices)-X(indices,:)*b(:,run);
            end
            btmp=b(1,:);
            btmp=[btmp b(2,:)];
            for run=1:length(params.stim)
                btmp=[btmp b((2+run),run)];
            end
            b=btmp;
            rss=norm(residual).^2;
        else
            b    = pinv(X)*vData;
            rss  = norm(vData-X*b).^2;
        end
    else
        Xv = params.analysis.X-outParams(1);
        Yv = params.analysis.Y-outParams(2);
        rf = exp( (Yv.*Yv + Xv.*Xv) ./ (-2.*(outParams(3).^2)) );
        X  = [params.analysis.allstimimages * rf trends];
        if isfield(params, 'seperateRunBetas') && params.seperateRunBetas
            framesPerStim=length(vData)/length(params.stim);
            b=zeros(size(X,2), length(params.stim));
            for run=1:length(params.stim)
                indices=((run-1)*framesPerStim+1):(run*framesPerStim);
                b(:,run) = pinv(X(indices,:))*vData(indices);
                residual(indices)=vData(indices)-X(indices,:)*b(:,run);
            end
            btmp=b(1,:);
            for run=1:length(params.stim)
                btmp=[btmp b((1+run),run)];
            end
            b=btmp;
            rss=norm(residual).^2;
        else
            b    = pinv(X)*vData;
            rss  = norm(vData-X*b).^2;
        end
    end
    if isfield(params.analysis, 'exp')
        periods=[0.05:0.05:1 2.1 0.1 0.6 0.9 1 0.15 0.65 0.85 1 0.2 0.7 0.8 1 0.25 0.75 1 0.3 0.7 0.8 1 0.35 0.65 0.85 1 0.4 0.6 0.9 1 0.45 0.55 0.95 1 0.5 1 0.55 1 0.6 1 0.65 1 0.7 1 0.75 1 0.8 1 0.85 1 0.9  1 0.95 1 1 2.1];
        freq=5./periods;
        scale=1./(freq'.^outParams(end));
        scale=scale.*freq';
        rf=rf./scale;
        X  = [params.analysis.allstimimages * rf trends];
        b    = pinv(X)*vData;
        rss  = norm(vData-X*b).^2;
    end

    % store results only if the first beta is positive, somehow fmincon
    % outputs negative fits. If the fit is negative keep old (grid) fit. We
    % do adjust the rss, so it won't be accidentally counted as a 'good'
    % fit. 
    if b(1)>=0,
        model.x0(vi)   = outParams(1);
        model.y0(vi)   = outParams(2);
        model.s(vi)    = outParams(3);
        model.s_major(vi)    = outParams(3);
        model.s_minor(vi)    = outParams(3);
        model.s_theta(vi)    = 0;
        model.rss(vi)  = rss;
        
        if isfield(params, 'seperateRunBetas') && params.seperateRunBetas && (strcmp(model.desc, 'compressive monotonic x plus compressive monotonic y(xcomp,ycomp, positive only)') || strcmp(model.desc, '1D pRF fit plus compressive monotonic Y(x,ycomp,sigma, positive only)'))
            model.b([1:2*length(params.stim) t_id], vi)=b;
        elseif strcmp(model.desc, 'compressive monotonic x plus compressive monotonic y(xcomp,ycomp, positive only)') || strcmp(model.desc, '1D pRF fit plus compressive monotonic Y(x,ycomp,sigma, positive only)')
            model.b([1 2 t_id],vi)  = b;
        elseif isfield(params, 'seperateRunBetas') && params.seperateRunBetas
            model.b([1:length(params.stim) t_id], vi)=b;
        else
            model.b([1 t_id],vi)  = b;
        end
        if isfield(params.analysis, 'exp')
            model.exp(vi) = outParams(end);
        end
    else
        % change the percent variance explained to be just under the
        % current vethresh. So it counts as a 'coarse'-fit but can still be
        % included in later 'fine'-fits
        model.rss(vi)  = (1-max((vethresh-0.01),0)).*model.rawrss(vi);
        nNegFit = nNegFit + 1;
    end
end

% end time monitor
et  = toc;
if floor(et/3600)>0,
    fprintf(1,'Done [%d hours].\n',ceil(et/3600));
else
    fprintf(1,'Done [%d minutes].\n',ceil(et/60));
end;
fprintf(1,'[%s]:Removed negative fits: %d (%.1f%%).\n',...
    mfilename,nNegFit,nNegFit./numel(wProcess).*100);
return;



%-----------------------------------
% make sure that the pRF can only be moved "step" away from original
% position "startParams" - for the one Gaussian model
function [C, Ceq]=distanceCon(x,startParams,step)
Ceq = [];
dist = x([1 2])-startParams([1 2]);
C = norm(dist) - step;
return;
%-----------------------------------

