function model = rmGridFit_twoPredictors(model,prediction,prediction2,data,params,t)
% rmGridFit_oneGaussian - core of one Gaussian fit
%
% model = rmGridFit_oneGaussian(model,prediction,data,params);
%
% 2008/01 SOD: split of from rmGridFit.

% input check 
if nargin < 4,
    error('Not enough arguments');
end

% some variables we need
%rssinf         = inf(size(data(1,:)),'single');
trends         = t.trends;
if length(params.stim)>1 && params.seperateRunBetas;
    t_id           = t.dcid+(2*length(params.stim));
else
    t_id           = t.dcid+2; 
end

% we compute mean rss but we need sum rss (old convention)
model.rss=single(model.rss./(size(prediction,1)-size(trends,2)+2));  

%-----------------------------------
%--- fit different receptive fields profiles
%--- another loop --- and a slow one too
%-----------------------------------
tic; progress = 0;

warning('off', 'MATLAB:lscov:RankDefDesignMat')
if ~isfield(params, 'seperateRunBetas')
    params.seperateRunBetas=0;
end

for n=1:numel(params.analysis.x0),
    %-----------------------------------
    % progress monitor (10 dots) and time indicator
    %-----------------------------------
    if floor(n./numel(params.analysis.x0).*10)>progress,
        if progress==0,
            % print out estimated time left
            esttime = toc.*10;
            if floor(esttime./3600)>0,
                fprintf(1,'[%s]:Estimated processing time: %d hours.\t(%s)\n',...
                    mfilename, ceil(esttime./3600), datestr(now));
            else
                fprintf(1, '[%s]:Estimated processing time: %d minutes.\t(%s)\n',...
                    mfilename, ceil(esttime./60), datestr(now));
            end;
            fprintf(1,'[%s]:Grid (x,y,sigma) fit:',mfilename);drawnow;
        end;
        % progress monitor
        fprintf(1,'.');drawnow;
        progress = progress + 1;
    end;

    %-----------------------------------
    %--- now apply glm to fit RF
    %-----------------------------------
    % minimum RSS fit
   nkeep=zeros(1,size(data,2));
   if length(params.stim)>1 && params.seperateRunBetas;
       
       framesPerScan=size(trends,1)/length(params.stim);
%        predictors=zeros(size(prediction,1),length(params.stim));
%        for stimRun=1:length(params.stim)
%            predictors(((stimRun-1)*framesPerScan+1):(framesPerScan*stimRun), stimRun)=prediction(((stimRun-1)*framesPerScan+1):(framesPerScan*stimRun),n);
%        end
%        X    = [predictors trends];
%        [b,ci,rss]    = lscov(X,data); 
%        nkeep   = false(size(b(1,:))); %b(1,:)<0 | b(2,:)<0 | b(3,:)<0 | b(4,:)<0 | b(5,:)<0;
       b=zeros(length(params.stim)*2+size(trends,2), size(data,2));
       rsstotal=zeros(length(params.stim), size(data,2));
       nkeep=zeros(1, size(data,2));
       for stimRun=1:length(params.stim)
          X    = [prediction(((stimRun-1)*framesPerScan+1):(framesPerScan*stimRun),n) prediction2(((stimRun-1)*framesPerScan+1):(framesPerScan*stimRun),n) trends(((stimRun-1)*framesPerScan+1):(framesPerScan*stimRun), stimRun)];
          [btmp,ci,rsstmp]    = lscov(X,data(((stimRun-1)*framesPerScan+1):(framesPerScan*stimRun),:));
          %set negative betas to zero and evaluate
          negb1=btmp(1,:)<0 & btmp(2,:)>0;
          if any(negb1)
              X= [prediction2(((stimRun-1)*framesPerScan+1):(framesPerScan*stimRun),n) trends(((stimRun-1)*framesPerScan+1):(framesPerScan*stimRun), stimRun)];
              [btmp2,ci,rsstmp2]    = lscov(X,data(((stimRun-1)*framesPerScan+1):(framesPerScan*stimRun),negb1));
              btmp(2:end, negb1)=btmp2;
              rsstmp(negb1)=rsstmp2;
          end
          negb2=btmp(1,:)>0 & btmp(2,:)<0;
          if any(negb2)
              X= [prediction(((stimRun-1)*framesPerScan+1):(framesPerScan*stimRun),n) trends(((stimRun-1)*framesPerScan+1):(framesPerScan*stimRun), stimRun)];
              [btmp2,ci,rsstmp2]    = lscov(X,data(((stimRun-1)*framesPerScan+1):(framesPerScan*stimRun),negb2));
              btmp([1 3:end], negb2)=btmp2;
              rsstmp(negb2)=rsstmp2;
          end
          negb3=btmp(1,:)<0 & btmp(2,:)<0;  % EVI this is the issue you were thinking about probably: if you set b(1,:) to 0, e.g., in negb1, it is not supposed to find it as negb3...
          btmp(1,negb1)=0;
          btmp(2,negb2)=0;
          if any(negb3)
              btmp(1:2,negb3)=0;
              X= [trends(((stimRun-1)*framesPerScan+1):(framesPerScan*stimRun), stimRun)];
              [btmp2,ci,rsstmp2]    = lscov(X,data(((stimRun-1)*framesPerScan+1):(framesPerScan*stimRun),negb3));
              btmp(3:end, negb3)=btmp2;
              rsstmp(negb3)=rsstmp2;
          end
          rsstotal(stimRun,:)=rsstmp;
          b(stimRun,:)=btmp(1,:);
          b((stimRun+length(params.stim)),:)=btmp(2,:);
          b((stimRun+(2*length(params.stim))),:)=btmp(3,:);
       end
       rss=mean(rsstotal,1);
       %nkeep   = b(1,:)<0 | b(2,:)<0 | b(3,:)<0 | b(4,:)<0 | b(5,:)<0; false(size(b(1,:))); %
   else
       X   = [prediction(:,n)   prediction2(:,n)    trends];
       
       % This line takes up 30% of the time
       % lscov takes as long as the pinv method but provides the rss as well...
       
       [b,ci,rss]    = lscov(X,data);
       %b     = pinv(X)*data;
       
       negb1=b(1,:)<0 & b(2,:)>0;
       if any(negb1)
           X   = [prediction2(:,n) trends];
           [b2,ci,rss2] =lscov(X,data(:,negb1)); % ?
           b(2:end, negb1)=b2;
           rss(negb1)=rss2;
       end
       
       negb2=b(1,:)>0 & b(2,:)<0;
       if any(negb2)
           X   = [prediction(:,n)   trends];
           [b2,ci,rss2]    =lscov(X,data(:,negb2));
           b([1 3:end], negb2)=b2;
           rss(negb2)=rss2;
       end
       
       negb3=b(1,:)<0 & b(2,:)<0;
       b(1,negb1)=0;
       b(2,negb2)=0;
       if any(negb3)
           b(1:2,negb3)=0;
           X   = [trends];
           [b2,ci,rss2]    = lscov(X,data(:,negb3));
           b(3:end, negb3)=b2;
           rss(negb3)=rss2;
       end
       
       % Compute RSS only for positive fits. The basic problem is
       % that if you have two complementary locations, you
       % could fit with a postive beta on the one that drives the signal or a
       % negative beta on the portion of the visual field that never sees the
       % stimulus. This would produce the same prediction. We don't like that
      % nkeep   = b(1,:)<0 |  b(2,:)<0; % Now we only set the negative fits to inf.
       
       
    end
    % To save time limit the rss computation to those we care about.
    % This line is takes up 60% of the time.... (replaced by lscov)
    if sum(nkeep)>0
        rss(nkeep) = inf('single');
    end
    
    %-----------------------------------
    %--- store data with lower rss
    %-----------------------------------
    minRssIndex = rss < model.rss;
    
    % now update
    model.x0(minRssIndex)       = params.analysis.x0(n);
    model.y0(minRssIndex)       = params.analysis.y0(n);
    model.s(minRssIndex)        = params.analysis.sigmaMajor(n);
    model.s_major(minRssIndex)        = params.analysis.sigmaMajor(n);
    model.s_minor(minRssIndex)        = params.analysis.sigmaMajor(n);
    model.s_theta(minRssIndex)        = params.analysis.theta(n);
    model.rss(minRssIndex)      = rss(minRssIndex);
    if length(params.stim)>1 && params.seperateRunBetas;
        model.b([1:(2*length(params.stim)) t_id],minRssIndex) = b(:,minRssIndex); 
    else
        model.b([1 2 t_id],minRssIndex) = b(:,minRssIndex);
    end
end;

%warning('on', 'MATLAB:lscov:RankDefDesignMat')

% Under some conditions, the grid fit never returns an acceptable fit, For
% example for onegaussian fits with data driven DC component, when the DC
% is artificially high. In this case some of the rss values remain Inf,
% which fails to interpolate and compute correct variance explained values.
% So we check it here and reset any Inf (bad fits) to rawrss, so the
% variance explained will be 0.
model.rss(model.rss==Inf)=model.rawrss(model.rss==Inf);

% Correct lscov. It returns the mean rss. To maintain compatibility with the
% sum rss this function expects, we have to multiply by the divisor. See
% the lscov docs for details.
model.rss=single(model.rss.*(size(prediction,1)-size(trends,2)+2));  %I think it's +2, but check

% end time monitor
et  = toc;
if floor(esttime/3600)>0,
    fprintf(1,'Done[%d hours].\t(%s)\n', ceil(et/3600), datestr(now));
else
    fprintf(1,'Done[%d minutes].\t(%s)\n', ceil(et/60), datestr(now));
end;
drawnow;
return;


