function e = rmModelSearchFit_oneOvalGaussian_seperateBetas(p,Y,Xv,Yv,stim,t,numRuns)
% rmModelSearchFit_oneOvalGaussian - actual fit function of rmSearchFit
%
% error = rmModelSearchFit_oneOvalGaussian(p,Y,trends,Xgrid,YGrid,stimulusMatrix);
%
% Basic barebones fit of a single time-series. Error is returned in
% percentage: 100% is RSS of unfitted time-series. This way we can quantify
% the improvement of the fit independend of the variation in the raw
% time-series.
%
% 2006/06 SOD: wrote it.
% 2006/12 SOD: modifications for fmincon, this is litterally called >10000
% times so we cut every corner possible. 

% make RF (taken from rfGaussian2d)
Xv = Xv - p(1);   % positive x0 moves center right
Yv = Yv - p(2);   % positive y0 moves center up

Xold = Xv;
Yold = Yv;
Xv = Xold .* cos(p(5)) - Yold .* sin(p(5));
Yv = Xold .* sin(p(5)) + Yold .* cos(p(5));

% make gaussian on current grid
RF = exp( -.5 * ((Yv ./ p(3)).^2 + (Xv ./ p(4)).^2));

% make prediction (taken from rfMakePrediction)
X = [stim*RF t];

% fit - inlining pinv
%b = pinv(X)*Y; 
[U,S,V] = svd(X,0);
s = diag(S); 
tol = numel(X) * eps(max(s));
r = sum(s > tol);
if (r == 0)
    pinvX = zeros(size(X'));
else
    s = diag(ones(r,1)./s(1:r));
    pinvX = V(:,1:r)*s*U(:,1:r)';
end

%To fit separate betas per run
residual=zeros(size(Y));
framesPerStim=length(Y)/numRuns;
b=zeros(size(X,2), numRuns);
for n=1:numRuns
    indices=((n-1)*framesPerStim+1):(n*framesPerStim);
    b(:,n) = pinvX(:,indices)*Y(indices);
    residual(indices)=Y(indices)-X(indices,:)*b(:,n);
end

% compute residual sum of squares (e)
% e = norm(Y - X*abs(b));
if all(b(1,:)>0)
    e = norm(residual);
else
    e = norm(Y).*(1+sum(abs(b(1,:))));
end
return;
