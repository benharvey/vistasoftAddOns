function e = rmModelSearchFit_xGaussianYCompressive(p,Y,Xv,Yv,stim,t)
% rmModelSearchFit_oneGaussian - actual fit function of rmSearchFit
%
% error = rmModelSearchFit(p,Y,trends,Xgrid,YGrid,stimulusMatrix);
%
% Basic barebones fit of a single time-series. Error is returned in
% percentage: 100% is RSS of unfitted time-series. This way we can quantify
% the improvement of the fit independend of the variation in the raw
% time-series.
%
% 2006/06 SOD: wrote it.
% 2006/12 SOD: modifications for fmincon, this is litterally called >10000
% times so we cut every corner possible. 
% 2010/02 SOD: evaluated lscov this did not improve performance here (using
% profiler)

% make RF (taken from rfGaussian2d)
Xv = Xv - p(1);   % positive x0 moves center right
RF = exp( -.5 *((Xv.*Xv ./ p(3)).^2)).*((Yv.^p(2))./Yv);
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
b = pinvX*Y;

% compute residual sum of squares (e)
% e = norm(Y - X*abs(b));
if b(1)>0,
    e = norm(Y - X*b);
else
    e = norm(Y).*(1+sum(abs(b(1))));
end
return;
