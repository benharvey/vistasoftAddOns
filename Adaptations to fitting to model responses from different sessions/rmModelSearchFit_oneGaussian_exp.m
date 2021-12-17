function e = rmModelSearchFit_oneGaussian_exp(p,Y,Xv,Yv,stim,t)
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
Yv = Yv - p(2);   % positive y0 moves center up
RF = exp( (Yv.*Yv + Xv.*Xv) ./ (-2.*(p(3).^2)) );

%Add stimulus exponent
periods=[0.05:0.05:1 2.1 0.1 0.6 0.9 1 0.15 0.65 0.85 1 0.2 0.7 0.8 1 0.25 0.75 1 0.3 0.7 0.8 1 0.35 0.65 0.85 1 0.4 0.6 0.9 1 0.45 0.55 0.95 1 0.5 1 0.55 1 0.6 1 0.65 1 0.7 1 0.75 1 0.8 1 0.85 1 0.9  1 0.95 1 1 2.1];
freq=5./periods;
scale=1./(freq'.^p(end));
scale=scale.*freq';
RF=RF./scale;

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
