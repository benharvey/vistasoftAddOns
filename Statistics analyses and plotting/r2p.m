function [p] = r2p(r, n)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
m=[2];
lowerhalf = (tril(ones(m),-1)>0);
   rv = r;%(lowerhalf);
   if length(n)>1
      nv = n(lowerhalf);
   else
      nv = n;
   end

   % Tstat = +/-Inf and p = 0 if abs(r) == 1, NaN if r == NaN.
   Tstat = rv .* sqrt((nv-2) ./ (1 - rv.^2));
   %p = zeros(m);
   p = 2*tpvalue(-abs(Tstat),nv-2);
   %p = p + p' + diag(diag(r)); % Preserve NaNs on diag.
end


function p = tpvalue(x,v)
%TPVALUE Compute p-value for t statistic.

normcutoff = 1e7;
if length(x)~=1 && length(v)==1
   v = repmat(v,size(x));
end

% Initialize P.
p = NaN(size(x));
nans = (isnan(x) | ~(0<v)); % v == NaN ==> (0<v) == false

% First compute F(-|x|).
%
% Cauchy distribution.  See Devroye pages 29 and 450.
cauchy = (v == 1);
p(cauchy) = .5 + atan(x(cauchy))/pi;

% Normal Approximation.
normal = (v > normcutoff);
p(normal) = 0.5 * erfc(-x(normal) ./ sqrt(2));

% See Abramowitz and Stegun, formulas 26.5.27 and 26.7.1.
gen = ~(cauchy | normal | nans);
p(gen) = betainc(v(gen) ./ (v(gen) + x(gen).^2), v(gen)/2, 0.5)/2;

% Adjust for x>0.  Right now p<0.5, so this is numerically safe.
reflect = gen & (x > 0);
p(reflect) = 1 - p(reflect);

% Make the result exact for the median.
p(x == 0 & ~nans) = 0.5;
end