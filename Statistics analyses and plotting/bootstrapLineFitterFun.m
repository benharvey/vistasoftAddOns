function [ypoints, yfitParams, b_upper, b_lower]=bootstrapLineFitterFun(x,y,ve, xpoints)
    x = x(isfinite(ve)); y = y(isfinite(ve)); ve = ve(isfinite(ve));
    iis=1:length(x);
    [B] = bootstrp(1000,@(iis) lineFitterFun(iis,x,y,ve),[1:length(x)]');
    B = B';
    roi.p=B;
    pct1 = 100*0.05/2;
    pct2 = 100-pct1;
    b_lower = prctile(B',pct1);
    b_upper = prctile(B',pct2);
    slopeCI=[b_lower(1) b_upper(1)]
    yfitParams=prctile(B', 50);
    %y2fit = polyval(roi.p2(:),x2fit);
    keep1 = B(1,:)>b_lower(1) &  B(1,:)<b_upper(1);
    keep2 = B(2,:)>b_lower(2) &  B(2,:)<b_upper(2);
    keep = keep1 & keep2;

    ypoints=yfitParams(2)+xpoints.*yfitParams(1);
    fits = [xpoints' ones(size(xpoints'))]*B(:,keep);
    b_upper = max(fits,[],2);
    b_lower = min(fits,[],2);
return



