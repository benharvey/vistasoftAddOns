function [data] = PlotRoiDistanceSizeDataAllBin(data, veThresh, showPlots, binSteps, meanThresh, meanData)
%Makes various plots from the output of RoiDistanceRatio

if ~exist('veThresh', 'var') || isempty(veThresh)
    veThresh=0;
end
if ~exist('showPlots', 'var') || isempty(showPlots)
    showPlots=[1 1 0 1 1 0];
end
if ~exist('binSteps', 'var') || isempty(binSteps)
    binSteps=28;
end
if ~exist('meanThresh', 'var') || isempty(meanThresh)
    meanThresh=0;
end
if ~exist('meanData', 'var') || isempty(meanData)
    if isfield(data{1}, 'meanSignal')
        meanData=data{1}.meanSignal;
    else
        meanData=ones(size(data{1}.x0s));
    end
end
colors={'k', 'r', 'g'};

for n=1:length(data)
    data{n}.veIndices=data{n}.ves>=veThresh & data{n}.x0s>(0.05) & data{n}.x0s<1.4 & data{n}.sigmas<1.29 & meanData>=meanThresh;
    %data{n}.x0fit=linreg(data{n}.ratio(data{n}.veIndices), data{n}.x0s(data{n}.veIndices), data{n}.ves(data{n}.veIndices));
    data{n}.sigmafit=linreg(data{n}.ratio(data{n}.veIndices), data{n}.sigmas(data{n}.veIndices), data{n}.ves(data{n}.veIndices));
end
b=linspace(0.05, 1.4, binSteps);



%size vs mean distance

%xpoints=0:0.01:data{1}.meanDist;
bins=1:2:ceil(data{1}.meanDist);
xpoints=min(bins):0.01:max(bins);
for n=1:length(data)
%     %Quick and dirty way, evenly-spaced bins at odd intervals
%     data{n}.x{4}=b.*data{n}.meanDist;
%     data{n}.y{4}=data{n}.y{3};
%     data{n}.ysterr{4}=data{n}.ysterr{3};
%     [ypoints{n}, data{n}.logLineFit, bupper{n}, blower{n}]=bootstrapLogLineFitter(data{n}.x{4},data{n}.y{4},1./data{n}.ysterr{4}(1,:), xpoints);

%   Slower, but with bins at regular distance intervals
    data{n}.y{4}=nan(size(bins));
    data{n}.ysterr{4}=nan(2,length(bins));
    for bCount=1:length(bins)
        bii = data{n}.ratio.*data{1}.meanDist> bins(bCount)-(bins(2)-bins(1))./2 & ...
            data{n}.ratio.*data{1}.meanDist< bins(bCount)+(bins(2)-bins(1))./2 & data{n}.veIndices;
        if any(bii) && sum(data{n}.ves(bii))>0
            s=wstat((data{n}.x0s(bii)), data{n}.ves(bii), 1.78.^2);
            data{n}.y{4}(bCount)=s.mean;
            data{n}.ysterr{4}(:,bCount)=s.sterr;
        end
    end
    data{n}.x{4}=bins;

    [ypoints{n}, data{n}.linLineFit, bupper{n}, blower{n}]=bootstrapLineFitter(data{n}.x{4},data{n}.y{4},1./data{n}.ysterr{4}(1,:), xpoints);
    
    %Permutation testing
    ydat=data{n}.y{4}(~isnan(data{n}.y{4}));
    xdat=data{n}.x{4}(~isnan(data{n}.y{4}));
    fitDist=zeros(2,10000);
    for m=1:10000
        yshuffle=ydat(randperm(length(ydat)));
        fitDist(:,m)=linreg(xdat, yshuffle);
    end
    CI=[prctile(fitDist(2,:), 2.5) prctile(fitDist(2,:), 97.5)]
    max(fitDist(2,:))
    measure=linreg(xdat, ydat);
    measure(2)
    
    %To compare log and linear fits
%     if n==1
%     [tmp1, data{n}.logLineFit, tmp2,tmp3]=bootstrapLogLineFitter(data{n}.x{4},data{n}.y{4},1./data{n}.ysterr{4}(1,:), xpoints);
%     data{n}.rawrss=var(data{n}.y{4});
%     data{n}.rsslin=var(data{n}.y{4}-(data{n}.x{4}.*data{n}.linLineFit(1)+data{n}.linLineFit(2)));
%     data{n}.rsslog=var(data{n}.y{4}-exp(data{n}.x{4}.*data{n}.logLineFit(1)+data{n}.logLineFit(2)));
%     data{n}.linresid=(data{n}.y{4}-(data{n}.x{4}.*data{n}.linLineFit(1)+data{n}.linLineFit(2)));
%     data{n}.logresid=(data{n}.y{4}-exp(data{n}.x{4}.*data{n}.logLineFit(1)+data{n}.logLineFit(2)));
%     end
end

if showPlots(1)==1
    figure;
    hold on;
    for n=1:length(data)
        errorbar(data{n}.x{4},data{n}.y{4},data{n}.ysterr{4}(1,:),data{n}.ysterr{4}(1,:),'ko',...
            'MarkerFaceColor',colors{n},'MarkerSize',8);
        %plot(data{n}.ratio(data{n}.veIndices).*data{n}.meanDist, exp(data{n}.x0s(data{n}.veIndices)), strcat(colors{n}, 'o'), 'MarkerSize',3, 'MarkerFaceColor', colors{n})
        hold on; plot(xpoints, ypoints{n}, colors{n},'LineWidth',4);
        hold on; plot(xpoints, bupper{n}, colors{n});
        hold on; plot(xpoints, blower{n}, colors{n});
        %hold on; plot(xpoints, ypoints{n}, colors{n});
    end
    axis([0 5*ceil(data{n}.meanDist/5) 0 1.401]);
    axis square
end

%Summary plot, linked points, no errorbars
if showPlots(2)==1
    figure;
    hold on;
    for n=1:length(data)
        if n==1
        errorbar(data{n}.x{4},data{n}.y{4},data{n}.ysterr{4}(1,:),data{n}.ysterr{4}(1,:),'ko',...
            'MarkerFaceColor',colors{n},'MarkerSize',8);
        %plot(data{n}.ratio(data{n}.veIndices).*data{n}.meanDist, exp(data{n}.x0s(data{n}.veIndices)), strcat(colors{n}, 'o'), 'MarkerSize',3, 'MarkerFaceColor', colors{n})
        hold on; plot(xpoints, ypoints{n}, colors{n},'LineWidth',4);
        hold on; plot(xpoints, bupper{n}, colors{n});
        hold on; plot(xpoints, blower{n}, colors{n});

        else
            hold on; plot(data{n}.x{4},data{n}.y{4}, colors{n},'LineWidth',4);    
        end
    end
    axis([0 5*ceil(data{n}.meanDist/5) 0 1.401]);
    axis square;
end

%Cortical magnification factor

if showPlots(3)==1
    for n=1%:length(data)
            figure;
            hold on;
            xaxis=mean([ypoints{n}(1:end-1); ypoints{n}(2:end)],1);
     %         %Smooth curves, derivitive of dist vs num
            plot(xaxis, 1./(diff(ypoints{n})*100), colors{n},'LineWidth',4);
            
            xp2=-25:0.01:50;
            allPlots2=exp([xp2' ones(size(xp2'))]*B{n});
            diffs2=diff(allPlots2)*100;
            xaxis2=(allPlots2(1:end-1,:)+allPlots2(2:end,:))/2;
            minCI=zeros(size(xaxis));
            maxCI=minCI;
            for index=1:length(xaxis)
                bii=xaxis2>ypoints{n}(index) & xaxis2<ypoints{n}(index+1);
                minCI(index)=prctile(diffs2(bii),2.5);
                maxCI(index)=prctile(diffs2(bii),97.5);
            end
            plot(xaxis, 1./minCI, colors{n});
            plot(xaxis, 1./maxCI, colors{n});
            axis([0 6 0 20])
            axis square
            
            figure; hold on;
            plot(xaxis, (diff(ypoints{n})*100), colors{n},'LineWidth',4);
            plot(xaxis, minCI, colors{n});
            plot(xaxis, maxCI, colors{n});   
            
            
            
    end
end

%fwhm vs number
if showPlots(4)==1;
xpoints=[];
for n=1%:length(data)
    %     fwhms{n}=data{n}.sigmas.*(2*sqrt(2*log(2)));
    %     fwhms{n}=exp(data{n}.x0s+fwhms./2)-exp(data{n}.x0s-fwhms./2);
    if isfield(data{n}, 'sigmas2')
        fwhms{n}=rmGetDoGFWHM_size(data{n});
    else
        fwhms{n}=data{n}.sigmas.*(2*sqrt(2*log(2)));
    end
%     fwhmLog{n}=(data{n}.x0s+fwhms{n}./2)-(data{n}.x0s-fwhms{n}./2);
%     fwhms{n}=exp(data{n}.x0s+fwhms{n}./2)-exp(data{n}.x0s-fwhms{n}./2);
    b=linspace(0.1, 0.9, 17);
    data{n}.x{9}=b;
    data{n}.y{9}=nan(size(b));
    data{n}.ysterr{9}=nan(2,length(b));
%     data{n}.yLog{9}=nan(size(b));
%     data{n}.ysterrLog{9}=nan(2,length(b));
    for bCount=1:length(b)
        bii = (data{n}.x0s)> b(bCount)-(b(2)-b(1))./2 & ...
            (data{n}.x0s)< b(bCount)+(b(2)-b(1))./2 & ...
            fwhms{n}<=20 & data{n}.veIndices;
        if any(bii) && sum(data{n}.ves(bii))>0
            s=wstat(fwhms{n}(bii), data{n}.ves(bii), 1.77^2);
            if s.sterr<1 && s.sterr>0
                data{n}.y{9}(bCount)=s.mean;
                data{n}.ysterr{9}(:,bCount)=s.sterr;
            end
%             s=wstat(fwhmLog{n}(bii), data{n}.ves(bii), 4);
%             data{n}.yLog{9}(bCount)=s.mean;
%             data{n}.ysterrLog{9}(:,bCount)=s.sterr;
        end
    end
    xpoints{n}=min(b(isfinite(data{n}.y{9}))):0.01:max(data{n}.x{9}(isfinite(data{n}.y{9})));%max(exp(data{n}.x0s));
    
    %data{n}.fwhmnumfit=linreg(data{n}.x{9}(isfinite(1./data{n}.ysterr{9}(1,:))), data{n}.y{9}(isfinite(1./data{n}.ysterr{9}(1,:))), 1./data{n}.ysterr{9}(1,isfinite(1./data{n}.ysterr{9}(1,:))));
    [ypoints{n}, data{n}.fwhmNumFit, bupper{n}, blower{n}]=bootstrapLineFitter(data{n}.x{9},data{n}.y{9},1./data{n}.ysterr{9}(1,:), xpoints{n});
    %[ypoints{n}, data{n}.fwhmNumFit, bupper{n}, blower{n}]=bootstrapLineFitter(data{n}.x{9},data{n}.yLog{9},1./data{n}.ysterrLog{9}(1,:), xpoints{n});

    
    %Permutation testing
%     ydat=data{n}.y{9}(~isnan(data{n}.y{9}));
%     xdat=data{n}.x{9}(~isnan(data{n}.y{9}));
%     fitDist=zeros(2,10000);
%     for m=1:10000
%         yshuffle=ydat(randperm(length(ydat)));
%         fitDist(:,m)=linreg(xdat, yshuffle);
%     end
%     CI=[prctile(fitDist(2,:), 2.5) prctile(fitDist(2,:), 97.5)]
%     min(fitDist(2,:))
%     measure=linreg(xdat, ydat);
%     measure(2)
end
end
if showPlots(5)==1
    figure;
    hold on;
    for n=1%:length(data)
        errorbar(b,data{n}.y{9},data{n}.ysterr{9}(1,:),data{n}.ysterr{9}(1,:),'ko',...
            'MarkerFaceColor',colors{n},'MarkerSize',8);
        
        %plot(exp(data{n}.x0s(data{n}.veIndices)), fwhms(data{n}.veIndices), strcat(colors{n}, 'o'), 'MarkerSize',3, 'MarkerFaceColor', colors{n})
        %fitIndices=data{n}.veIndices & fwhms{n}<30 & data{n}.x0s<=7;

        hold on; plot(xpoints{n}, ypoints{n}, colors{n},'LineWidth',4);
        hold on; plot(xpoints{n}, bupper{n}, colors{n});
        hold on; plot(xpoints{n}, blower{n}, colors{n});
        %hold on; plot([0 7], [data{n}.fwhmnumfit(1), data{n}.fwhmnumfit(1)+data{n}.fwhmnumfit(2)*7], colors{n});
    end
    %     xpoints=0:0.01:1;
    %     ypoints=x0fit(1)+xpoints.*x0fit(2);
    %     hold on; plot(xpoints, exp(ypoints));
    axis([0 max(b)+0.1 0 2]);
    axis square
    xlabel('Preferred Object Size (deg diameter)')
    ylabel('Population Tuning Width (deg diameter)')
    
    %axis([0 max(b)+0.1 0 1.5]);
end

return

function [ypoints, yfitParams, b_upper, b_lower]=bootstrapLogLineFitter(x,y,ve, xpoints)
    x = x(isfinite(ve)); y = y(isfinite(ve)); ve = ve(isfinite(ve));
    iis=1:length(x);
    [B] = bootstrp(1000,@(iis) logLineFitter(iis,x,y,ve),[1:length(x)]');
    B = B';
    roi.p=B;
    pct1 = 100*0.05/2;
    pct2 = 100-pct1;
    b_lower = prctile(B',pct1);
    b_upper = prctile(B',pct2);
    yfitParams=prctile(B', 50);
    %y2fit = polyval(roi.p2(:),x2fit);
    keep1 = B(1,:)>b_lower(1) &  B(1,:)<b_upper(1);
    keep2 = B(2,:)>b_lower(2) &  B(2,:)<b_upper(2);
    keep = keep1 & keep2;

    ypoints=exp(yfitParams(2)+xpoints.*yfitParams(1));
    fits = exp([xpoints' ones(size(xpoints'))]*B(:,keep));
    b_upper = max(fits,[],2);
    b_lower = min(fits,[],2);
return

function [B,e]=logLineFitter(iis, x,y,ve)
%x = x(isfinite(ve)); y = y(isfinite(ve)); ve = ve(isfinite(ve));
x=x(iis); 
y=y(iis); 
ve=ve(iis);
options = optimset('MaxFunEvals',10000000);
[B,e] = fminsearch(@(z) mylogfit(z,x,y,ve),[0.2;1.3], options);
return

function e=mylogfit(z,x,y,ve)
e=sum(ve.*(y-(exp(z(1).*x+z(2)))).^2)./sum(ve);
return

function [ypoints, yfitParams, b_upper, b_lower]=bootstrapLineFitter(x,y,ve, xpoints)
    x = x(isfinite(ve)); y = y(isfinite(ve)); ve = ve(isfinite(ve));
    iis=1:length(x);
    [B] = bootstrp(1000,@(iis) lineFitter(iis,x,y,ve),[1:length(x)]');
    B = B';
    roi.p=B;
    pct1 = 100*0.05/2;
    pct2 = 100-pct1;
    b_lower = prctile(B',pct1);
    b_upper = prctile(B',pct2);
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

function [B,e]=lineFitter(iis, x,y,ve)
%x = x(isfinite(ve)); y = y(isfinite(ve)); ve = ve(isfinite(ve));
x=x(iis); 
y=y(iis); 
ve=ve(iis);
options = optimset('MaxFunEvals',10000000);
[B,e] = fminsearch(@(z) mylinfit(z,x,y,ve),[0.2;1.3], options);
return

function e=mylinfit(z,x,y,ve)
e=sum(ve.*(y-((z(1).*x+z(2)))).^2)./sum(ve);
return


function [ypoints, yfitParams, b_upper, b_lower] = bootstrapCmfFitter(x, y, ve, xpoints)
    x = x(isfinite(ve)); y = y(isfinite(ve)); ve = ve(isfinite(ve));
    iis=1:length(x);
    [B] = bootstrp(1000,@(iis) cmfLineFitter(iis,x,y,ve),[1:length(x)]');
    B = B';
    roi.p=B;
    pct1 = 100*0.05/2;
    pct2 = 100-pct1;
    b_lower = prctile(B',pct1);
    b_upper = prctile(B',pct2);
    yfitParams=prctile(B', 50);
    %y2fit = polyval(roi.p2(:),x2fit);
    keep1 = B(1,:)>b_lower(1) &  B(1,:)<b_upper(1);
    keep2 = B(2,:)>b_lower(2) &  B(2,:)<b_upper(2);
    keep = keep1 & keep2;

    ypoints=1./(xpoints.*yfitParams(1)+yfitParams(2));
    fits = 1./([xpoints' ones(size(xpoints'))]*B(:,keep));
    b_upper = max(fits,[],2);
    b_lower = min(fits,[],2);
return


function [B, e]=cmfLineFitter(ii,x,y,ve)
x = x(ii); y = y(ii); ve = ve(ii);
options = optimset('MaxFunEvals',10000000);
[B, e] = fminsearch(@(z) mycmffit(z,x,y,ve),[0.05;0.2], options);
return

function e=mycmffit(z,x,y,ve)
e=sum(ve.*(y-(1./(z(1).*x+z(2)))).^2)./sum(ve);
return