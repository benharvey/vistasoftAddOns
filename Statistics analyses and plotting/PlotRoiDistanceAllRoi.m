function [dataAll] = PlotRoiDistanceAllRoi(data, veThresh, showPlots, binSteps, meanThresh)
%Makes various plots from the output of RoiDistanceRatio. THis version is
%made to use Log2Lin models. It's a simpler version of
%PlotRoiDistanceDataAllBin, giving fewer and simpler outputs

if ~exist('veThresh', 'var') || isempty(veThresh)
    veThresh=0;
end
if ~exist('showPlots', 'var') || isempty(showPlots)
    showPlots=[0 0 0 1 1 0];
end
if ~exist('binSteps', 'var') || isempty(binSteps)
    binSteps=11;
end
if ~exist('meanThresh', 'var') || isempty(meanThresh)
    meanThresh=zeros(size(data,1), 1);
end

colors={'k', 'r', 'g', 'b', 'c','m'};

dataAll=data;

for whichRoi=1:size(dataAll, 1)
    data=dataAll(whichRoi,:);
    
    
     for n=1:length(data)
        %Just run this part if data are in log space
        data(n).x0sLog=data(n).x0s;
        data(n).sigmasLog=data(n).sigmas;
        fwhms=data(n).sigmas.*(2*sqrt(2*log(2)));
        data(n).sigmas=exp(data(n).x0s+fwhms./2)-exp(data(n).x0s-fwhms./2);
        data(n).x0s=exp(data(n).x0s);
        
        
        
        data(n).veIndices=data(n).ves>=veThresh & data(n).x0s>1.02 & data(n).x0s<7.0 & data(1).meanSignal>=meanThresh(whichRoi);
        %data(n).x0fit=linreg(data(n).ratio(data(n).veIndices), data(n).x0s(data(n).veIndices), data(n).ves(data(n).veIndices));
        %data(n).sigmafit=linreg(data(n).ratio(data(n).veIndices), data(n).sigmas(data(n).veIndices), data(n).ves(data(n).veIndices));
    end
    b=linspace(0.5/binSteps, 1-(0.5/binSteps), binSteps);
    
    %number vs mean distance
    
    %xpoints=0:0.01:data{1}.meanDist;
    bins=(1:2:ceil(data(n).meanDist))./data(n).meanDist;
    for n=1%:length(data)
        %     %Quick and dirty way, evenly-spaced bins at odd intervals
        %     data(n).x{4}=b.*data(n).meanDist;
        %     data(n).y{4}=data(n).y{3};
        %     data(n).ysterr{4}=data(n).ysterr{3};
        %     [data(n).ypoints, data(n).logLineFit, data(n).bupper, data(n).blower]=bootstrapLogLineFitter(data(n).x{4},data(n).y{4},1./data(n).ysterr{4}(1,:), xpoints);
        
        %   Slower, but with bins at regular distance intervals
        data(n).y{4}=nan(size(bins));
        data(n).ysterr{4}=nan(2,length(bins));
        for bCount=1:length(bins)
            bii = data(n).ratio> bins(bCount)-0.5/binSteps & ...
                data(n).ratio< bins(bCount)+0.5/binSteps & data(n).veIndices;
            if any(bii)
                s=wstat(data(n).x0s(bii), data(n).ves(bii), 1.77^2);
                data(n).y{4}(bCount)=s.mean;
                data(n).ysterr{4}(:,bCount)=s.sterr;
            end
        end
        data(n).x{4}=bins;
        
        xpoints=min(bins(~isnan(data(n).y{4}))):0.01:max(bins(~isnan(data(n).y{4})));
        data(n).xpoints=xpoints;
        try
        [data(n).ypoints, data(n).logLineFit, data(n).bupper, data(n).blower, data(n).B]=bootstrapLogLineFitter(data(n).x{4},data(n).y{4},1./data(n).ysterr{4}(1,:), xpoints);
        catch
           n 
        end
        
        %Statistical test of progression with distance
            permutations=10000;
            ydat=log(data(n).y{4}(~isnan(data(n).y{4})));
            xdat=data(n).x{4}(~isnan(data(n).y{4}));
            fitDist=zeros(2,permutations);
            for m=1:permutations
                yshuffle=ydat(randperm(length(ydat)));
                fitDist(:,m)=linreg(xdat, yshuffle);
            end
            CI=[prctile(fitDist(2,:), 2.5) prctile(fitDist(2,:), 97.5)];
            max(fitDist(2,:));
            measure=linreg(xdat, ydat);
            measure(2);
            [whichRoi, n]
            if measure(2)<max(fitDist(2,:))
               tmp= fitDist(2,:);
               tmp=sort(tmp);
               tmp=tmp>measure(2);
               p=(permutations-find(tmp, 1, 'first')+1)./permutations
            else
                p=0
            end
        
        %To compare log and linear fits
        if n==1
            [tmp1, data(n).linLineFit, tmp2,tmp3]=bootstrapLineFitter(data(n).x{4},data(n).y{4},1./data(n).ysterr{4}(1,:), xpoints);
            data(n).rawrss=var(data(n).y{4});
            data(n).rsslin=var(data(n).y{4}-(data(n).x{4}.*data(n).linLineFit(1)+data(n).linLineFit(2)));
            data(n).rsslog=var(data(n).y{4}-exp(data(n).x{4}.*data(n).logLineFit(1)+data(n).logLineFit(2)));
            data(n).linresid=(data(n).y{4}-(data(n).x{4}.*data(n).linLineFit(1)+data(n).linLineFit(2)));
            data(n).logresid=(data(n).y{4}-exp(data(n).x{4}.*data(n).logLineFit(1)+data(n).logLineFit(2)));
        end
    end
    
    if showPlots(1)==1
        figure;
        hold on;
        for n=1:length(data)
            errorbar(data(n).x{4}.*data(n).meanDist,data(n).y{4},data(n).ysterr{4}(1,:),data(n).ysterr{4}(1,:),strcat(colors{n},'o'),...
                'MarkerFaceColor',colors{n},'MarkerEdgeColor','k','MarkerSize',8);
            %plot(data(n).ratio(data(n).veIndices).*data(n).meanDist, exp(data(n).x0s(data(n).veIndices)), strcat(colors{n}, 'o'), 'MarkerSize',3, 'MarkerFaceColor', colors{n})

            hold on; plot(data(n).xpoints.*data(n).meanDist, data(n).ypoints, colors{n},'LineWidth',2);
            hold on; plot(data(n).xpoints.*data(n).meanDist, data(n).bupper, colors{n});
            hold on; plot(data(n).xpoints.*data(n).meanDist, data(n).blower, colors{n});
            
            %hold on; plot(xpoints, data(n).ypoints, colors{n});
        end
        axis([0 5*ceil(data(n).meanDist/5) 0 10]);
        axis square
    end
    
    %Summary plot, linked points, no errorbars
    if showPlots(2)==1
        figure;
        hold on;
        for n=1:length(data)
            if n==1
                errorbar(data(n).x{4}.*data(n).meanDist,data(n).y{4},data(n).ysterr{4}(1,:),data(n).ysterr{4}(1,:),strcat(colors{n},'o'),...
                'MarkerFaceColor',colors{n},'MarkerEdgeColor','k','MarkerSize',8);
                %plot(data(n).ratio(data(n).veIndices).*data(n).meanDist, exp(data(n).x0s(data(n).veIndices)), strcat(colors{n}, 'o'), 'MarkerSize',3, 'MarkerFaceColor', colors{n})

                hold on; plot(data(n).xpoints.*data(n).meanDist, data(n).ypoints, colors{n},'LineWidth',2);
                hold on; plot(data(n).xpoints.*data(n).meanDist, data(n).bupper, colors{n});
                hold on; plot(data(n).xpoints.*data(n).meanDist, data(n).blower, colors{n});
                
                
            else
                hold on; plot(data(n).x{4}.*data(n).meanDist,data(n).y{4}, colors{n},'LineWidth',2);
            end
        end
        axis([0 5*ceil(data(n).meanDist/5) 0 10]);
        axis square;
    end
    dataOut(whichRoi,:)=data;
end

dataAll=dataOut;

colors={[128 0 255]./255, [255 0 255]./255, [255 128 0]./255, [128 64 0]./255, [128 128 128]./255, [0 128 0]./255};

%fwhm vs number
if showPlots(4)==1
    for n=1:size(dataAll,1)
        dataAll(n,1).xpointsFwhm=[];
        %     fwhms=data(n).sigmas.*(2*sqrt(2*log(2)));
        %     fwhms=exp(data(n).x0s+fwhms./2)-exp(data(n).x0s-fwhms./2);
        
        b=linspace(1, 5, 17);
        dataAll(n,1).x{9}=b;
        dataAll(n,1).y{9}=nan(size(b));
        dataAll(n,1).ysterr{9}=nan(2,length(b));
%         dataAll(n,1).yLog{9}=nan(size(b));
%         dataAll(n,1).ysterrLog{9}=nan(2,length(b));
        for bCount=1:length(b)
            bii = dataAll(n,1).x0s> b(bCount)-(b(2)-b(1))./2 & ...
                dataAll(n,1).x0s< b(bCount)+(b(2)-b(1))./2 & ...
                dataAll(n,1).sigmas<=20 & dataAll(n,1).veIndices;
            if any(bii) && sum(dataAll(n,1).ves(bii))>0;
                s=wstat(dataAll(n,1).sigmas(bii), dataAll(n,1).ves(bii), 1.77^2);
                if s.sterr<3 && s.sterr>0
                    dataAll(n,1).y{9}(bCount)=s.mean;
                    dataAll(n,1).ysterr{9}(:,bCount)=s.sterr;
                    s=wstat(dataAll(n,1).sigmasLog(bii), dataAll(n,1).ves(bii), 1.77^2);
                    dataAll(n,1).y{10}(bCount)=s.mean;
                    dataAll(n,1).ysterr{10}(:,bCount)=s.sterr;
                end
                
            end
        end
        dataAll(n,1).xpointsFwhm=min(b(isfinite(dataAll(n,1).y{9}))):0.1:max(b(isfinite(dataAll(n,1).y{9})));%max(exp(data(n).x0s));
        
        
        %data(n).fwhmnumfit=linreg(dataAll(n,1).x{9}(isfinite(1./dataAll(n,1).y{sterr9}(1,:))), dataAll(n,1).y{9}(isfinite(1./dataAll(n,1).y{sterr9}(1,:))), 1./dataAll(n,1).y{sterr9}(1,isfinite(1./dataAll(n,1).y{sterr9}(1,:))));
        [dataAll(n,1).ypointsFwhm, dataAll(n,1).fwhmNumFit, dataAll(n,1).fwhmbupper, dataAll(n,1).fwhmblower]=bootstrapLineFitter(dataAll(n,1).x{9},dataAll(n,1).y{9},1./dataAll(n,1).ysterr{9}(1,:), dataAll(n,1).xpointsFwhm);
        
%         permutations=10000;
%         ydat=dataAll(n,1).y{9}(~isnan(dataAll(n,1).y{9}));
%         xdat=b(~isnan(dataAll(n,1).y{9}));
%         fitDist=zeros(2,permutations);
%         for m=1:permutations
%             yshuffle=ydat(randperm(length(ydat)));
%             fitDist(:,m)=linreg(xdat, yshuffle);
%         end
%         CI=[prctile(fitDist(2,:), 2.5) prctile(fitDist(2,:), 97.5)];
%         max(fitDist(2,:));
%         measure=linreg(xdat, ydat);
%         measure(2);
%         n
%         if measure(2)<max(fitDist(2,:))
%             tmp= fitDist(2,:);
%             tmp=sort(tmp);
%             tmp=tmp>measure(2);
%             p=(permutations-find(tmp, 1, 'first')+1)./permutations
%         else
%             p=0
%         end
    end
    
    if showPlots(5)==1
        figure;
        hold on;
        for n=1:size(dataAll,1)%:length(data)
            errorbar(b,dataAll(n,1).y{9},dataAll(n,1).ysterr{9}(1,:),dataAll(n,1).ysterr{9}(1,:),strcat('o'),...
                'Color', colors{n},'MarkerFaceColor',colors{n},'MarkerEdgeColor','k','MarkerSize',8);
            
            %plot(exp(data(n).x0s(data(n).veIndices)), fwhms(data(n).veIndices), strcat(colors{n}, 'o'), 'MarkerSize',3, 'MarkerFaceColor', colors{n})
            %fitIndices=data(n).veIndices & dataAll(n,1).sigmas<30 & data(n).x0s<=7;
            
            hold on; plot(dataAll(n,1).xpointsFwhm, dataAll(n,1).ypointsFwhm, 'Color', colors{n},'LineWidth',2);
            hold on; plot(dataAll(n,1).xpointsFwhm, dataAll(n,1).fwhmbupper, 'Color', colors{n});
            hold on; plot(dataAll(n,1).xpointsFwhm, dataAll(n,1).fwhmblower, 'Color', colors{n});
            %hold on; plot([0 7], [data(n).fwhmnumfit(1), data(n).fwhmnumfit(1)+data(n).fwhmnumfit(2)*7], colors{n});
        end
        %     xpoints=0:0.01:1;
        %     ypoints=x0fit(1)+xpoints.*x0fit(2);
        %     hold on; plot(xpoints, exp(ypoints));
        axis([0 6 0 20]);
        axis square
        
        %Summary plot
        figure;
        hold on;
        for n=1:size(dataAll,1)
            dataAll(n,1).x{9}=b(isfinite(dataAll(n,1).y{9}));
            plot(dataAll(n,1).x{9},dataAll(n,1).y{9}(isfinite(dataAll(n,1).y{9})), 'Color', colors{n},'LineWidth',2);
        end
        axis([0 6 0 20]);
        axis square
    end
    
end
    

% CMF
if showPlots(3)==1
    figure;
    hold on;
    for n=1:size(dataAll,1)

            xaxis{n}=mean([dataAll(n,1).ypoints(1:end-1); dataAll(n,1).ypoints(2:end)],1);
     %         %Smooth curves, derivitive of dist vs num
            plot(xaxis{n}, 1./(diff(dataAll(n,1).ypoints)*100/dataAll(n,1).meanDist), 'Color', colors{n},'LineWidth',2);

            xp2=-25:0.01:50;
            allPlots2=exp([xp2' ones(size(xp2'))]*dataAll(n,1).B);
            diffs2=diff(allPlots2)*100/dataAll(n,1).meanDist;
            xaxis2=(allPlots2(1:end-1,:)+allPlots2(2:end,:))/2;
            minCI{n}=zeros(size(xaxis{n}));
            maxCI{n}=minCI{n};
            for index=1:length(xaxis{n})
                bii=xaxis2>dataAll(n,1).ypoints(index) & xaxis2<dataAll(n,1).ypoints(index+1);
                minCI{n}(index)=prctile(diffs2(bii),2.5);
                maxCI{n}(index)=prctile(diffs2(bii),97.5);
            end
            fit=linreg(xaxis{n}, minCI{n});
            minCI{n}=fit(1)+fit(2)*xaxis{n};
            fit=linreg(xaxis{n}, maxCI{n});
            maxCI{n}=fit(1)+fit(2)*xaxis{n};
            plot(xaxis{n}, 1./minCI{n}, 'Color', colors{n});
            plot(xaxis{n}, 1./maxCI{n}, 'Color', colors{n});






            %plot(mean([data(n).bupper(1:end-1); data(n).bupper(2:end)],1), 1./(diff(data(n).bupper)*100), colors{n});
            %plot(mean([data(n).blower(1:end-1); data(n).blower(2:end)],1), 1./(diff(data(n).blower)*100), colors{n});
%         plot(data(n).x{6},data(n).y{6},'ko','MarkerFaceColor',colors{n},'MarkerSize',8);
%         hold on; plot(xpoints, ypointsDist{n}, colors{n},'LineWidth',4);
%         hold on; plot(xpoints, bupperDist{n}, colors{n});
%         hold on; plot(xpoints, blowerDist{n}, colors{n});
    end
    axis([0 7 0 20])
    axis square
    figure; hold on;
    for n=1:size(dataAll,1)    
            plot(xaxis{n}, (diff(dataAll(n,1).ypoints)*100/dataAll(n,1).meanDist), 'Color', colors{n},'LineWidth',2);
            plot(xaxis{n}, minCI{n}, 'Color', colors{n});
            plot(xaxis{n}, maxCI{n}, 'Color', colors{n});
    end
    axis square
    axis([0 7 0 0.8])
end

%Number of voxels preferring each number
if showPlots(6)==1
for n=1:size(dataAll,1)
    figure; hold on;
    for m=1:12
        count(m)=sum(dataAll(n,1).x0s>((m+1)/2) & dataAll(n,1).x0s<((m+2)/2));
    end
    bar((((1:12)+1.5)./2), count./(1.77^2*5))
    axis square
end

count=[]
for n=1:size(dataAll,1)
    figure; hold on;
    for m=1:6
        count(m)=sum(dataAll(n,1).x0s>((m)) & dataAll(n,1).x0s<((m+1)));
    end
    bar((1:6)+0.5, count./(1.77^2*5))
    plot((1:6)+0.5, count./(1.77^2*5), 'k')
    axis square
    axis([1 7 0 70])
end

%Separate subjects error bars.
if isfield(dataAll, 'subjectVoxels')
figure; hold on;
for n=1:size(dataAll,1)
    count=[];
    dataAll(n,1).subjectVoxels=[0 dataAll(n,1).subjectVoxels];
    dataAll(n,1).subjectVoxels=cumsum(dataAll(n,1).subjectVoxels);
    for sub=2:(length(dataAll(n,1).subjectVoxels));
        vox=(dataAll(n,1).subjectVoxels(sub-1)+1):dataAll(n,1).subjectVoxels(sub);
        for m=1:6
            count(m, sub-1)=sum(dataAll(n,1).x0s(vox)>((m)) & dataAll(n,1).x0s(vox)<((m+1)));
        end
    end
    dataAll(n,1).voxelCountY=count;
    dataAll(n,1).voxelCountX=1:6+0.5; 
    
    count=(count./repmat(sum(count, 1), [6 1])).*mean(sum(count, 1));
    errorbar((1:6)+0.5, mean(count, 2)./1.77^2, std(count, [], 2)./(1.77^2*sqrt(length(dataAll(n,1).subjectVoxels)-1)), 'k')
    axis square
    axis([1 7 0 100])
end
end




    %All ROIs num vs dist plot, normalized for distance
    %b=linspace(0.5/binSteps, 1-(0.5/binSteps), binSteps);
%     b=(1:2:ceil(data(n).meanDist))./data(n).meanDist;
%     xpoints=min(b):0.01:max(b);
        figure;
        hold on;
        for n=1:size(dataAll,1)
            xdata=dataAll(n,1).x{4};
            xpoints=dataAll(n,1).xpoints;
            xpoints=xpoints-min(xdata(~isnan(dataAll(n,1).y{4})));
            xdata=xdata-min(xdata(~isnan(dataAll(n,1).y{4})));
            xpoints=xpoints./max(xdata(~isnan(dataAll(n,1).y{4})));
            xdata=xdata./max(xdata(~isnan(dataAll(n,1).y{4})));
            xpoints=xpoints.*0.9;
            xdata=xdata.*0.9;
            xpoints=xpoints+0.05;
            xdata=xdata+0.05;
            
%             0ffset=-min(dataAll(n,1).x{4})+0.05;
%             scale=(max(dataAll(n,1).x{4})-min(dataAll(n,1).x{4}))./1.1;
            errorbar(xdata,dataAll(n,1).y{4},dataAll(n,1).ysterr{4}(1,:),dataAll(n,1).ysterr{4}(1,:),strcat('o'),...
                'Color', colors{n}, 'MarkerFaceColor',colors{n},'MarkerEdgeColor','k','MarkerSize',8);
            %plot(data(n).ratio(data(n).veIndices).*data(n).meanDist, exp(data(n).x0s(data(n).veIndices)), strcat(colors{n}, 'o'), 'MarkerSize',3, 'MarkerFaceColor', colors{n})
            hold on; plot(xpoints, dataAll(n,1).ypoints, 'Color', colors{n},'LineWidth',2);
            hold on; plot(xpoints, dataAll(n,1).bupper, 'Color', colors{n});
            hold on; plot(xpoints, dataAll(n,1).blower, 'Color', colors{n});
            %hold on; plot(xpoints, data(n).ypoints, colors{n});
        end
        axis([0 1 0 10]);
        axis square
        
        %Summary Plot
        figure;
        hold on;
        for n=1:size(dataAll,1)
            xdata=dataAll(n,1).x{4};
            xdata=xdata-min(xdata(~isnan(dataAll(n,1).y{4})));
            xdata=xdata./max(xdata(~isnan(dataAll(n,1).y{4})));
            xdata=xdata.*0.9;
            xdata=xdata+0.05;
            plot(xdata,dataAll(n,1).y{4}, 'Color', colors{n},'LineWidth',2);
        end
        axis([0 1 0 10]);
        axis square
        
        %All ROIs num vs dist plot, normalized for range
        figure;
        hold on;
        highestX=0;
        lowestX=1;
        for n=1:size(dataAll,1)
            errorbar(dataAll(n,1).x{4}/(1/dataAll(n,1).logLineFit(1))+dataAll(n,1).logLineFit(2),dataAll(n,1).y{4},dataAll(n,1).ysterr{4}(1,:),dataAll(n,1).ysterr{4}(1,:),strcat('o'),...
                'Color', colors{n}, 'MarkerFaceColor',colors{n},'MarkerEdgeColor','k','MarkerSize',8);
             hold on; plot(dataAll(n,1).xpoints/(1/dataAll(n,1).logLineFit(1))+dataAll(n,1).logLineFit(2), dataAll(n,1).ypoints, 'Color',colors{n},'LineWidth',2);
             hold on; plot(dataAll(n,1).xpoints/(1/dataAll(n,1).logLineFit(1))+dataAll(n,1).logLineFit(2), dataAll(n,1).bupper, 'Color',colors{n});
             hold on; plot(dataAll(n,1).xpoints/(1/dataAll(n,1).logLineFit(1))+dataAll(n,1).logLineFit(2), dataAll(n,1).blower, 'Color',colors{n});
             if max(dataAll(n,1).x{4}(~isnan(dataAll(n,1).y{4}))/(1/dataAll(n,1).logLineFit(1))+dataAll(n,1).logLineFit(2))>highestX,
                 highestX=max(dataAll(n,1).x{4}(~isnan(dataAll(n,1).y{4}))/(1/dataAll(n,1).logLineFit(1))+dataAll(n,1).logLineFit(2));
             end
             if min(dataAll(n,1).x{4}(~isnan(dataAll(n,1).y{4}))/(1/dataAll(n,1).logLineFit(1))+dataAll(n,1).logLineFit(2))<lowestX,
                 lowestX=min(dataAll(n,1).x{4}(~isnan(dataAll(n,1).y{4}))/(1/dataAll(n,1).logLineFit(1))+dataAll(n,1).logLineFit(2));
             end
        end
        axis([lowestX-(highestX-lowestX)/(binSteps*2) highestX+(highestX-lowestX)/(binSteps*2) 0 10]);
        axis square
        
                %All ROIs num vs dist plot, normalized for range, summary
        figure;
        hold on;
%         highestX=0;
%         lowestX=1;
        for n=1:size(dataAll,1)
            plot(dataAll(n,1).x{4}/(1/dataAll(n,1).logLineFit(1))+dataAll(n,1).logLineFit(2),dataAll(n,1).y{4},'Color',colors{n},'LineWidth',2);

%              if max(xpoints/(1/dataAll(n,1).logLineFit(1))+dataAll(n,1).logLineFit(2))>highestX,
%                  highestX=max(xpoints/(1/dataAll(n,1).logLineFit(1))+dataAll(n,1).logLineFit(2));
%              end
%              if min(xpoints/(1/dataAll(n,1).logLineFit(1))+dataAll(n,1).logLineFit(2))<lowestX,
%                  lowestX=min(xpoints/(1/dataAll(n,1).logLineFit(1))+dataAll(n,1).logLineFit(2));
%              end
        end
        axis([lowestX-(highestX-lowestX)/(binSteps*2) highestX+(highestX-lowestX)/(binSteps*2) 0 10]);
        axis square
end
end


function [ypoints, yfitParams, b_upper, b_lower,B]=bootstrapLogLineFitter(x,y,ve, xpoints)
x = x(isfinite(ve)); y = y(isfinite(ve)); ve = ve(isfinite(ve));
ve=ones(size(ve));
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
end

function [B,e]=logLineFitter(iis, x,y,ve)
%x = x(isfinite(ve)); y = y(isfinite(ve)); ve = ve(isfinite(ve));
x=x(iis);
y=y(iis);
ve=ve(iis);
options = optimset('MaxFunEvals',10000000);
[B,e] = fminsearch(@(z) mylogfit(z,x,y,ve),[0.2;1.3], options);
end

function e=mylogfit(z,x,y,ve)
e=sum(ve.*(y-(exp(z(1).*x+z(2)))).^2)./sum(ve);
end

function [ypoints, yfitParams, b_upper, b_lower]=bootstrapLineFitter(x,y,ve, xpoints)
x = x(isfinite(ve)); y = y(isfinite(ve)); ve = ve(isfinite(ve));
ve=ones(size(ve));
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
end

function [B,e]=lineFitter(iis, x,y,ve)
%x = x(isfinite(ve)); y = y(isfinite(ve)); ve = ve(isfinite(ve));
x=x(iis);
y=y(iis);
ve=ve(iis);
options = optimset('MaxFunEvals',10000000);
[B,e] = fminsearch(@(z) mylinfit(z,x,y,ve),[0.2;1.3], options);
end

function e=mylinfit(z,x,y,ve)
e=sum(ve.*(y-((z(1).*x+z(2)))).^2)./sum(ve);
end


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
end


function [B, e]=cmfLineFitter(ii,x,y,ve)
x = x(ii); y = y(ii); ve = ve(ii);
options = optimset('MaxFunEvals',10000000);
[B, e] = fminsearch(@(z) mycmffit(z,x,y,ve),[0.05;0.2], options);
end

function e=mycmffit(z,x,y,ve)
e=sum(ve.*(y-(1./(z(1).*x+z(2)))).^2)./sum(ve);
end