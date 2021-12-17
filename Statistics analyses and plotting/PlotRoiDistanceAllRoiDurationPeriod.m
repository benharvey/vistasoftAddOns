function [stats] = PlotRoiDistanceAllRoiDurationPeriod(data, veThresh, showPlots, binSteps, meanThresh)
%Makes various plots from the output of RoiDistanceRatio. THis version is
%made to use Log2Lin models. It's a simpler version of
%PlotRoiDistanceDataAllBin, giving fewer and simpler outputs

if ~exist('veThresh', 'var') || isempty(veThresh)
    veThresh=0;
end
if ~exist('showPlots', 'var') || isempty(showPlots)
    showPlots=[1 1 1 1 1 1];
end
if ~exist('binSteps', 'var') || isempty(binSteps)
    binSteps=11;
end
if ~exist('meanThresh', 'var') || isempty(meanThresh)
    meanThresh=zeros(size(data,1), 1);
end

colors={'k', 'r', 'g', 'b', 'm', 'c'};

dataAll=data;
clear data;

fig1=figure; 
hold on;
for whichRoi=1:size(dataAll, 1)
    if ~isempty(dataAll{whichRoi,1})
        data{1}=dataAll{whichRoi,1};
        data{2}=dataAll{whichRoi,2};
        data{3}=dataAll{whichRoi,3};
        
        for n=1:length(data)
            %         %Just run this part the first time through
            %
            %         data{n}.x0sLog=data{n}.x0s;
            %         data{n}.sigmasLog=data{n}.sigmas;
            %         fwhms=data{n}.sigmas.*(2*sqrt(2*log(2)));
            %         data{n}.sigmas=exp(data{n}.x0s+fwhms./2)-exp(data{n}.x0s-fwhms./2);
            %         data{n}.x0s=exp(data{n}.x0s);
            %
            %n=1;
            data{n}.veIndices=data{n}.ves>=veThresh & data{n}.x0s>0.05 & data{n}.x0s<1 & data{n}.y0s>0.05 & data{n}.y0s<1 & data{1}.meanSignal>=meanThresh(whichRoi);
            dataAll{whichRoi, n}.veIndices=data{n}.ves>=veThresh & data{n}.x0s>0.05 & data{n}.x0s<1 & data{n}.y0s>0.05 & data{n}.y0s<1 & data{1}.meanSignal>=meanThresh(whichRoi);
            %         %data{n}.x0fit=linreg(data{n}.ratio(data{n}.veIndices), data{n}.x0s(data{n}.veIndices), data{n}.ves(data{n}.veIndices));
            %         %data{n}.sigmafit=linreg(data{n}.ratio(data{n}.veIndices), data{n}.sigmas(data{n}.veIndices), data{n}.ves(data{n}.veIndices));
        end

        b=linspace(0.5/binSteps, 1-(0.5/binSteps), binSteps);
        
        %number vs mean distance
        
        %xpoints=0:0.01:data{1}.meanDist;
        %if showPlots(1)==1 || showPlots(2)==1
        bins=(1:2:ceil(data{n}.meanDist))./data{n}.meanDist;
        for n=1:length(data)
            if sum(data{n}.veIndices)>1
                tmp=corrcoef(data{n}.ratio(data{n}.veIndices),data{n}.x0s(data{n}.veIndices));
                rDistDurCorr(whichRoi, n)=tmp(1,2);
                tmp=corrcoef(data{n}.ratio(data{n}.veIndices),data{n}.y0s(data{n}.veIndices));
                rDistPerCorr(whichRoi, n)=tmp(1,2);
                ns(whichRoi, n)=floor(sum(data{n}.veIndices)./1.7^2);
                pDistDurCorr(whichRoi, n)=r2p(rDistDurCorr(whichRoi, n), ns(whichRoi, n));
                pDistPerCorr(whichRoi, n)=r2p(rDistPerCorr(whichRoi, n), ns(whichRoi, n));
                
                tmp=corrcoef(data{n}.x0s(data{n}.veIndices), data{n}.y0s(data{n}.veIndices));
                rDurPerCorr(whichRoi,n)=tmp(1,2);
                pDurPerCorr(whichRoi, n)=r2p(rDurPerCorr(whichRoi, n), ns(whichRoi, n));
            end
            %     %Quick and dirty way, evenly-spaced bins at odd intervals
            %     data{n}.x{4}=b.*data{n}.meanDist;
            %     data{n}.y{4}=data{n}.y{3};
            %     data{n}.ysterr{4}=data{n}.ysterr{3};
            %     [data{n}.ypoints, data{n}.logLineFit, data{n}.bupper, data{n}.blower]=bootstrapLogLineFitter(data{n}.x{4},data{n}.y{4},1./data{n}.ysterr{4}(1,:), xpoints);
            
            %   Slower, but with bins at regular distance intervals
            if showPlots(1)==1 || showPlots(2)==1
                data{n}.y{4}=nan(size(bins));
                data{n}.ysterr{4}=nan(2,length(bins));
                data{n}.y{5}=nan(size(bins));
                data{n}.ysterr{5}=nan(2,length(bins));
                for bCount=1:length(bins)
                    bii = data{n}.ratio> bins(bCount)-0.5/binSteps & ...
                        data{n}.ratio< bins(bCount)+0.5/binSteps & data{n}.veIndices;
                    if any(bii)
                        s=wstat(data{n}.x0s(bii), data{n}.ves(bii), 1.77^2);
                        data{n}.y{4}(bCount)=s.mean;
                        data{n}.ysterr{4}(:,bCount)=s.sterr;
                        s=wstat(data{n}.y0s(bii), data{n}.ves(bii), 1.77^2);
                        data{n}.y{5}(bCount)=s.mean;
                        data{n}.ysterr{5}(:,bCount)=s.sterr;
                    end
                end
                data{n}.x{4}=bins;
                data{n}.x{5}=bins;
                xpoints=min(bins(~isnan(data{n}.y{4}))):0.01:max(bins(~isnan(data{n}.y{4})));
                data{n}.xpoints=xpoints;
                try
                    [data{n}.ypointsX, data{n}.logLineFitX, data{n}.bupperX, data{n}.blowerX, data{n}.B]=bootstrapLogLineFitter(data{n}.x{4},data{n}.y{4},1./data{n}.ysterr{4}(1,:), xpoints);
                    [data{n}.ypointsY, data{n}.logLineFitY, data{n}.bupperY, data{n}.blowerY, data{n}.B]=bootstrapLogLineFitter(data{n}.x{5},data{n}.y{5},1./data{n}.ysterr{5}(1,:), xpoints);
                    
                catch
                    n
                end
                %Statistical test of progression with distance (superceded
                %by correlations above)
%                 if showPlots(1)==1
%                     permutations=10000;
%                     ydat=log(data{n}.y{4}(~isnan(data{n}.y{4})));
%                     xdat=data{n}.x{4}(~isnan(data{n}.y{4}));
%                     fitDist=zeros(2,permutations);
%                     for m=1:permutations
%                         yshuffle=ydat(randperm(length(ydat)));
%                         fitDist(:,m)=linreg(xdat, yshuffle);
%                     end
%                     CI=[prctile(fitDist(2,:), 2.5) prctile(fitDist(2,:), 97.5)];
%                     max(fitDist(2,:));
%                     measure=linreg(xdat, ydat);
%                     measure(2);
%                     [whichRoi, n]
%                     if measure(2)<max(fitDist(2,:))
%                         tmp= fitDist(2,:);
%                         tmp=sort(tmp);
%                         tmp=tmp>measure(2);
%                         pXprog(whichRoi, n)=(permutations-find(tmp, 1, 'first')+1)./permutations;
%                         nYprog=sum(~isnan(data{n}.y{4}));
%                     else
%                         pXprog(whichRoi, n)=0;
%                         nXprog=sum(~isnan(data{n}.y{4}));
%                     end
%                     
%                     ydat=log(data{n}.y{5}(~isnan(data{n}.y{5})));
%                     xdat=data{n}.x{5}(~isnan(data{n}.y{5}));
%                     fitDist=zeros(2,permutations);
%                     for m=1:permutations
%                         yshuffle=ydat(randperm(length(ydat)));
%                         fitDist(:,m)=linreg(xdat, yshuffle);
%                     end
%                     CI=[prctile(fitDist(2,:), 2.5) prctile(fitDist(2,:), 97.5)];
%                     max(fitDist(2,:));
%                     measure=linreg(xdat, ydat);
%                     measure(2);
%                     [whichRoi, n]
%                     if measure(2)<max(fitDist(2,:))
%                         tmp= fitDist(2,:);
%                         tmp=sort(tmp);
%                         tmp=tmp>measure(2);
%                         pYprog(whichRoi, n)=(permutations-find(tmp, 1, 'first')+1)./permutations;
%                         nYprog=sum(~isnan(data{n}.y{5}));
%                     else
%                         pYprog(whichRoi, n)=0;
%                         nYprog=sum(~isnan(data{n}.y{5}));
%                     end
%                 end

                %To compare log and linear fits
                %         if n==1
                %             [tmp1, data{n}.linLineFit, tmp2,tmp3]=bootstrapLineFitter(data{n}.x{4},data{n}.y{4},1./data{n}.ysterr{4}(1,:), xpoints);
                %             data{n}.rawrss=var(data{n}.y{4});
                %             data{n}.rsslin=var(data{n}.y{4}-(data{n}.x{4}.*data{n}.linLineFit(1)+data{n}.linLineFit(2)));
                %             data{n}.rsslog=var(data{n}.y{4}-exp(data{n}.x{4}.*data{n}.logLineFit(1)+data{n}.logLineFit(2)));
                %             data{n}.linresid=(data{n}.y{4}-(data{n}.x{4}.*data{n}.linLineFit(1)+data{n}.linLineFit(2)));
                %             data{n}.logresid=(data{n}.y{4}-exp(data{n}.x{4}.*data{n}.logLineFit(1)+data{n}.logLineFit(2)));
                %         end
            end
        end
        
        %     if showPlots(1)==1
        %         figure;
        %         hold on;
        %         for n=1:length(data)
        %             errorbar(data{n}.x{4}.*data{n}.meanDist,data{n}.y{4},data{n}.ysterr{4}(1,:),data{n}.ysterr{4}(1,:),strcat(colors{n},'o'),...
        %                 'MarkerFaceColor',colors{n},'MarkerEdgeColor','k','MarkerSize',8);
        %             %plot(data{n}.ratio(data{n}.veIndices).*data{n}.meanDist, exp(data{n}.x0s(data{n}.veIndices)), strcat(colors{n}, 'o'), 'MarkerSize',3, 'MarkerFaceColor', colors{n})
        %
        %             hold on; plot(data{n}.xpoints.*data{n}.meanDist, data{n}.ypoints, colors{n},'LineWidth',2);
        %             hold on; plot(data{n}.xpoints.*data{n}.meanDist, data{n}.bupper, colors{n});
        %             hold on; plot(data{n}.xpoints.*data{n}.meanDist, data{n}.blower, colors{n});
        %
        %             %hold on; plot(xpoints, data{n}.ypoints, colors{n});
        %         end
        %         axis([0 5*ceil(data{n}.meanDist/5) 0 1]);
        %         axis square
        %     end
        
        %Summary plot, linked points, no errorbars
        if showPlots(2)==1
           % try
                plotNum=floor((whichRoi+1)/2);
                if mod(whichRoi, 2)==0
                    plotNum=plotNum+size(dataAll, 1)/2;
                end
                ax{whichRoi}=subplot(2, size(dataAll, 1)/2, plotNum);

            for n=[3 2 1]%1:length(data)
                try
                if n==1

                    %plot(data{n}.ratio(data{n}.veIndices).*data{n}.meanDist, exp(data{n}.x0s(data{n}.veIndices)), strcat(colors{n}, 'o'), 'MarkerSize',3, 'MarkerFaceColor', colors{n})
                    hold on; plot(data{n}.xpoints.*data{n}.meanDist, data{n}.ypointsY, 'Color', [0.5 0.5 0.5],'LineWidth',2);
                    hold on; plot(data{n}.xpoints.*data{n}.meanDist, data{n}.bupperY, '--', 'Color', [0.5 0.5 0.5]);
                    hold on; plot(data{n}.xpoints.*data{n}.meanDist, data{n}.blowerY, '--', 'Color', [0.5 0.5 0.5]);
                    
                    hold on; plot(data{n}.xpoints.*data{n}.meanDist, data{n}.ypointsX, colors{n},'LineWidth',2);
                    hold on; plot(data{n}.xpoints.*data{n}.meanDist, data{n}.bupperX, [colors{n}, '--'],'LineWidth',0.5);
                    hold on; plot(data{n}.xpoints.*data{n}.meanDist, data{n}.blowerX, [colors{n}, '--'],'LineWidth',0.5);
                    
                    hold on; errorbar(data{n}.x{5}.*data{n}.meanDist,data{n}.y{5},data{n}.ysterr{5}(1,:),data{n}.ysterr{5}(1,:),strcat(colors{n},'.'),...
                        'MarkerFaceColor',[0.5 0.5 0.5],'LineWidth',0.5, 'Color', [0.5 0.5 0.5]);
                    %plot(data{n}.ratio(data{n}.veIndices).*data{n}.meanDist, exp(data{n}.x0s(data{n}.veIndices)), strcat(colors{n}, 'o'), 'MarkerSize',3, 'MarkerFaceColor', colors{n})
                    hold on; errorbar(data{n}.x{4}.*data{n}.meanDist,data{n}.y{4},data{n}.ysterr{4}(1,:),data{n}.ysterr{4}(1,:),strcat(colors{n},'.'),...
                        'MarkerFaceColor',[0 0 0],'LineWidth',0.5);
                    
                    hold on; plot(data{n}.x{5}.*data{n}.meanDist,data{n}.y{5},strcat(colors{n},'o'),...
                        'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','k','MarkerSize',3);
                    hold on; plot(data{n}.x{4}.*data{n}.meanDist,data{n}.y{4},strcat(colors{n},'o'),...
                        'MarkerFaceColor',[0 0 0],'MarkerEdgeColor','k','MarkerSize',3);
                else
                    hold on; plot(data{n}.x{4}.*data{n}.meanDist,data{n}.y{4}, colors{n},'LineWidth',1);
                    hold on; plot(data{n}.x{5}.*data{n}.meanDist,data{n}.y{5}, colors{n+3},'LineWidth',1);
                end
                end
            end
            axis square;
            axis([0 5*ceil(max([max(data{n}.x{4}.*data{n}.meanDist) max(data{n}.x{5}.*data{n}.meanDist)])/5) 0 1]);
%            xticks(0:5:(5*ceil(max([max(data{n}.x{4}.*data{n}.meanDist) max(data{n}.x{5}.*data{n}.meanDist)])/5)))
            axis square;
            drawnow;
            %saveas(gcf, ['DurationDistance' data{n}.ROItitle], 'epsc');
            %catch
            %end
        end

           
        if showPlots(3)==1;
            figure; plot(data{1}.x0s(data{1}.veIndices), data{1}.y0s(data{1}.veIndices), 'k.');
            line=linreg(data{1}.x0s(data{1}.veIndices), data{1}.y0s(data{1}.veIndices))
            hold on; plot([0 1], [line(1) line(1)+line(2)], 'k');
            axis([0 1 0 1])
            axis square;
            title(['PeriodDuration' data{n}.ROItitle]);
            saveas(gcf, ['PeriodDuration' data{n}.ROItitle], 'epsc');
        end
        
        try
        indices=data{2}.veIndices & data{3}.veIndices;
        tmp=corrcoef(data{2}.x0s(indices), data{3}.x0s(indices));
        rDurOddEven(whichRoi)=tmp(1,2);
        tmp=corrcoef(data{2}.y0s(indices), data{3}.y0s(indices));
        rPerOddEven(whichRoi)=tmp(1,2);
        nsOddEven(whichRoi)=length(data{2}.x0s(indices))./1.7^2;
        pDurOddEven(whichRoi)=r2p(rDurOddEven(whichRoi), nsOddEven(whichRoi));
        pPerOddEven(whichRoi)=r2p(rPerOddEven(whichRoi), nsOddEven(whichRoi));
        
        if showPlots(4)==1;
            figure; plot(data{2}.x0s(indices), data{3}.x0s(indices), 'k.');
            hold on; plot(data{2}.y0s(indices), data{3}.y0s(indices), 'b.');
            line=linreg(data{2}.x0s(indices), data{3}.x0s(indices))
            hold on; plot([0 1], [line(1) line(1)+line(2)], 'k');
            line=linreg(data{2}.y0s(indices), data{3}.y0s(indices))
            hold on; plot([0 1], [line(1) line(1)+line(2)], 'b');
            axis([0 1 0 1])
            axis square;
            title(['OddEvenCor' data{n}.ROItitle]);
            saveas(gcf, ['OddEvenCor' data{n}.ROItitle], 'epsc');
        end
        
        if showPlots(5)==1;
            [X,Y] = meshgrid(0:0.001:1,0:0.001:1);
            veIndices=find(data{1}.veIndices);
            rf=zeros([length(X(:)), length(veIndices)]) ;
            for whichRF=1:length(veIndices)
                rftmp   = rfGaussian2d(X(:), Y(:),...
                    data{1}.sigmas(veIndices(whichRF)), ...
                    data{1}.sigmaMinor(veIndices(whichRF)), ...
                    data{1}.sigmaTheta(veIndices(whichRF)), ...
                    data{1}.x0s(veIndices(whichRF)), ...
                    data{1}.y0s(veIndices(whichRF)));
                rf(:,whichRF)=rftmp;
            end
            
            vol = data{1}.sigmas(data{1}.veIndices).*data{1}.sigmaMinor(data{1}.veIndices).*data{1}.ves(data{1}.veIndices);
            vol = vol * (2 * pi);
            rf=rf./(ones(size(rf)) * vol');
            [~, rfPeaks]=max(rf, [], 1);
            rfall=max(rf, [], 2);
            rfall(rfPeaks)=min(rfall(:));
            rfall=reshape(rfall, size(X));
            figure; imagesc(flipud(rfall)); colormap hot; axis square
            saveas(gcf, ['Coverage' data{n}.ROItitle], 'epsc');
        end
        catch 
        end
        dataOut{whichRoi}=data;
    end
    
end

%False discovery rate correction

for whichDT=1:3
    [~,~,~,pDistDurCorr(1:2:size(pDistDurCorr, 1),whichDT)]=fdr_bh(pDistDurCorr(1:2:size(pDistDurCorr, 1),whichDT));
    [~,~,~,pDistDurCorr(2:2:size(pDistDurCorr, 1),whichDT)]=fdr_bh(pDistDurCorr(2:2:size(pDistDurCorr, 1),whichDT));
    [~,~,~,pDistPerCorr(1:2:size(pDistPerCorr, 1),whichDT)]=fdr_bh(pDistPerCorr(1:2:size(pDistPerCorr, 1),whichDT));
    [~,~,~,pDistPerCorr(2:2:size(pDistPerCorr, 1),whichDT)]=fdr_bh(pDistPerCorr(2:2:size(pDistPerCorr, 1),whichDT));
end
 %Generating significance stars for the title
 if showPlots(2)==1
     for whichRoi=1:20
         for whichSig=1:3
             if whichSig==1
                 sigDur{whichSig}='\color{black}';
             elseif whichSig==2
                 sigDur{whichSig}='\color{red}';
             elseif whichSig==3
                 sigDur{whichSig}='\color{green}';
             end
             
             if rDistDurCorr(whichRoi, whichSig)>0
                 if pDistDurCorr(whichRoi, whichSig)<0.0001
                     sigDur{whichSig}=[sigDur{whichSig} '****'];
                 elseif pDistDurCorr(whichRoi, whichSig)<0.001
                     sigDur{whichSig}=[sigDur{whichSig} '***'];
                 elseif pDistDurCorr(whichRoi, whichSig)<0.01
                     sigDur{whichSig}=[sigDur{whichSig} '**'];
                 elseif pDistDurCorr(whichRoi, whichSig)<0.05
                     sigDur{whichSig}=[sigDur{whichSig} '*'];
                 end
             end
             if whichSig==1
                 sigPer{whichSig}='\color{gray}';
             elseif whichSig==2
                 sigPer{whichSig}='\color{magenta}';
             elseif whichSig==3
                 sigPer{whichSig}='\color{cyan}';
             end
             if rDistPerCorr(whichRoi, whichSig)>0
                 if pDistPerCorr(whichRoi, whichSig)<0.0001
                     sigPer{whichSig}=[sigPer{whichSig} '****'];
                 elseif pDistPerCorr(whichRoi, whichSig)<0.001
                     sigPer{whichSig}=[sigPer{whichSig} '***'];
                 elseif pDistPerCorr(whichRoi, whichSig)<0.01
                     sigPer{whichSig}=[sigPer{whichSig} '**'];
                 elseif pDistPerCorr(whichRoi, whichSig)<0.05
                     sigPer{whichSig}=[sigPer{whichSig} '*'];
                 end
             end
         end
         plotNum=floor((whichRoi+1)/2);
         if mod(whichRoi, 2)==0
             plotNum=plotNum+size(dataAll, 1)/2;
         end
         subplot(2, size(dataAll, 1)/2, plotNum)
         hold on;
         t=title([sigDur sigPer]);
         set(t, 'horizontalAlignment', 'left')
         drawnow;
     end
     
     %Moving subplots around
     for whichRoi=[2:20 1]
         if ~isempty(ax{whichRoi})
             plotNum=floor((whichRoi+1)/2);
             if mod(whichRoi, 2)==0
                 plotNum=plotNum+size(dataAll, 1)/2;
             end
             %ax{whichRoi}=subplot(2, 10, plotNum);
             %Specific to 10 maps
             xpos=plotNum/10;
             ypos=plotNum>10;
             if ypos==0
                 ypos=0.7;
             else
                 ypos=0.2;
             end
             if xpos>1; xpos=xpos-1; end
             set(ax{whichRoi}, 'position', [xpos-0.135 ypos 0.17 0.17]);
         end
     end
 end

if showPlots(1)==1
    stats.pXprog=pXprog;
    stats.pYprog=pYprog;
    stats.nXprog=nXprog;
    stats.nYprog=nYprog;  
end

stats.rDistDurCorr=rDistDurCorr;
stats.rDistPerCorr=rDistPerCorr;
stats.ns=ns;
stats.pDistDurCorr=pDistDurCorr;
stats.pDistPerCorr=pDistPerCorr;
stats.rDurPerCorr=rDurPerCorr;
stats.pDurPerCorr=pDurPerCorr;

stats.rDurOddEven=rDurOddEven;
stats.rPerOddEven=rPerOddEven;
stats.nsOddEven=nsOddEven;
stats.pDurOddEven=pDurOddEven;
stats.pPerOddEven=pPerOddEven;

%All left then all right
%colors={[0 0 255]./255, [128 0 255]./255, [255 0 255]./255, [255 0 0]./255, [255 128 0]./255, [128 64 0]./255, [128 128 128]./255, [0 128 0]./255, ...
%   [0 0 255]./255, [128 0 255]./255, [255 0 255]./255, [255 0 0]./255, [255 128 0]./255, [128 64 0]./255, [128 128 128]./255, [0 128 0]./255 };
%Alternating hemispheres
colors={[0 0 255]./255, [0 0 255]./255, [128 0 255]./255, [128 0 255]./255, [255 0 255]./255, [255 0 255]./255, [255 0 0]./255, [255 0 0]./255, [255 128 0]./255, [255 128 0]./255, [128 64 0]./255, [128 64 0]./255, [128 128 128]./255, [128 128 128]./255, [0 128 0]./255, [0 128 0]./255 , [0 255 0]./255, [0 255 0]./255};
colors={[0.625 0.625 0.625], [0.375 0.375 0.375],[0.625 0.625 0.625], [0.375 0.375 0.375],[0.625 0.625 0.625], [0.375 0.375 0.375],[0.625 0.625 0.625], [0.375 0.375 0.375],[0.625 0.625 0.625], [0.375 0.375 0.375],[0.625 0.625 0.625], [0.375 0.375 0.375],[0.625 0.625 0.625], [0.375 0.375 0.375],[0.625 0.625 0.625], [0.375 0.375 0.375],[0.625 0.625 0.625], [0.375 0.375 0.375],[0.625 0.625 0.625], [0.375 0.375 0.375],[0.625 0.625 0.625], [0.375 0.375 0.375],[0.625 0.625 0.625], [0.375 0.375 0.375]};

if showPlots(6)==1;
    for n=1:size(dataAll,1)
        
            if ~isempty(dataAll{n,1})
                [~, ~, T, P, df]=lineFitterV(1:sum(dataAll{n,1}.veIndices), dataAll{n,1}.x0s(dataAll{n,1}.veIndices),dataAll{n,1}.sigmas(dataAll{n,1}.veIndices),ones(size(dataAll{n,1}.sigmas(dataAll{n,1}.veIndices))), 1.7.^2);
                tTuningWidthMajor(n,:)=T;
                pTuningWidthMajor(n,:)=P;
                nTuningWidth(n)=df+3;
                nUnderOverTuningWidth(n,:)=[sum(dataAll{n,1}.x0s(dataAll{n,1}.veIndices)<0.5)./1.7^2 sum(dataAll{n,1}.x0s(dataAll{n,1}.veIndices)>0.5)./1.7^2];
                [~, ~, T, P, df]=lineFitterV(1:sum(dataAll{n,1}.veIndices), dataAll{n,1}.x0s(dataAll{n,1}.veIndices),dataAll{n,1}.sigmaMinor(dataAll{n,1}.veIndices),ones(size(dataAll{n,1}.sigmas(dataAll{n,1}.veIndices))), 1.7.^2);
                tTuningWidthMinor(n,:)=T;
                pTuningWidthMinor(n,:)=P;
                
                if showPlots(6)==1
                    
                    
                    b=linspace(0, 1, 21);
                    dataAll{n,1}.x{8}=b;
                    dataAll{n,1}.y{8}=nan(size(b));
                    dataAll{n,1}.ysterr{8}=nan(2,length(b));
                    dataAll{n,1}.x{9}=b;
                    dataAll{n,1}.y{9}=nan(size(b));
                    dataAll{n,1}.ysterr{9}=nan(2,length(b));
                    %         dataAll{n,1}.yLog{9}=nan(size(b));
                    %         dataAll{n,1}.ysterrLog{9}=nan(2,length(b));
                    for bCount=1:length(b)
                        bii = dataAll{n,1}.x0s> b(bCount)-(b(2)-b(1))./2 & ...
                            dataAll{n,1}.x0s< b(bCount)+(b(2)-b(1))./2 & ...
                            dataAll{n,1}.sigmas<=20 & dataAll{n,1}.veIndices;
                        if any(bii) && sum(dataAll{n,1}.ves(bii))>0 && numel(dataAll{n,1}.sigmas(bii))>1.77^2;
                            s=wstat(dataAll{n,1}.sigmas(bii), dataAll{n,1}.ves(bii), 1.77^2);
                            dataAll{n,1}.y{8}(bCount)=s.mean;
                            dataAll{n,1}.ysterr{8}(:,bCount)=s.sterr;
                            s=wstat(dataAll{n,1}.sigmaMinor(bii), dataAll{n,1}.ves(bii), 1.77^2);
                            dataAll{n,1}.y{9}(bCount)=s.mean;
                            dataAll{n,1}.ysterr{9}(:,bCount)=s.sterr;
                            
                        end
                    end
                    dataAll{n,1}.xpointsSigma=min(b(isfinite(dataAll{n,1}.y{9}))):0.01:max(b(isfinite(dataAll{n,1}.y{9})));%max(exp(data(n).x0s));
                    
                    
                    %data(n).fwhmnumfit=linreg(dataAll{n,1}.x{9}(isfinite(1./dataAll{n,1}.y{sterr9}(1,:))), dataAll{n,1}.y{9}(isfinite(1./dataAll{n,1}.y{sterr9}(1,:))), 1./dataAll{n,1}.y{sterr9}(1,isfinite(1./dataAll{n,1}.y{sterr9}(1,:))));
                    [dataAll{n,1}.ypointsSigMajor, dataAll{n,1}.SigMajorNumFit, dataAll{n,1}.SigMajorbupper, dataAll{n,1}.SigMajorblower]=bootstrapLineFitterV(dataAll{n,1}.x{8},dataAll{n,1}.y{8},1./dataAll{n,1}.ysterr{8}(1,:), dataAll{n,1}.xpointsSigma);
                    [dataAll{n,1}.ypointsSigMinor, dataAll{n,1}.SigMinorNumFit, dataAll{n,1}.SigMinorbupper, dataAll{n,1}.SigMinorblower]=bootstrapLineFitterV(dataAll{n,1}.x{9},dataAll{n,1}.y{9},1./dataAll{n,1}.ysterr{9}(1,:), dataAll{n,1}.xpointsSigma);
                    
                    %Balanced bootstrap of individual recording sites
                    xpoints=dataAll{n,1}.x0s(dataAll{n,1}.veIndices);
                    sigmaMajors=dataAll{n,1}.sigmas(dataAll{n,1}.veIndices);
                    sigmaMinors=dataAll{n,1}.sigmaMinor(dataAll{n,1}.veIndices);
                    ves=dataAll{n,1}.ves(dataAll{n,1}.veIndices);
                    lowIndices=xpoints<0.5;
                    highIndices=xpoints>=0.5;
                    highList=find(highIndices);
                    lowList=find(lowIndices);
                    if sum(lowIndices)>3*1.75^2 && sum(highIndices)>3*1.75^2
                        for permutation=1:1000
                            if length(lowList)>length(highList)
                                lowListPerm=Shuffle(lowList);
                                lowListPerm=lowListPerm(1:length(highList));
                                fullList=[lowListPerm highList];
                            else
                                highListPerm=Shuffle(highList);
                                highListPerm=highListPerm(1:length(lowList));
                                fullList=[lowList highListPerm];
                            end
                            [B, e, Tmajor(permutation,:), P, df]=lineFitterV(1:length(fullList), xpoints(fullList),sigmaMajors(fullList),ones(size(fullList)), 1.75^2);
                            [B, e, Tminor(permutation,:), P, df]=lineFitterV(1:length(fullList), xpoints(fullList),sigmaMinors(fullList),ones(size(fullList)), 1.75^2);
                        end
                        TsigMajorMean(n,:)=mean(Tmajor(:,1:2),1);
                        TsigMinorMean(n,:)=mean(Tminor(:,1:2),1);
                        for whichP=1:2
                            PsigMajorMean(n, whichP) = 2*tpvalue(-abs(TsigMajorMean(n,whichP)),df);
                            PsigMinorMean(n, whichP) = 2*tpvalue(-abs(TsigMinorMean(n,whichP)),df);
                        end
                        Nsigma(n)=floor(length(fullList)./1.75^2);
                    end
                    
                    %Permutation test of bins
                    if showPlots(7)==1
                        permutations=10000;
                        ydat=dataAll{n,1}.y{8}(~isnan(dataAll{n,1}.y{8}));
                        xdat=b(~isnan(dataAll{n,1}.y{8}));
                        fitDist=zeros(3,permutations);
                        for m=1:permutations
                            yshuffle=ydat(randperm(length(ydat)));
                            fitDist(:,m)=lineFitterV(1:length(xdat), xdat,yshuffle,ones(size(xdat)));
                        end
                        CI=[prctile(fitDist(1:2,:)', 2.5); prctile(fitDist(1:2,:)', 97.5)];
                        measure=lineFitterV(1:length(xdat), xdat,ydat,ones(size(xdat)));
                        
                        if measure(1)<max(fitDist(1,:))
                            tmp= fitDist(1,:);
                            tmp=sort(tmp);
                            tmp=tmp>measure(1);
                            pSMajProg(n,1)=(permutations-find(tmp, 1, 'first')+1)./permutations;
                        else
                            pSMajProg(n,1)=0;
                        end
                        if measure(2)<max(fitDist(2,:))
                            tmp= fitDist(2,:);
                            tmp=sort(tmp);
                            tmp=tmp>measure(2);
                            pSMajProg(n,2)=(permutations-find(tmp, 1, 'first')+1)./permutations;
                        else
                            pSMajProg(n,2)=0;
                        end
                        
                        ydat=dataAll{n,1}.y{9}(~isnan(dataAll{n,1}.y{9}));
                        xdat=b(~isnan(dataAll{n,1}.y{9}));
                        fitDist=zeros(3,permutations);
                        for m=1:permutations
                            yshuffle=ydat(randperm(length(ydat)));
                            fitDist(:,m)=lineFitterV(1:length(xdat), xdat,yshuffle,ones(size(xdat)));
                        end
                        CI=[prctile(fitDist(1:2,:)', 2.5); prctile(fitDist(1:2,:)', 97.5)];
                        measure=lineFitterV(1:length(xdat), xdat,ydat,ones(size(xdat)));
                        
                        if measure(1)<max(fitDist(1,:))
                            tmp= fitDist(1,:);
                            tmp=sort(tmp);
                            tmp=tmp>measure(1);
                            pSMinProg(n,1)=(permutations-find(tmp, 1, 'first')+1)./permutations;
                        else
                            pSMinProg(n,1)=0;
                        end
                        if measure(2)<max(fitDist(2,:))
                            tmp= fitDist(2,:);
                            tmp=sort(tmp);
                            tmp=tmp>measure(2);
                            pSMinProg(n,2)=(permutations-find(tmp, 1, 'first')+1)./permutations;
                        else
                            pSMinProg(n,2)=0;
                        end
                    end
                    
                    
                    %Linear fit only
                    %         if showPlots(7)==1
                    %             permutations=10000;
                    %             ydat=dataAll{n,1}.y{8}(~isnan(dataAll{n,1}.y{8}));
                    %             xdat=b(~isnan(dataAll{n,1}.y{8}));
                    %             fitDist=zeros(2,permutations);
                    %             for m=1:permutations
                    %                 yshuffle=ydat(randperm(length(ydat)));
                    %                 fitDist(:,m)=linreg(xdat, yshuffle);
                    %             end
                    %             CI=[prctile(fitDist(2,:), 2.5) prctile(fitDist(2,:), 97.5)];
                    %             max(fitDist(2,:));
                    %             measure=linreg(xdat, ydat);
                    %             measure(2);
                    %             if measure(2)<max(fitDist(2,:))
                    %                 tmp= fitDist(2,:);
                    %                 tmp=sort(tmp);
                    %                 tmp=tmp>measure(2);
                    %                 pSMajProg(n)=(permutations-find(tmp, 1, 'first')+1)./permutations;
                    %             else
                    %                 pSMajProg(n)=0;
                    %             end
                    %
                    %             ydat=dataAll{n,1}.y{9}(~isnan(dataAll{n,1}.y{9}));
                    %             xdat=b(~isnan(dataAll{n,1}.y{9}));
                    %             fitDist=zeros(2,permutations);
                    %             for m=1:permutations
                    %                 yshuffle=ydat(randperm(length(ydat)));
                    %                 fitDist(:,m)=linreg(xdat, yshuffle);
                    %             end
                    %             CI=[prctile(fitDist(2,:), 2.5) prctile(fitDist(2,:), 97.5)];
                    %             max(fitDist(2,:));
                    %             measure=linreg(xdat, ydat);
                    %             measure(2);
                    %             if measure(2)<max(fitDist(2,:))
                    %                 tmp= fitDist(2,:);
                    %                 tmp=sort(tmp);
                    %                 tmp=tmp>measure(2);
                    %                 pSMinProg(n)=(permutations-find(tmp, 1, 'first')+1)./permutations;
                    %             else
                    %                 pSMinProg(n)=0;
                    %             end
                    %         end
                end
            end

    end
    
    
    figure;
    hold on;
    for n=1:2:size(dataAll,1)%:length(data)
        plotNum=floor((n+1)/2);
        %For hemispheres on different lines
        %                 if mod(whichRoi, 2)==0
        %                     plotNum=plotNum+size(dataAll, 1)/2;
        %                 end
        ax{n}=subplot(1, size(dataAll, 1)/2, plotNum);
        
        if ~isempty(dataAll{n,1}) && showPlots(6)==1
            hold on; plot(dataAll{n,1}.xpointsSigma, dataAll{n,1}.SigMajorbupper, '--', 'Color',colors{n},'LineWidth',0.5);
            hold on; plot(dataAll{n,1}.xpointsSigma, dataAll{n,1}.SigMajorblower, '--', 'Color',colors{n},'LineWidth',0.5);
            hold on; plot(dataAll{n+1,1}.xpointsSigma, dataAll{n+1,1}.SigMajorbupper, '--', 'Color',colors{n+1},'LineWidth',0.5);
            hold on; plot(dataAll{n+1,1}.xpointsSigma, dataAll{n+1,1}.SigMajorblower, '--', 'Color',colors{n+1},'LineWidth',0.5);
            
            hold on; plot(dataAll{n,1}.xpointsSigma, dataAll{n,1}.SigMinorbupper, '--', 'Color',colors{n},'LineWidth',0.5);
            hold on; plot(dataAll{n,1}.xpointsSigma, dataAll{n,1}.SigMinorblower, '--', 'Color',colors{n},'LineWidth',0.5);
            hold on; plot(dataAll{n+1,1}.xpointsSigma, dataAll{n+1,1}.SigMinorbupper, '--', 'Color',colors{n+1},'LineWidth',0.5);
            hold on; plot(dataAll{n+1,1}.xpointsSigma, dataAll{n+1,1}.SigMinorblower, '--', 'Color',colors{n+1},'LineWidth',0.5);
            
            hold on; plot(dataAll{n,1}.xpointsSigma, dataAll{n,1}.ypointsSigMajor, 'Color', colors{n},'LineWidth',1);
            hold on; plot(dataAll{n+1,1}.xpointsSigma, dataAll{n+1,1}.ypointsSigMajor, 'Color', colors{n+1},'LineWidth',1);
            hold on; plot(dataAll{n,1}.xpointsSigma, dataAll{n,1}.ypointsSigMinor, 'Color', colors{n},'LineWidth',0.5);
            hold on; plot(dataAll{n+1,1}.xpointsSigma, dataAll{n+1,1}.ypointsSigMinor, 'Color', colors{n+1},'LineWidth',0.5);
            
            errorbar(b,dataAll{n,1}.y{8},dataAll{n,1}.ysterr{8}(1,:),dataAll{n,1}.ysterr{8}(1,:),strcat('k','o'),...
                'MarkerFaceColor',colors{n},'MarkerEdgeColor','k','MarkerSize',5, 'LineWidth',0.5, 'Color', colors{n});
            errorbar(b,dataAll{n+1,1}.y{8},dataAll{n+1,1}.ysterr{8}(1,:),dataAll{n+1,1}.ysterr{8}(1,:),strcat('k','o'),...
                'MarkerFaceColor',colors{n+1},'MarkerEdgeColor','k','MarkerSize',5, 'LineWidth',0.5, 'Color', colors{n+1});
                    
            errorbar(b,dataAll{n,1}.y{9},dataAll{n,1}.ysterr{9}(1,:),dataAll{n,1}.ysterr{9}(1,:),strcat('k','o'),...
                'MarkerFaceColor',colors{n},'MarkerEdgeColor','k','MarkerSize',3, 'LineWidth',0.5, 'Color', colors{n});
            errorbar(b,dataAll{n+1,1}.y{9},dataAll{n+1,1}.ysterr{9}(1,:),dataAll{n+1,1}.ysterr{9}(1,:),strcat('k','o'),...
                'MarkerFaceColor',colors{n+1},'MarkerEdgeColor','k','MarkerSize',3, 'LineWidth',0.5, 'Color', colors{n+1});
              
            axis([0 1 0 2]);
            axis square;
            drawnow;%
        end
    end
            
for whichDT=1:3
    [~,~,pTuningWidthMajor(1:2:size(pTuningWidthMajor, 1),whichDT)]=fdr_bh(pTuningWidthMajor(1:2:size(pTuningWidthMajor, 1),whichDT));
    [~,~,pTuningWidthMajor(2:2:size(pTuningWidthMajor, 1),whichDT)]=fdr_bh(pTuningWidthMajor(2:2:size(pTuningWidthMajor, 1),whichDT));
    [~,~,pTuningWidthMinor(1:2:size(pTuningWidthMinor, 1),whichDT)]=fdr_bh(pTuningWidthMinor(1:2:size(pTuningWidthMinor, 1),whichDT));
    [~,~,pTuningWidthMinor(2:2:size(pTuningWidthMinor, 1),whichDT)]=fdr_bh(pTuningWidthMinor(2:2:size(pTuningWidthMinor, 1),whichDT));
end
for whichDT=1:2
    [~,~,PsigMajorMean(1:2:size(PsigMajorMean, 1),whichDT)]=fdr_bh(PsigMajorMean(1:2:size(PsigMajorMean, 1),whichDT));
    [~,~,PsigMajorMean(2:2:size(PsigMajorMean, 1),whichDT)]=fdr_bh(PsigMajorMean(2:2:size(PsigMajorMean, 1),whichDT));
    [~,~,PsigMinorMean(1:2:size(PsigMinorMean, 1),whichDT)]=fdr_bh(PsigMinorMean(1:2:size(PsigMinorMean, 1),whichDT));
    [~,~,PsigMinorMean(2:2:size(PsigMinorMean, 1),whichDT)]=fdr_bh(PsigMinorMean(2:2:size(PsigMinorMean, 1),whichDT));
end
     
     %Moving subplots around
     for whichRoi=[2:10 1]
         plotNum=whichRoi;%floor((whichRoi+1)/2);
         
         ax{whichRoi}=subplot(2, size(dataAll, 1)/2, plotNum);
         %Specific to 10 maps
         xpos=plotNum/10;
         floor((whichRoi+1)/2)
         
         ypos=0.5;
         if xpos>1; xpos=xpos-1; end
         set(ax{floor((whichRoi*2-1))}, 'position', [xpos-0.135 ypos 0.17 0.17]);
     end

    
    stats.tTuningWidthMajor=tTuningWidthMajor;
    stats.pTuningWidthMajor=pTuningWidthMajor;
    stats.nTuningWidth=nTuningWidth;
    stats.nUnderOverTuningWidth=nUnderOverTuningWidth;
    stats.tTuningWidthMinor=tTuningWidthMinor;
    stats.pTuningWidthMinor=pTuningWidthMinor;
    
    stats.TsigMajorMean=TsigMajorMean;
    stats.TsigMinorMean=TsigMinorMean;
    stats.PsigMajorMean=PsigMajorMean;
    stats.PsigMinorMean=PsigMinorMean;
    stats.Nsigma=Nsigma;
end

if showPlots(7)==1
    stats.pSMajProg=pSMajProg;
    stats.pSMinProg=pSMinProg;
    
    for n=1:18;
        try
            tuningSlopesMajor(n,:)=dataAll{n,1}.SigMajorNumFit;
            tuningSlopesMinor(n,:)=dataAll{n,1}.SigMinorNumFit;
            tuningXs(n,:)=dataAll{n,1}.x{8};
            tuningYmajor(n,:)=dataAll{n,1}.y{8};
            tuningYminor(n,:)=dataAll{n,1}.y{9};
        end
    end
    stats.tuningSlopesMajor=tuningSlopesMajor;
    stats.tuningSlopesMinor=tuningSlopesMinor;
    stats.tuningXs=tuningXs;
    stats.tuningYmajor=tuningYmajor;
    stats.tuningYminor=tuningYminor;
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

function [ypoints, yfitParams, b_upper, b_lower, B]=bootstrapLineFitter(x,y,ve, xpoints)
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


function [ypoints, yfitParams, b_upper, b_lower, B]=bootstrapLineFitterV(x,y,ve, xpoints)
x = x(isfinite(ve)); y = y(isfinite(ve)); ve = ve(isfinite(ve));
ve=ones(size(ve));
iis=1:length(x);
[B] = bootstrp(1000,@(iis) lineFitterV(iis,x,y,ve),[1:length(x)]');
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
keep3 = B(3,:)>b_lower(3) &  B(3,:)<b_upper(3);
keep = keep1 & keep2 & keep3;

ypoints=yfitParams(3)+xpoints.*yfitParams(1)+abs(0.5-xpoints).*yfitParams(2);
fits = [xpoints' abs(0.5-xpoints') ones(size(xpoints'))]*B(:,keep);
b_upper = max(fits,[],2);
b_lower = min(fits,[],2);
end

function [B, e, T, P, df]=lineFitterV(iis, x,y,ve, upsample)
%x = x(isfinite(ve)); y = y(isfinite(ve)); ve = ve(isfinite(ve));
if ~exist('upsample', 'var') || isempty(upsample)
    upsample=1;
end
x=x(iis);
y=y(iis);
ve=ve(iis);
designMatrix=[x' abs(0.5-x)' ones(size(x))'];
B=pinv(designMatrix)*y';
pred=designMatrix*B;
pred=pred';

df=floor(length(iis)./upsample)-size(designMatrix, 2);
T=zeros(1, size(designMatrix,2));
for n=1:size(designMatrix, 2)
    c=zeros(1, size(designMatrix,2));
    c(n)=1;
    
    SE=sqrt((sum((y-pred).^2)./df)*(c*pinv(designMatrix'*designMatrix)*c'));
    T(n)=c*B./SE;
    P(n) = 2*tpvalue(-abs(T(n)),df);
end

e=sum(ve.*(y-pred).^2)./sum(ve);
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