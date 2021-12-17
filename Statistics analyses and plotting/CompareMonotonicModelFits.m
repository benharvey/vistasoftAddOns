%Plots bar chart and pairwise differences between candidate model fits.

%Determine where to find data in structure
for whichSub=1:length(subjectOrder)
    for whichMap=1:length(mapNames)
        for whichHemi=1:2;

        for modelN=1:length(modelNamesAll)
            targetDataOdd{whichSub, whichMap, whichHemi, modelN}=char(strcat(subjectOrder{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNamesAll{modelN}, '.',DTname, 'Odd.vesXval'));
            targetDataEven{whichSub, whichMap,whichHemi, modelN}=char(strcat(subjectOrder{whichSub}, '.', mapNames{whichMap}, '.',hemispheres{whichHemi},'.', modelNamesAll{modelN}, '.',DTname, 'Even.vesXval'));
        end
        end
    end
end

%Get data from these locations (cross validated variance explained)
barDataMeans=nan([length(subjectOrder), length(mapNames), 2,2,length(modelNamesAll)]);
for whichSub=1:length(subjectOrder)
    for whichMap=1:length(mapNames)
        if isfield(eval(subjectOrder{whichSub}), char(mapNames{whichMap}))
            %mapOK(whichSub,whichMap)=1;
            for whichHemi=1:2;
                if isfield(eval([char(subjectOrder{whichSub}), '.', char(mapNames{whichMap})]), char(hemispheres{whichHemi}))
                    for modelN=1:length(modelNamesAll);
                        barData{whichSub, whichMap, whichHemi, 1}(:,modelN)=eval(targetDataOdd{whichSub, whichMap, whichHemi, modelN});
                        barData{whichSub, whichMap, whichHemi, 2}(:,modelN)=eval(targetDataEven{whichSub, whichMap, whichHemi, modelN});
                    end
                    for oddEven=1:2
                        indices=max(barData{whichSub, whichMap, whichHemi, oddEven},[],2)>veThresh;
                        barDataMeans(whichSub, whichMap, whichHemi, oddEven,:)=mean(barData{whichSub, whichMap, whichHemi, oddEven}(indices,:), 1);
                    end
                end
            end
%         else
%             mapOK(whichSub,whichMap)=0;
        end
    end
end

%Compute and plot means and confidence intervals of model fits
barPoints=[];
whichBars=1:length(modelNamesAll);
for n=1:length(whichBars)
    tmp=barDataMeans(:, :, :, :,whichBars(n));
    tmp=(tmp(:));
    tmp=tmp(~isnan(tmp));
    barMeans(n)=mean(tmp);
    barStd(n)=std(tmp);
    barSerr(n)=std(tmp)/sqrt(length(tmp));
    barPoints(:,n)=tmp;
    CI95(n,:) = tinv([0.025 0.975], length(tmp)-1);
    CI95(n,:) =bsxfun(@times, barSerr(n), CI95(n,:));
end
figure; subplot(1,2,1);
bar(barMeans);
hold on; errorbar(1:size(barMeans,2), barMeans, CI95(:,1), CI95(:,2), 'k.')
axis square;
xlim([0.5 length(barMeans)+0.5])
ylim([0 0.35])

%Compute t-statistics and corresponding probabilities of pairwise
%differences between model fits
for x=1:length(whichBars)
    for y=1:length(whichBars)
        [~,p,ci,stats] = ttest(barPoints(:,x), barPoints(:,y), 'tail', 'both');
        pvals(x,y)=p;
        tMatrix(x,y)=stats.tstat;
    end
end
tMatrix(tMatrix>200)=200;
tMatrix(tMatrix<-200)=-200;

%This matrix uses 10*10 pixel cells, because some viewing software
%interpolates pixel edges
tMatrixImg=zeros(size(tMatrix,1)*10, size(tMatrix,2)*10);
for x=1:size(tMatrix,1)
    for y=1:size(tMatrix,1)
        tMatrixImg((x-1)*10+1:x*10, (y-1)*10+1:y*10)=tMatrix(x,y);
    end
end
tMatrixImg(isnan(tMatrixImg))=0;
%Add image of resulting t-statistics
figure; subplot(1,2,2);
imagesc(tMatrixImg);