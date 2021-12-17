function AFNIclusterCountVsThreshold(selSizes, selVarExps, paths, whichSubs)

%Plots cluster count as a function of variance explained, to find knee
%point of count (where count flattens out for a large range ov variance
%explained thresholds. Choose the threshold at the bottom of the flat
%range, by eye for now.

%01/2020, BMH & AF

allSubjLeftQuant = zeros([length(selSizes), length(selVarExps), length(paths)]);
allSubjRightQuant = zeros([length(selSizes), length(selVarExps), length(paths)]);
for nSubj = whichSubs %1:length(paths)
    fname = dir(sprintf('dataClustering_Subject_%s_Date*', num2str(nSubj)));
    load(fname.name)
    allSubjLeftQuant(:,:,nSubj) = outputClusteringLeftQuant;
    allSubjRightQuant(:,:,nSubj) = outputClusteringRightQuant;
end

% takes the average clustering result between hemispheres
clustOutput = (allSubjLeftQuant + allSubjRightQuant) / 2;
%clustOutput=allSubjRightQuant;
% takes the average clustering result between participants (or optionally
% choose one participant)
clustOutputMean = clustOutput(:,:,nSubj);  %mean(clustOutput, 3);
% takes the standard error of the mean clustering result
clustOutputSe = std(clustOutput, 0, 3) ./ sqrt( length( paths ) );

% #clustOutput <- abind( outputClusteringSubjLeft, outputClusteringSubjRight, along=3 ) 
% #clustOutputMean <- apply( clustOutput, c(1,2), mean )
% #clustOutputSe <- apply( clustOutput, c(1,2), sd ) / sqrt( dim( clustOutput )[3] )

% performs the linear and piecewise regression 
clusterSize = selSizes; % 86, cluster size (or set of sizes) in mm^2 of cortical surface
% set.seed(123)
% graphics.off()
% x11(width = 5, height = 4.5)
% library(segmented)
%this is the selected clustering size corresponding to selSizes[31] = 86 squared millimeters
lineCode = find(selSizes==clusterSize); 
% y axis (cluster number) limit
maxY = 16; 
%X axis candidates to plot
plotVals = 1:(size(clustOutputMean,2));
figure; plot(selVarExps(plotVals), clustOutputMean(lineCode,plotVals), 'k', 'LineWidth', 2);
axis([selVarExps(plotVals(1)) selVarExps(plotVals(end)) 0 maxY]);
axis square;
ylabel('Number of clusters');
xlabel('Variance explained');
% axis(1, selVarExps[plotVals[seq(1,length(plotVals),15)]], round(selVarExps[plotVals[seq(1,length(plotVals),15)]],3) )
% axis(2, seq(0,maxY,10), seq(0,maxY,10), las=1  )
hold on; plot(selVarExps(plotVals), clustOutputMean(lineCode,plotVals)-clustOutputSe(lineCode,plotVals), 'Color', [0.5 0.5 0.5]);
hold on; plot(selVarExps(plotVals), clustOutputMean(lineCode,plotVals)+clustOutputSe(lineCode,plotVals), 'Color', [0.5 0.5 0.5]);
end