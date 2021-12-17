function [dataAll] = MotionIntegrationPlotterAcrossVsAlong(dataAll,SD, BK, MS)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
for roiCounter=1:size(dataAll,2);
    sigmasIn=[dataAll{roiCounter}.sigma(dataAll{roiCounter}.ii,4)-dataAll{roiCounter}.sigma(dataAll{roiCounter}.ii,2) dataAll{roiCounter}.sigma(dataAll{roiCounter}.ii,7)-dataAll{roiCounter}.sigma(dataAll{roiCounter}.ii,6)];
    varexpIn=[dataAll{roiCounter}.varexp(dataAll{roiCounter}.ii,4)+dataAll{roiCounter}.varexp(dataAll{roiCounter}.ii,2)./2 dataAll{roiCounter}.varexp(dataAll{roiCounter}.ii,7)+dataAll{roiCounter}.varexp(dataAll{roiCounter}.ii,6)./2];
    eccsIn=[dataAll{roiCounter}.ecc(dataAll{roiCounter}.ii,4)-dataAll{roiCounter}.ecc(dataAll{roiCounter}.ii,2) dataAll{roiCounter}.ecc(dataAll{roiCounter}.ii,7)-dataAll{roiCounter}.ecc(dataAll{roiCounter}.ii,6)];
    s2In=[dataAll{roiCounter}.s2(dataAll{roiCounter}.ii,4)-dataAll{roiCounter}.s2(dataAll{roiCounter}.ii,2) dataAll{roiCounter}.s2(dataAll{roiCounter}.ii,7)-dataAll{roiCounter}.s2(dataAll{roiCounter}.ii,6)];
    varexpCompIn=[dataAll{roiCounter}.varexp(dataAll{roiCounter}.ii,4)-dataAll{roiCounter}.varexp(dataAll{roiCounter}.ii,2) dataAll{roiCounter}.varexp(dataAll{roiCounter}.ii,7)-dataAll{roiCounter}.varexp(dataAll{roiCounter}.ii,6)];
    
[dataAll{roiCounter}.sigmaP, dataAll{roiCounter}.sigmaMeanDif, dataAll{roiCounter}.sigmaMean, dataAll{roiCounter}.sigmaStErr] = compareParameterSignificance(sigmasIn, varexpIn);
[dataAll{roiCounter}.eccP, dataAll{roiCounter}.eccMeanDif, dataAll{roiCounter}.eccMean, dataAll{roiCounter}.eccStErr] = compareParameterSignificance(eccsIn, varexpIn);
[dataAll{roiCounter}.s2P, dataAll{roiCounter}.s2MeanDif, dataAll{roiCounter}.s2Mean, dataAll{roiCounter}.s2StErr] = compareParameterSignificance(s2In, varexpIn);
[dataAll{roiCounter}.veP, dataAll{roiCounter}.veMeanDif, dataAll{roiCounter}.veMean, dataAll{roiCounter}.veStErr] = compareParameterSignificance(varexpCompIn, varexpIn);

sigmaMeans(roiCounter,1)=mean(dataAll{roiCounter}.sigma(dataAll{roiCounter}.ii,4)+dataAll{roiCounter}.sigma(dataAll{roiCounter}.ii,2)./2);
sigmaMeans(roiCounter,2)=mean(dataAll{roiCounter}.sigma(dataAll{roiCounter}.ii,7)+dataAll{roiCounter}.sigma(dataAll{roiCounter}.ii,6)./2);

sigmaErrs(roiCounter,1)=std(dataAll{roiCounter}.sigma(dataAll{roiCounter}.ii,4)+dataAll{roiCounter}.sigma(dataAll{roiCounter}.ii,2)./2)/sqrt(length(dataAll{roiCounter}.sigma(dataAll{roiCounter}.ii,2)./2)/8);
sigmaErrs(roiCounter,2)=std(dataAll{roiCounter}.sigma(dataAll{roiCounter}.ii,7)+dataAll{roiCounter}.sigma(dataAll{roiCounter}.ii,6)./2)/sqrt(length(dataAll{roiCounter}.sigma(dataAll{roiCounter}.ii,6)./2)/8);

s=wstat(dataAll{roiCounter}.sigma(dataAll{roiCounter}.ii,2), dataAll{roiCounter}.varexp(dataAll{roiCounter}.ii,2),8);
sigmaMeanConditions(roiCounter,1)=s.mean;
sigmaErrConditions(roiCounter,1)=s.sterr;
s=wstat(dataAll{roiCounter}.sigma(dataAll{roiCounter}.ii,4), dataAll{roiCounter}.varexp(dataAll{roiCounter}.ii,4),8);
sigmaMeanConditions(roiCounter,2)=s.mean;
sigmaErrConditions(roiCounter,2)=s.sterr;
s=wstat(dataAll{roiCounter}.sigma(dataAll{roiCounter}.ii,6), dataAll{roiCounter}.varexp(dataAll{roiCounter}.ii,6),8);
sigmaMeanConditions(roiCounter,3)=s.mean;
sigmaErrConditions(roiCounter,3)=s.sterr;
s=wstat(dataAll{roiCounter}.sigma(dataAll{roiCounter}.ii,7), dataAll{roiCounter}.varexp(dataAll{roiCounter}.ii,7),8);
sigmaMeanConditions(roiCounter,4)=s.mean;
sigmaErrConditions(roiCounter,4)=s.sterr;

s=wstat(dataAll{roiCounter}.ecc(dataAll{roiCounter}.ii,2), dataAll{roiCounter}.varexp(dataAll{roiCounter}.ii,2),8);
eccMeanConditions(roiCounter,1)=s.mean;
eccErrConditions(roiCounter,1)=s.sterr;
s=wstat(dataAll{roiCounter}.ecc(dataAll{roiCounter}.ii,4), dataAll{roiCounter}.varexp(dataAll{roiCounter}.ii,4),8);
eccMeanConditions(roiCounter,2)=s.mean;
eccErrConditions(roiCounter,2)=s.sterr;
s=wstat(dataAll{roiCounter}.ecc(dataAll{roiCounter}.ii,6), dataAll{roiCounter}.varexp(dataAll{roiCounter}.ii,6),8);
eccMeanConditions(roiCounter,3)=s.mean;
eccErrConditions(roiCounter,3)=s.sterr;
s=wstat(dataAll{roiCounter}.ecc(dataAll{roiCounter}.ii,7), dataAll{roiCounter}.varexp(dataAll{roiCounter}.ii,7),8);
eccMeanConditions(roiCounter,4)=s.mean;
eccErrConditions(roiCounter,4)=s.sterr;
end

for n=1:size(dataAll,2);
    sigmaMeanAll(n,:)=dataAll{n}.sigmaMean(:);
    sigmaErrAll(n,:)=dataAll{n}.sigmaStErr(:);
    eccMeanAll(n,:)=dataAll{n}.eccMean(:);
    eccErrAll(n,:)=dataAll{n}.eccStErr(:);
    s2MeanAll(n,:)=dataAll{n}.s2Mean(:);
    s2ErrAll(n,:)=dataAll{n}.s2StErr(:);
    veMeanAll(n,:)=dataAll{n}.veMean(:);
    veErrAll(n,:)=dataAll{n}.veStErr(:);  
end

figure; bar(sigmaMeanAll); colormap gray; hold on; errorbar(sigmaMeanAll, sigmaErrAll, 'k.');
axis([0.5 size(dataAll,2)-0.5 -1.5 2.5]);
axis square
figure; bar(eccMeanAll); colormap gray; hold on; errorbar(eccMeanAll, eccErrAll, 'k.');
axis([0.5 size(dataAll,2)-0.5 -1.5 2.5]);
axis square
% figure; bar(s2MeanAll); colormap gray; hold on; errorbar(s2MeanAll, s2ErrAll, 'k.');
% axis([0.5 size(dataAll,2)-0.5 -10 10]);
% axis square
figure; bar(veMeanAll); colormap gray; hold on; errorbar(veMeanAll, veErrAll, 'k.');
axis([0.5 size(dataAll,2)-0.5 -0.1 0.2]);
axis square

%normalized to mean pRF size
figure; bar(sigmaMeanAll./sigmaMeans); colormap gray; hold on; errorbar(sigmaMeanAll./sigmaMeans, sigmaErrAll./sigmaMeans, 'k.');
axis([0.5 size(dataAll,2)-0.5 -0.5 0.5]);
axis square
figure; bar(eccMeanAll./sigmaMeans); colormap gray; hold on; errorbar(eccMeanAll./sigmaMeans, eccErrAll./sigmaMeans, 'k.');
axis([0.5 size(dataAll,2)-0.5 -0.5 0.5]);
axis square

%Raw pRF sizes, not differences
figure; bar(sigmaMeanConditions(:,3:4)); colormap gray; hold on; errorbar(sigmaMeanConditions(:,3:4), sigmaErrConditions(:,3:4), 'k.');
axis([0.5 size(dataAll,2)-0.5 0 5]);
axis square
figure; bar(eccMeanConditions(:,3:4)); colormap gray; hold on; errorbar(eccMeanConditions(:,3:4), eccErrConditions(:,3:4), 'k.');
axis([0.5 size(dataAll,2)-0.5 0 5]);
axis square




% figure; bar(eccMeanAll./repmat(mean(eccMeanAll,2), [1 4])); colormap gray; hold on; errorbar(eccMeanAll./repmat(mean(eccMeanAll,2), [1 4]), eccErrAll./repmat(mean(eccMeanAll,2), [1 4]), 'k.');
% axis([0.5 size(dataAll,2)-0.5 0.5 1.5]);
% axis square
% figure; bar(s2MeanAll./repmat(mean(s2MeanAll,2), [1 4])); colormap gray; hold on; errorbar(s2MeanAll./repmat(mean(s2MeanAll,2), [1 4]), s2ErrAll./repmat(mean(s2MeanAll,2), [1 4]), 'k.');
% axis([0.5 size(dataAll,2)-0.5 0.5 1.5]);
% axis square

if exist('SD','var') && ~isempty(SD)
for n=1:size(dataAll,2);
    sigmaMeanAll(n,:)=SD{n}.sigmaMean([1:4]);
    sigmaErrAll(n,:)=SD{n}.sigmaStErr([1:4]);
    eccMeanAll(n,:)=SD{n}.eccMean([1:4]);
    eccErrAll(n,:)=SD{n}.eccStErr([1:4]);
    s2MeanAll(n,:)=SD{n}.s2Mean([1:4]);
    s2ErrAll(n,:)=SD{n}.s2StErr([1:4]);
end
h1=figure; subplot(2,2,1); bar(sigmaMeanAll); colormap gray; hold on; errorbar(sigmaMeanAll, sigmaErrAll, 'k.');
axis([0.5 size(dataAll,2)-0.5 0 8]);
axis square

h2=figure; subplot(2,2,1); bar(eccMeanAll); colormap gray; hold on; errorbar(eccMeanAll, eccErrAll, 'k.');
axis([0.5 size(dataAll,2)-0.5 0 8]);
axis square

end


if exist('BK','var') && ~isempty(BK)
for n=1:size(dataAll,2);
    sigmaMeanAll(n,:)=BK{n}.sigmaMean([1:4]);
    sigmaErrAll(n,:)=BK{n}.sigmaStErr([1:4]);
    eccMeanAll(n,:)=BK{n}.eccMean([1:4]);
    eccErrAll(n,:)=BK{n}.eccStErr([1:4]);
    s2MeanAll(n,:)=BK{n}.s2Mean([1:4]);
    s2ErrAll(n,:)=BK{n}.s2StErr([1:4]);
end
figure(h1); hold on; subplot(2,2,2); bar(sigmaMeanAll); colormap gray; hold on; errorbar(sigmaMeanAll, sigmaErrAll, 'k.');
axis([0.5 size(dataAll,2)-0.5 0 8]);
axis square

figure(h2); hold on; subplot(2,2,2); bar(eccMeanAll); colormap gray; hold on; errorbar(eccMeanAll, eccErrAll, 'k.');
axis([0.5 size(dataAll,2)-0.5 0 8]);
axis square

if exist('MS','var') && ~isempty(MS)
for n=1:size(dataAll,2);
    sigmaMeanAll(n,:)=MS{n}.sigmaMean([1:4]);
    sigmaErrAll(n,:)=MS{n}.sigmaStErr([1:4]);
    eccMeanAll(n,:)=MS{n}.eccMean([1:4]);
    eccErrAll(n,:)=MS{n}.eccStErr([1:4]);
    s2MeanAll(n,:)=MS{n}.s2Mean([1:4]);
    s2ErrAll(n,:)=MS{n}.s2StErr([1:4]);
end
figure(h1); hold on; subplot(2,2,3); bar(sigmaMeanAll); colormap gray; hold on; errorbar(sigmaMeanAll, sigmaErrAll, 'k.');
axis([0.5 size(dataAll,2)-0.5 0 8]);
axis square

figure(h2); hold on; subplot(2,2,3); bar(eccMeanAll); colormap gray; hold on; errorbar(eccMeanAll, eccErrAll, 'k.');
axis([0.5 size(dataAll,2)-0.5 0 8]);
axis square
end

if exist('WZ','var') && ~isempty(WZ)
for n=1:size(dataAll,2);
    sigmaMeanAll(n,:)=WZ{n}.sigmaMean([1:4]);
    sigmaErrAll(n,:)=WZ{n}.sigmaStErr([1:4]);
    eccMeanAll(n,:)=WZ{n}.eccMean([1:4]);
    eccErrAll(n,:)=WZ{n}.eccStErr([1:4]);
    s2MeanAll(n,:)=WZ{n}.s2Mean([1:4]);
    s2ErrAll(n,:)=WZ{n}.s2StErr([1:4]);
end
figure(h1); hold on; subplot(2,2,4); bar(sigmaMeanAll); colormap gray; hold on; errorbar(sigmaMeanAll, sigmaErrAll, 'k.');
axis([0.5 size(dataAll,2)-0.5 0 8]);
axis square

figure(h2); hold on; subplot(2,2,4); bar(eccMeanAll); colormap gray; hold on; errorbar(eccMeanAll, eccErrAll, 'k.');
axis([0.5 size(dataAll,2)-0.5 0 8]);
axis square
end

end


