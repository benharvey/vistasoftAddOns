function [correlations, sds, means, ns, dfs] = numPrefCorrelator(data, meanThresh)

indices=data{1}.meanSignal>meanThresh;
for n=1:length(data)
%     data{n}.x0s=log(data{n}.x0s);
    indices=indices & data{n}.x0s>0.01 & data{n}.x0s<2.5 & data{n}.ves>0.3;
end

%indices=indices

correlations=zeros(length(data));
sds=zeros(length(data));
means=zeros(length(data));
ns=zeros(length(data));
dfs=zeros(length(data));

for n=1:length(data)-1
    %correlations(n,n)=1;
    for m=n+1:length(data)
        cor=corrcoef(data{n}.x0s(indices), data{m}.x0s(indices));
        correlations(n,m)=cor(1,2);
        
%         %x0mean=(data{n}.x0s+data{m}.x0s)./2;
%         %rawrss=sum((x0mean(indices)-mean(x0mean(indices))).^2);
%         rawrss1=sum((data{n}.x0s(indices)-mean(data{m}.x0s(indices))).^2);
%         rawrss2=sum((data{m}.x0s(indices)-mean(data{n}.x0s(indices))).^2);
%         rawrss=(rawrss1+rawrss2)/2;
%         rss=sum((data{n}.x0s(indices)-data{m}.x0s(indices)).^2);
%         correlations(n,m)=(1-(rss./rawrss1))+(1-(rss./rawrss2));
%         %correlations(m,n)=1-(rss./rawrss2);
%         
%         correlations=zeros(length(data));
        
        diffs=data{n}.x0s(indices)-data{m}.x0s(indices);
        stats=wstat(diffs, [], 4);
        %[tmp, p, ci, stats]=ttest(diffs);
        
        correlations(m,n)=abs(stats.tval);
        sds(m,n)=stats.stdev;
        means(m,n)=stats.mean;
        ns(m,n)=stats.n;
        dfs(m,n)=stats.df;
        

        
    end
end
%correlations(length(data),length(data))=1;
end

