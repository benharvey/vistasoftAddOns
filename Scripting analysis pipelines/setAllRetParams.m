function setAllRetParams(scan, destScans)
%sets retinotopy params for all datatypes to be the same as those in 'scan'
%Useful when running many models from the same parameters

global dataTYPES

if ~exist('destScans','var') || isempty(destScans)
    destScans=1:length(dataTYPES);
end
for n=destScans
    if isstruct(scan)
        dataTYPES(n).retinotopyModelParams=scan;
    elseif double(n)~=double(scan)
        dataTYPES(n).retinotopyModelParams=dataTYPES(scan).retinotopyModelParams;
    end
end
saveSession;
end

