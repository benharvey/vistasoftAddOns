%function combineOneG_DoG_All(path)
a=dir('*.mat');
for n=1:length(a)
    a(n).name=strcat(pwd, '/', a(n).name);
end
if mod(length(a), 2)>0
    disp('error, wrong number of files')
else
    for n=1:length(a)/2
        combineOneG_DoG(a(n*2-1).name, a(n*2).name, VOLUME{1});
        n
        a(n*2).name
        a(n*2-1).name
    end
end


%end

