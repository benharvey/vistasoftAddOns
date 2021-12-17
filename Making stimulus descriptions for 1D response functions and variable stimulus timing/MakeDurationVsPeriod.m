%MAKE MODEL PARAMETER AND IMAGE DEFINITIONS FOR TIMING-SELECTIVE PRF MODELS

%Set Variable
ConLumPeriod=[0.05:0.05:1 2.1];
ConLumDuration=[0.05:0.05:1 2];
ConLumOffDur=ConLumPeriod-ConLumDuration;
ConLumOccupancy=ConLumDuration./ConLumPeriod;
ConDurPeriod=[0.05:0.05:1 2.1];
ConDurDuration=repmat(0.05, size(ConDurPeriod));
ConDurOffDur=ConDurPeriod-ConDurDuration;
ConDurOccupancy=ConDurDuration./ConDurPeriod;
ConPerPeriod=repmat(1, size(ConDurPeriod));
ConPerPeriod(end)=2.1;
ConPerDuration=[0.05:0.05:1 2];
ConPerOffDur=ConPerPeriod-ConPerDuration;
ConPerOccupancy=ConPerDuration./ConPerPeriod;
GapsDuration=[0.05:0.05:0.5 0.05:0.05:0.5 0.05];
GapsPeriod=[0.95:-0.05:0.50 0.55:0.05:1 2.1]; 
GapsOffDur=GapsPeriod-GapsDuration;
GapsOccupancy=GapsDuration./GapsPeriod;


%Plot linear figure
figure; hold on; plot(ConLumDuration, ConLumPeriod, 'ro')
axis square
hold on; plot(ConPerDuration, ConPerPeriod, 'ko')
hold on; plot(ConDurDuration, ConDurPeriod, 'go')
hold on; plot(GapsDuration, GapsPeriod, 'bo')

%Plot log figure
figure; hold on; plot(log(ConLumDuration), log(ConLumPeriod), 'ro')
axis square
hold on; plot(log(ConPerDuration), log(ConPerPeriod), 'ko')
hold on; plot(log(ConDurDuration), log(ConDurPeriod), 'go')
hold on; plot(log(GapsDuration), log(GapsPeriod), 'bo')
ax=gca;
ax.XTick=log(ConLumDuration)
ax.YTick=log(ConLumDuration)
axis([-3.1 log(2) -3.1 log(2)])

%Plot linear on time-off time figure
figure; hold on; plot(ConLumDuration, ConLumOffDur, 'ro')
axis square
hold on; plot(ConPerDuration, ConPerOffDur, 'ko')
hold on; plot(ConDurDuration, ConDurOffDur, 'go')
hold on; plot(GapsDuration, GapsOffDur, 'bo')

%Plot linear occupancy period figure
figure; hold on; plot(ConLumOccupancy, ConLumPeriod, 'ro')
axis square
hold on; plot(ConPerOccupancy, ConPerPeriod, 'ko')
hold on; plot(ConDurOccupancy, ConDurPeriod, 'go')
hold on; plot(GapsOccupancy, GapsPeriod, 'bo')

%Plot linear occupancy duration figure (NOT a useful space)
figure; hold on; plot(1./ConLumOccupancy, ConLumDuration,  'ro')
axis square
hold on; plot(1./ConPerOccupancy, ConPerDuration,  'ko')
hold on; plot(1./ConDurOccupancy, ConDurDuration, 'go')
hold on; plot(1./GapsOccupancy, GapsDuration, 'bo')
axis([0 1 0 1])

%MAKES STIMULUS IMAGES
%Make modelling images (to run 2d pRF models)
ConLumImages=uint8(zeros(40, 40, size(ConDurDuration,2)));
ConDurImages=uint8(zeros(40, 40, size(ConDurDuration,2)));
ConPerImages=uint8(zeros(40, 40, size(ConDurDuration,2)));
GapsImages=uint8(zeros(40, 40, size(ConDurDuration,2)));

for n=1:size(ConDurDuration,2)
    ConLumImages(uint8(ConLumDuration(n).*20), uint8(ConLumPeriod(n).*20), n)=1;
    ConDurImages(uint8(ConDurDuration(n).*20), uint8(ConDurPeriod(n).*20), n)=1;
    ConPerImages(uint8(ConPerDuration(n).*20), uint8(ConPerPeriod(n).*20), n)=1;
    GapsImages(uint8(GapsDuration(n).*20), uint8(GapsPeriod(n).*20), n)=1;
end

for n=1:size(ConDurDuration,2)
    ConLumImages(:,:,n)=ConLumImages(:,:,n)';
    ConDurImages(:,:,n)=ConDurImages(:,:,n)';
    ConPerImages(:,:,n)=ConPerImages(:,:,n)';
    GapsImages(:,:,n)=GapsImages(:,:,n)';
end

%Test the result
figure;
for n=1:size(ConDurDuration,2)
    imagesc([flipud(ConLumImages(:,:,n)) ones(40,1) flipud(ConPerImages(:,:,n)) ones(40,1) flipud(ConDurImages(:,:,n)) ones(40,1) flipud(GapsImages(:,:,n))])
    axis image;
    drawnow;
    pause(0.4)
end

images=cat(3, ConLumImages,  ConPerImages, ConDurImages, GapsImages);
images=images([1:20 40],[1:20 40], :);
images(:,:,85)=zeros([21,21]);
images=images.*255;

%IF NEEDED: Changes to params and stimuli for mix modelling (started from separate
%scans)
%First, make the prescan period a mixture of the ends of the three possible
%preceeding scans
images=images.*3;
images(:,:,86)=(images(:,:,21)+images(:,:,21)+images(:,:,42))./3;
images(:,:,87)=(images(:,:,21)+images(:,:,42)+images(:,:,42))./3;
%Then set some fields in params to avoid the assumption of three cycles
params.numCycles=1;
params.scanDuration=params.period+params.prescanDuration;
cyclePrescanFrames=20*2.1*(56+6);
stimulus.seq=stimulus.seq(1:cyclePrescanFrames);
stimulus.seqtiming=stimulus.seqtiming(1:cyclePrescanFrames);
%Then set the prescan stimulus frames correctly
prescanSeq=stimulus.seq(1:252);
prescanSeq(prescanSeq==21)=87;
prescanSeq(prescanSeq==42)=86;
prescanSeq(prescanSeq==63)=87;
prescanSeq(prescanSeq==84)=86;
stimulus.seq(1:252)=prescanSeq;

%MAKE STIMULUS PARAMETER FILES FOR MODELLING
%To convert stimulus sequence to mrVista input
%Cut prescan (if needed, can take sequence before prescan is added)
seq=seq(253:end);
seqTiming=seqTiming(253:end);
newSet=newSet(253:end);

%Convert to onsets and offsets
onsets=zeros(size(newSet));
offsets=zeros(size(newSet));
for n=1:length(onsets)
    if newSet(n)>0
        if newSet(n)==2000
            onsets(n)=21;
            offsets(n+39)=21;
        elseif newSet(n)==50 && sum(newSet(n+1:n+41))==0 && ((n+41)==length(newSet) || newSet(n+42)>0)
            %Watch out, hard coded
            
                
                %(n>7476 && n<7602) || (n>8022 && n<8129) || (n>8568 && n<8654) || n<9114
               onsets(n)=21;
               offsets(n+1)=21;
        else
            onsets(n)=newSet(n)/50;
            offsets(n+newSet(n)/50-1)=newSet(n)/50;
        end
    end
end

onsetsConLum=onsets(1:2352);
onsetsConPer=onsets(2353:4704);
onsetsConDur=onsets(4705:7056);
onsetsGaps=onsets(7057:9408);
SecondSweep=(onsetsGaps>0 & onsetsGaps<21);
SecondSweep(1:546)=0;
SecondSweep(1638:end)=0;
onsetsGaps(SecondSweep)=onsetsGaps(SecondSweep)+10;


offsetsConLum=offsets(1:2352);
offsetsConPer=offsets(2353:4704);
offsetsConDur=offsets(4705:7056);
offsetsGaps=offsets(7057:9408);
SecondSweep=(offsetsGaps>0 & offsetsGaps<21);
SecondSweep(1:546)=0;
SecondSweep(1638:end)=0;
offsetsGaps(SecondSweep)=offsetsGaps(SecondSweep)+10;

%Check result
figure; 
axis image;
colormap gray;
for n=7056:length(onsets)
    if n<=2352
        if onsetsConLum(n)>0
            imagesc(flipud(ConLumImages(:,:,onsetsConLum(n))), [0 1])
        else
            imagesc(zeros(40,40), [0 1])
        end
    elseif n<=4704
        if onsetsConDur(n-2352)>0
            imagesc(flipud(ConDurImages(:,:,onsetsConDur(n-2352))), [0 1])
        else
            imagesc(zeros(40,40), [0 1])
        end     
    elseif n<=7056
        if onsetsConPer(n-4704)>0
            imagesc(flipud(ConPerImages(:,:,onsetsConPer(n-4704))), [0 1])
        else
            imagesc(zeros(40,40), [0 1])
        end 
    else
        if onsetsGaps(n-7056)>0
            imagesc(flipud(GapsImages(:,:,onsetsGaps(n-7056))), [0 1])
        else
            imagesc(zeros(40,40), [0 1])
        end 
    end
    drawnow;
    %pause(0.001)
end

%Convert onsets and offsets to full mrVista seq, SEPARATE SCANS
seqCLon=[onsetsConLum((end-251):end) onsetsConLum onsetsConLum onsetsConLum];
seqCPon=[onsetsConPer((end-251):end) onsetsConPer onsetsConPer onsetsConPer];
seqCDon=[onsetsConDur((end-251):end) onsetsConDur onsetsConDur onsetsConDur];
seqGPon=[onsetsGaps((end-251):end) onsetsGaps onsetsGaps onsetsGaps];
seqCLoff=[offsetsConLum((end-251):end) offsetsConLum offsetsConLum offsetsConLum];
seqCPoff=[offsetsConPer((end-251):end) offsetsConPer offsetsConPer offsetsConPer];
seqCDoff=[offsetsConDur((end-251):end) offsetsConDur offsetsConDur offsetsConDur];
seqGPoff=[offsetsGaps((end-251):end) offsetsGaps offsetsGaps offsetsGaps];

seqCPon(seqCPon>0)=seqCPon(seqCPon>0)+21;
seqCPoff(seqCPoff>0)=seqCPoff(seqCPoff>0)+21;
seqCDon(seqCDon>0)=seqCDon(seqCDon>0)+42;
seqCDoff(seqCDoff>0)=seqCDoff(seqCDoff>0)+42;
seqGPon(seqGPon>0)=seqGPon(seqGPon>0)+63;
seqGPoff(seqGPoff>0)=seqGPoff(seqGPoff>0)+63;

seqCLon(seqCLon==0)=85;
seqCLoff(seqCLoff==0)=85;
seqCPon(seqCPon==0)=85;
seqCPoff(seqCPoff==0)=85;
seqCDon(seqCDon==0)=85;
seqCDoff(seqCDoff==0)=85;
seqGPon(seqGPon==0)=85;
seqGPoff(seqGPoff==0)=85;

%Save as params files
load('params_lin_duration.mat', 'params')
stimulus.seqtiming=seqTiming(1:length(seqCLon));
stimulus.seq=seqCLon;
save('params_ConLumOn2D.mat', 'params', 'stimulus')
stimulus.seq=seqCLoff;
save('params_ConLumOff2D.mat', 'params', 'stimulus')

stimulus.seq=seqCDon;
save('params_ConDurOn2D.mat', 'params', 'stimulus')
stimulus.seq=seqCDoff;
save('params_ConDurOff2D.mat', 'params', 'stimulus')

stimulus.seq=seqCPon;
save('params_ConPerOn2D.mat', 'params', 'stimulus')
stimulus.seq=seqCPoff;
save('params_ConPerOff2D.mat', 'params', 'stimulus')

stimulus.seq=seqGPon;
save('params_TimeGapsOn2D.mat', 'params', 'stimulus')
stimulus.seq=seqGPoff;
save('params_TimeGapsOff2D.mat', 'params', 'stimulus')

%Convert onsets and offsets to full mrVista seq, MIXED SCANS
onsetsConPer(onsetsConPer>0)=onsetsConPer(onsetsConPer>0)+21;
onsetsConDur(onsetsConDur>0)=onsetsConDur(onsetsConDur>0)+42;
onsetsGaps(onsetsGaps>0)=onsetsGaps(onsetsGaps>0)+63;
offsetsConPer(offsetsConPer>0)=offsetsConPer(offsetsConPer>0)+21;
offsetsConDur(offsetsConDur>0)=offsetsConDur(offsetsConDur>0)+42;
offsetsGaps(offsetsGaps>0)=offsetsGaps(offsetsGaps>0)+63;

onsetsConLum(onsetsConLum==0)=85;
onsetsConPer(onsetsConPer==0)=85;
onsetsConDur(onsetsConDur==0)=85;
onsetsGaps(onsetsGaps==0)=85;
offsetsConLum(offsetsConLum==0)=85;
offsetsConPer(offsetsConPer==0)=85;
offsetsConDur(offsetsConDur==0)=85;
offsetsGaps(offsetsGaps==0)=85;

seqCLonSame=[onsetsConLum((end-251):end) onsetsConLum];
seqCPonSame=[onsetsConPer((end-251):end) onsetsConPer];
seqCDonSame=[onsetsConDur((end-251):end) onsetsConDur];
seqGPonSame=[onsetsGaps((end-251):end) onsetsGaps];
seqCLonDif=[onsetsConDur((end-251):end) onsetsConLum];
seqCPonDif=[onsetsConDur((end-251):end) onsetsConPer];
seqCDonDif=[onsetsConLum((end-251):end) onsetsConDur];
seqGPonDif=[onsetsConLum((end-251):end) onsetsGaps];

seqCLoffSame=[offsetsConLum((end-251):end) offsetsConLum];
seqCPoffSame=[offsetsConPer((end-251):end) offsetsConPer];
seqCDoffSame=[offsetsConDur((end-251):end) offsetsConDur];
seqGPoffSame=[offsetsGaps((end-251):end) offsetsGaps];
seqCLoffDif=[offsetsConDur((end-251):end) offsetsConLum];
seqCPoffDif=[offsetsConDur((end-251):end) offsetsConPer];
seqCDoffDif=[offsetsConLum((end-251):end) offsetsConDur];
seqGPoffDif=[offsetsConLum((end-251):end) offsetsGaps];

%Save as params files
load('params_lin_duration.mat', 'params')
stimulus.seqtiming=seqTiming(1:length(seqCLonSame));
stimulus.seq=seqCLonSame;
save('params_ConLumOnSame2D.mat', 'params', 'stimulus')
stimulus.seq=seqCLoffSame;
save('params_ConLumOffSame2D.mat', 'params', 'stimulus')

stimulus.seq=seqCDonSame;
save('params_ConDurOnSame2D.mat', 'params', 'stimulus')
stimulus.seq=seqCDoffSame;
save('params_ConDurOffSame2D.mat', 'params', 'stimulus')

stimulus.seq=seqCPonSame;
save('params_ConPerOnSame2D.mat', 'params', 'stimulus')
stimulus.seq=seqCPoffSame;
save('params_ConPerOffSame2D.mat', 'params', 'stimulus')

stimulus.seq=seqGPonSame;
save('params_TimeGapsOnSame2D.mat', 'params', 'stimulus')
stimulus.seq=seqGPoffSame;
save('params_TimeGapsOffSame2D.mat', 'params', 'stimulus')

save('params_ConLumOnDif2D.mat', 'params', 'stimulus')
stimulus.seq=seqCLoffDif;
save('params_ConLumOffDif2D.mat', 'params', 'stimulus')

stimulus.seq=seqCDonDif;
save('params_ConDurOnDif2D.mat', 'params', 'stimulus')
stimulus.seq=seqCDoffDif;
save('params_ConDurOffDif2D.mat', 'params', 'stimulus')

stimulus.seq=seqCPonDif;
save('params_ConPerOnDif2D.mat', 'params', 'stimulus')
stimulus.seq=seqCPoffDif;
save('params_ConPerOffDif2D.mat', 'params', 'stimulus')

stimulus.seq=seqGPonDif;
save('params_TimeGapsOnDif2D.mat', 'params', 'stimulus')
stimulus.seq=seqGPoffDif;
save('params_TimeGapsOffDif2D.mat', 'params', 'stimulus')

% To make temporal frequency transforms of stimulus sequence
load('params_ConLumOn2D.mat')
stimulusConLumOn=stimulus;
load('params_ConLumOff2D.mat')
stimulusConLumOff=stimulus;

stimOn=stimulusConLumOn.seq~=85;
stimOff=stimulusConLumOff.seq~=85;
stimOnUS=zeros(1,length(stimOn)*10);
stimOffUS=zeros(1,length(stimOff)*10);
stimOnUS(10*find(stimOn)-9)=1;
stimOffUS(10*find(stimOff)-1)=1;
stimOnUS=[stimOnUS stimOnUS];
stimOffUS=[stimOffUS stimOffUS];

onsets=find(stimOnUS);
offsets=find(stimOffUS);
%fftsOut=zeros(420,(length(onsets)/2));

for whichEvent=31:(length(onsets)/2)
    cycle=zeros(1,offsets(whichEvent+(length(onsets)/2))-offsets(whichEvent-1+(length(onsets)/2)));
    cycle((end-(offsets(whichEvent+(length(onsets)/2))-onsets(whichEvent+(length(onsets)/2)))):end)=1;
    tmpfft=abs(fft(cycle));
    if length(tmpfft)<420;
        tmpfft(end+1:420)=0;
    end
    if ~exist('fftsOut', 'var') || isempty(fftsOut)
        fftsOut(1,:)=tmpfft;
        transformList(whichEvent)=1;
    else
        matchFound=0;
        for n=1:size(fftsOut,1)
            if all(tmpfft==fftsOut(n,:))
                transformList(whichEvent)=n;
                matchFound=1;
                break;
            end
        end
        if matchFound==0
            transformList(whichEvent)=size(fftsOut,1)+1;
            fftsOut((size(fftsOut,1)+1),:)=tmpfft;
        end
    end
end
fftsOut=zeros(420,sum(stimConLumOn));
