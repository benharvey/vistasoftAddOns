function TimingANOVAs(dataS8, dataS9, dataS10, dataS11, dataS12, dataS13, dataS1, dataS3)
subjectOrder={'dataS8', 'dataS9', 'dataS10', 'dataS11', 'dataS12', 'dataS13', 'dataS1', 'dataS3'};
mapNames={'TLO', 'TTOP', 'TTOA', 'TPO', 'TLS', 'TPCI', 'TPCM', 'TPCS', 'TFI' 'TFS'};
hemispheres={'Left', 'Right'};

VEs=nan([length(subjectOrder) length(mapNames) length(hemispheres)]);
Exps=VEs;
SigmaMajor=VEs;
SigmaMinor=VEs;
SigmaRatio=VEs;
SigmaTheta=VEs;
Q1D=VEs;
Q2D=VEs;
Q3D=VEs;
IQRD=VEs;
Q1P=VEs;
Q2P=VEs;
Q3P=VEs;
IQRP=VEs;
linearSlopeMajor=VEs;
vSlopeMajor=VEs;
linearSlopeMinor=VEs;
vSlopeMinor=VEs;
nVoxels=VEs;
%figure; hold on;
for n=1:length(subjectOrder)
    for m=1:length(mapNames)
        for whichHemi=1:length(hemispheres)
            if eval(char(strcat('isfield(', subjectOrder{n}, ', ''', mapNames{m}, ''') && isfield(', subjectOrder{n}, '.', mapNames{m}, ',''',hemispheres{whichHemi}, ''')' ))) 
                
                eval(char(strcat('data=', subjectOrder{n},'.', mapNames{m},'.', hemispheres{whichHemi},'.DurationPeriodHRF.All',';')));
                
                veIndices=data.ves>0.1 & data.x0s>0.05 & data.x0s<1 ;
                if any(veIndices)
                    nVoxels(n,m,whichHemi)=sum(veIndices);
                end
                VEs(n,m,whichHemi)=mean(data.ves(veIndices));
                Exps(n,m,whichHemi)=mean(data.exp(veIndices));
                SigmaMajor(n,m,whichHemi)=mean(data.sigmas(veIndices));
                SigmaMinor(n,m,whichHemi)=mean(data.sigmaMinor(veIndices));
                SigmaRatio(n,m,whichHemi)=mean(data.sigmas(veIndices)./data.sigmaMinor(veIndices));
                SigmaTheta(n,m,whichHemi)=mean(data.sigmaTheta(veIndices));
                Q1D(n,m,whichHemi)=prctile(data.x0s(veIndices), 25);
                Q2D(n,m,whichHemi)=mean(data.x0s(veIndices));
                Q3D(n,m,whichHemi)=prctile(data.x0s(veIndices), 75);
                IQRD(n,m,whichHemi)=prctile(data.x0s(veIndices), 75)-prctile(data.x0s(veIndices), 25);
                
                
                Q1P(n,m,whichHemi)=prctile(data.y0s(veIndices), 25);
                Q2P(n,m,whichHemi)=mean(data.y0s(veIndices));
                Q3P(n,m,whichHemi)=prctile(data.y0s(veIndices), 75);
                IQRP(n,m,whichHemi)=prctile(data.y0s(veIndices), 75)-prctile(data.y0s(veIndices), 25);                                
            end
        end
    end
end

%Setting up ANOVA structures
hemisphereGroups=cat(3, ones(size(VEs(:,:,1))), ones(size(VEs(:,:,1)))*2);
tmp=1:length(subjectOrder);
subjectLabels=cat(3, repmat(tmp(:), [1,length(mapNames)]), repmat(tmp(:), [1,length(mapNames)]));
mapLabels=cat(3, repmat(1:length(mapNames), [length(subjectOrder),1]), repmat(1:length(mapNames), [length(subjectOrder),1]));
mapLabels=mapLabels(:);
subjectLabels=subjectLabels(:);
hemisphereLabels=hemisphereGroups(:);
subjectLabels=subjectOrder(subjectLabels);
mapLabels=mapNames(mapLabels);
hemiTmp={'L', 'R'};
hemisphereLabels=hemiTmp(hemisphereLabels);

%Stats and plot of exponents (general code for 3 factors, use when
%hemisphere difference)
[p, tmp, statsOut] = anovan(Exps(:),{hemisphereLabels subjectLabels mapLabels}, 'varnames', {'hemisphere', 'subject', 'map'});
%Multiple comparison tests on ANOVA output
figure; results=multcompare(statsOut, 'Dimension', [1 3])
axis([0 10 0.5 18.5])
axis square

[p, tmp, statsOut] = anovan(SigmaTheta(:),{subjectLabels mapLabels}, 'varnames', {'subject', 'map'});
%Multiple comparison tests on ANOVA output
figure; results=multcompare(statsOut, 'Dimension', [2])
axis([0 pi/2 0.5 18.5])
axis square

%Plot of VEs
[p, tmp, statsOut] = anovan(VEs(:),{hemisphereLabels subjectLabels mapLabels}, 'varnames', {'hemisphere', 'subject', 'map'});
%Multiple comparison tests on ANOVA output
figure; results=multcompare(statsOut, 'Dimension', [1 3])
axis([0 0.6 0.5 20.5])
axis square

%Plot of ROI voxel count
[p, tmp, statsOut] = anovan(nVoxels(:),{hemisphereLabels subjectLabels mapLabels}, 'varnames', {'hemisphere', 'subject', 'map'});
%Multiple comparison tests on ANOVA output
figure; results=multcompare(statsOut, 'Dimension', [1 3])
axis([0 800 0.5 20.5])
axis square

%Plot of sigma majors
[p, tmp, statsOut] = anovan(SigmaMajor(:),{hemisphereLabels subjectLabels mapLabels}, 'varnames', {'hemisphere', 'subject', 'map'});
%Multiple comparison tests on ANOVA output
figure; results=multcompare(statsOut, 'Dimension', [1 3])
axis([0 1.2 0.5 20.5])
axis square
%Stats, because there is no hemisphere effect
[p, tmp, statsOut] = anovan(SigmaMajor(:),{subjectLabels mapLabels}, 'varnames', {'subject', 'map'});

%Plot of sigma minors
[p, tmp, statsOut] = anovan(SigmaMinor(:),{hemisphereLabels subjectLabels mapLabels}, 'varnames', {'hemisphere', 'subject', 'map'});
%Multiple comparison tests on ANOVA output
figure; results=multcompare(statsOut, 'Dimension', [1 3])
axis([0 1.2 0.5 20.5])
axis square

%Plot of sigma ratios
[p, tmp, statsOut] = anovan(SigmaRatio(:),{hemisphereLabels subjectLabels mapLabels}, 'varnames', {'hemisphere', 'subject', 'map'});
%Multiple comparison tests on ANOVA output
figure; results=multcompare(statsOut, 'Dimension', [1 3])
axis([1 7 0.5 20.5])
axis square
%Stats, because there is no hemisphere effect
[p, tmp, statsOut] = anovan(SigmaRatio(:),{subjectLabels mapLabels}, 'varnames', {'subject', 'map'});

%Plot of sigma thetas
[p, tmp, statsOut] = anovan(0.5*pi-SigmaTheta(:),{hemisphereLabels subjectLabels mapLabels}, 'varnames', {'hemisphere', 'subject', 'map'});
%Multiple comparison tests on ANOVA output
figure; results=multcompare(statsOut, 'Dimension', [1 3])
axis([1/8*pi 3/8*pi 0.5 20.5])
axis square
set(gca,'XTick',[pi/8 2*pi/8 3*pi/8])
%Stats, because there is no hemisphere effect
[p, tmp, statsOut] = anovan(SigmaTheta(:),{subjectLabels mapLabels}, 'varnames', {'subject', 'map'});

%Plot of duration mean
[p, tmp, statsOut] = anovan(Q2D(:),{hemisphereLabels subjectLabels mapLabels}, 'varnames', {'hemisphere', 'subject', 'map'});
%Multiple comparison tests on ANOVA output
figure; results=multcompare(statsOut, 'Dimension', [1 3])
axis([0 1 0.5 20.5])
axis square

%Plot of duration inter quartile range
[p, tmp, statsOut] = anovan(IQRD(:),{hemisphereLabels subjectLabels mapLabels}, 'varnames', {'hemisphere', 'subject', 'map'});
%Multiple comparison tests on ANOVA output
figure; results=multcompare(statsOut, 'Dimension', [1 3])
axis([0 0.5 0.5 20.5])
axis square
%Stats, because there is no hemisphere effect
[p, tmp, statsOut] = anovan(IQRD(:),{subjectLabels mapLabels}, 'varnames', {'subject', 'map'});

%Plot of period inter quartile range
[p, tmp, statsOut] = anovan(IQRP(:),{hemisphereLabels subjectLabels mapLabels}, 'varnames', {'hemisphere', 'subject', 'map'});
%Multiple comparison tests on ANOVA output
figure; results=multcompare(statsOut, 'Dimension', [1 3])
axis([0 0.5 0.5 20.5])
axis square
%Stats, because there is no hemisphere effect
[p, tmp, statsOut] = anovan(IQRP(:),{subjectLabels mapLabels}, 'varnames', {'subject', 'map'});
end