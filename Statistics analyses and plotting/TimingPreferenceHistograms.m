

%Histograms of duration and period preferences
subjectOrder={'dataAllSub'};
mapNames={'TLO', 'TTOP', 'TTOA', 'TPO', 'TLS', 'TPCI', 'TPCM', 'TPCS', 'TFI' 'TFS'};
hemispheres={'Left', 'Right'};

%Get data locations
veThresh=0.2;
whichSub=1;
figure; hold on;
for whichMap=1:length(mapNames)
    for whichHemi=1:2;
            data=eval(strcat(subjectOrder{whichSub}, '.', mapNames{whichMap}, '.', hemispheres{whichHemi}, '.DurationPeriodHRF.All'));
            subplot(length(hemispheres),length(mapNames),(whichHemi-1)*length(mapNames)+whichMap);
            histogram(data.y0s(data.y0s>0.05 & data.y0s<1 & data.ves>veThresh), 0.05:0.05:1, 'FaceColor', 'none');
            axis square;
            title(strcat(hemispheres(whichHemi), ' ', mapNames(whichMap)))
    end
end