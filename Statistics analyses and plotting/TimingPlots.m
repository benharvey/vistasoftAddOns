%Put everything in format for timing vs distance plots
DTnames=["All", "Odd", "Even"];

for thisSub=whichSubs
    clear plotdata;
    eval(['data=', char(subjectOrder(thisSub)), ';'])
    for whichMap=mapList %1:length(mapNames)
        for whichHemi=1:2
            count=(whichMap-1)*2+whichHemi;
            dataPresent=char(strcat('isfield(data, ''', mapNames(whichMap), ''') && isfield(data.',  mapNames(whichMap), ',''',hemispheres(whichHemi), ''')' ));
            if eval(dataPresent)
                for whichDT=1:3;
                    dataName=char(strcat('data.', mapNames(whichMap), '.', hemispheres(whichHemi),'.DurationPeriod.', DTnames(whichDT)));
                    eval(['dataTmp=', dataName, ';'])
                    dataTmp.ROItitle=char(strcat(subjectOrder(thisSub),mapNames(whichMap), hemispheres(whichHemi)));
                    
                    if thisSub==9
                        dataTmp.meanDist=25;
                    end
                    
                    plotdata{count, whichDT}=dataTmp(1,1);
                end
            end
        end
    end
    %try
    stats{thisSub}=PlotRoiDistanceAllRoiDurationPeriod(plotdata, 0.1, whichPlots); %Plot everything [1 1 1 0 1 1 1] %Tuning widths only [0 0 0 0 0 1 1]
    %end
end