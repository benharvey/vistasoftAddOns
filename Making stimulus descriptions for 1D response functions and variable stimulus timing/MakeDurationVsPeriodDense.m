%MAKE MODEL PARAMETER AND IMAGE DEFINITIONS FOR TIMING-SELECTIVE PRF MODELS
%WITH DENSE SAMPLING OF THE STIMULUS SPACE

%DENSER COVERAGE OF DURATION PERIOD SPACE
TRcounter=0;
clear duration;
clear period;
for sweep=1:11
    for sweepTR=1:10
        TRcounter=TRcounter+1;
        if sweepTR==1
            periodTmp=(sweep-1)*100;
            durationTmp=0;
            movePeriod=0;
        end
        if durationTmp>=periodTmp && movePeriod==0
            movePeriod=1;
            durationTmp=durationTmp+50;
            periodTmp=periodTmp-50;
        end
        if movePeriod==0;
            durationTmp=durationTmp+100;
        else
            periodTmp=periodTmp+100;
        end
        duration(TRcounter)=durationTmp;
        period(TRcounter)=periodTmp;
    end
end
figure; plot(duration, period, '.-')
hold on;
for blank=1:6
    TRcounter=TRcounter+1;
    duration(TRcounter)=2000;
    period(TRcounter)=2000;
end
TRcountRev=TRcounter+1;
for sweep=1:10
    for sweepTR=1:11
        TRcounter=TRcounter+1;
        if sweepTR==1
            periodTmp=1100;
            durationTmp=1000-(sweep-1)*100;
            moveDuration=0;
        end
        if durationTmp>=periodTmp && moveDuration==0
            moveDuration=1;
            durationTmp=durationTmp+50;
            periodTmp=periodTmp-50;
        end
        if moveDuration==0;
            periodTmp=periodTmp-100;
        else
            durationTmp=durationTmp-100;
        end
        duration(TRcounter)=durationTmp;
        period(TRcounter)=periodTmp;
    end
end
plot(duration(TRcountRev:end)-5, period(TRcountRev:end), 'r.-')
for blank=1:6
    TRcounter=TRcounter+1;
    duration(TRcounter)=50;
    period(TRcounter)=2000;
end