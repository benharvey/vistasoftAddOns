%Run initial models, finding best fitting HRF
VOLUME{1} = rmRunNumbersScriptLogMonotonic(VOLUME{1}, 10:12, 'gray-Layer1',10);
VOLUME{1} = rmRunNumbersScriptLogMonotonic(VOLUME{1}, 13:14, 'gray-Layer1',10);

%Testing
VOLUME{1} = rmRunNumbersScriptLogMonotonic(VOLUME{1}, 14, 'gray-Layer1',10);

%Get resulting best fit HRFs from each model
hrfParamsLogMono1=viewGet(VOLUME{1},'rmhrf');
hrfParamsLogMono2=viewGet(VOLUME{1},'rmhrf');
hrfParamsLogMono3=viewGet(VOLUME{1},'rmhrf');
hrfParamsLogMono4=viewGet(VOLUME{1},'rmhrf');

%Average HRF parameters for different stimuli
hrfParamsLogMonoAverage{1}=hrfParamsLogMono1{1};
hrfParamsLogMonoAverage{2}=mean([hrfParamsLogMono1{2}; hrfParamsLogMono2{2}; hrfParamsLogMono3{2}; hrfParamsLogMono4{2}]);

times=0:1:24;
figure; plot(times, makeHRF(times, hrfParamsLogMonoAverage{2}), 'b');
hold on; plot(times, makeHRF(times, hrfParamsLogMonoAverage{2}), 'r');

figure; plot(times, makeHRF(times, hrfParamsLogMono1{2}), 'r');
hold on; plot(times, makeHRF(times, hrfParamsLogMono2{2}), 'g');
hold on; plot(times, makeHRF(times, hrfParamsLogMono3{2}), 'b');
hold on; plot(times, makeHRF(times, hrfParamsLogMono4{2}), 'c');

%Re-run models with common HRF
VOLUME{1} = rmRunNumbersScriptLogMonotonic(VOLUME{1}, 10:12, 'gray-Layer1',9, [], hrfParamsLogMonoAverage, sprintf('JointHRF-%s-Monotonic-Log-FullBlanks-DT0.5', datestr(now,'yyyymmdd-HHMMSS')));
VOLUME{1} = rmRunNumbersScriptLogMonotonic(VOLUME{1}, 13:14, 'gray-Layer1',9, [], hrfParamsLogMonoAverage, sprintf('JointHRF-%s-Monotonic-Log-FullBlanks-DT0.5', datestr(now,'yyyymmdd-HHMMSS')));

%Now load each model in turn, as model1, model2 etc,hrfParamsLogMono2=viewGet(VOLUME{1},'rmhrf'); and each tuned log2lin pRF
%model as modelTuned1, modelTuned2, etc
%Then load average monotonic model as model and params

%Then find all sites where tuned models fit better, in each individual configurations and set their VE to
%zero.
model1{1}.rss(model1{1}.rss>modelTuned1{1}.rss & modelTuned1{1}.x0>1.05)=model1{1}.rawrss(model1{1}.rss>modelTuned1{1}.rss & modelTuned1{1}.x0>1.05);
model2{1}.rss(model2{1}.rss>modelTuned2{1}.rss & modelTuned2{1}.x0>1.05)=model2{1}.rawrss(model2{1}.rss>modelTuned2{1}.rss & modelTuned2{1}.x0>1.05);
model3{1}.rss(model3{1}.rss>modelTuned3{1}.rss & modelTuned3{1}.x0>1.05)=model3{1}.rawrss(model3{1}.rss>modelTuned3{1}.rss & modelTuned3{1}.x0>1.05);
model4{1}.rss(model4{1}.rss>modelTuned4{1}.rss & modelTuned4{1}.x0>1.05)=model4{1}.rawrss(model4{1}.rss>modelTuned4{1}.rss & modelTuned4{1}.x0>1.05);

%Find minimum VE (conjunction analysis)
ve1= 1 - (model1{1}.rss ./ model1{1}.rawrss);
ve1(~isfinite(ve1)) = 0;
ve2= 1 - (model2{1}.rss ./ model2{1}.rawrss);
ve2(~isfinite(ve2)) = 0;
ve3= 1 - (model3{1}.rss ./ model3{1}.rawrss);
ve3(~isfinite(ve3)) = 0;
ve4= 1 - (model4{1}.rss ./ model4{1}.rawrss);
ve4(~isfinite(ve4)) = 0;

tmp=[ve1; ve2; ve3; ve4];
veAll=min(tmp,[],1);

model{1}.rss=(1-veAll).*model{1}.rawrss;

%To exclude tuned voxels based on mean model, less conservative than
%individual models
model{1}.rss(model{1}.rss>modelTuned5{1}.rss & modelTuned5{1}.x0>1.05)=model{1}.rawrss(model{1}.rss>modelTuned5{1}.rss & modelTuned5{1}.x0>1.05);

maxNumThreshold=6.5;
veTuned1= 1 - (modelTuned1{1}.rss ./ modelTuned1{1}.rawrss);
veTuned1(~isfinite(veTuned1)) = 0;
veTuned1(modelTuned1{1}.x0<maxNumThreshold)=0;

veTuned2= 1 - (modelTuned2{1}.rss ./ modelTuned2{1}.rawrss);
veTuned2(~isfinite(veTuned2)) = 0;
veTuned1(modelTuned2{1}.x0<maxNumThreshold)=0;

veTuned3= 1 - (modelTuned3{1}.rss ./ modelTuned3{1}.rawrss);
veTuned3(~isfinite(veTuned3)) = 0;
veTuned1(modelTuned3{1}.x0<maxNumThreshold)=0;

veTuned4= 1 - (modelTuned4{1}.rss ./ modelTuned4{1}.rawrss);
veTuned4(~isfinite(veTuned4)) = 0;
veTuned1(modelTuned4{1}.x0<maxNumThreshold)=0;

tmp=[veTuned1; veTuned2; veTuned3; veTuned4];
veTunedAll=min(tmp,[],1);
rssTunedAll=(1-veTunedAll).*model{1}.rawrss;

model{1}.x0(rssTunedAll<model{1}.rss)=ones(1, sum(rssTunedAll<model{1}.rss)).*3;
model{1}.rss(rssTunedAll<model{1}.rss)=rssTunedAll(rssTunedAll<model{1}.rss);
