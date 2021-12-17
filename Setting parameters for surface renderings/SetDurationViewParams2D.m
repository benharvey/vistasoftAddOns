VOLUME{1} = rmLoad(VOLUME{1}, 1, 'y0', 'amp');
VOLUME{1} = rmLoad(VOLUME{1}, 1, 'x0', 'map');

VOLUME{1} = setClipMode(VOLUME{1}, 'amp', [0 1.1]);
VOLUME{1} = setClipMode(VOLUME{1}, 'map', [0 1.1]);

VOLUME{1} = refreshScreen(VOLUME{1}, 1);
VOLUME{1}.ui.ampMode=setColormap(VOLUME{1}.ui.ampMode, 'hsvCmap'); 
VOLUME{1}.ui.mapMode=setColormap(VOLUME{1}.ui.mapMode, 'hsvCmap'); 
VOLUME{1}.ui.mapWinMin.sliderHandle.Value=0.05;
VOLUME{1}.ui.mapWinMax.sliderHandle.Value=2.00;
VOLUME{1}=refreshScreen(VOLUME{1}, 1);
%VOLUME{1} = setClipMode(VOLUME{1}, 'amp', [0 2]);
%VOLUME{1}=linearizeCmap(VOLUME{1});
VOLUME{1} = refreshScreen(VOLUME{1}, 1);