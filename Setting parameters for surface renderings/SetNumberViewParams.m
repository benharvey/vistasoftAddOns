VOLUME{1} = rmLoad(VOLUME{1}, 1, 'x0', 'map');

VOLUME{1} = setClipMode(VOLUME{1}, 'map', [1 8.5]);
VOLUME{1} = refreshScreen(VOLUME{1}, 1);
VOLUME{1}.ui.mapMode=setColormap(VOLUME{1}.ui.mapMode, 'hsvCmap'); 
VOLUME{1}=refreshScreen(VOLUME{1}, 1);
% VOLUME{1} = setClipMode(VOLUME{1}, 'amp', [0 20]);
% VOLUME{1} = refreshScreen(VOLUME{1}, 1);
VOLUME{1}.ui.mapWinMin.sliderHandle.Value=1.05;
VOLUME{1}.ui.mapWinMax.sliderHandle.Value=7.00;
VOLUME{1}=refreshScreen(VOLUME{1}, 1);