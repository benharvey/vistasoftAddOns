VOLUME{1} = setClipMode(VOLUME{1}, 'map', [0 1.6]);
VOLUME{1} = refreshScreen(VOLUME{1}, 1);
VOLUME{1}.ui.mapMode=setColormap(VOLUME{1}.ui.mapMode, 'hsvCmap'); 
VOLUME{1}=refreshScreen(VOLUME{1}, 1);
VOLUME{1} = setClipMode(VOLUME{1}, 'amp', [0 1.6]);
%VOLUME{1}=linearizeCmap(VOLUME{1});
VOLUME{1} = refreshScreen(VOLUME{1}, 1);