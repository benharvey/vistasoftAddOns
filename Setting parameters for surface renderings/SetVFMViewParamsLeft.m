VOLUME{1} = resetDisplayModes(VOLUME{1}); VOLUME{1} = setPhWindow(VOLUME{1}, [0 2*pi]);
VOLUME{1} = setClipMode(VOLUME{1}, 'map', [0 6.35]);
%VOLUME{1}.ui.mapWinMax.sliderHandle.Value=8;
% VOLUME{1} = setClipMode(VOLUME{1}, 'map', [0.3 1]);
% VOLUME{1} = setClipMode(VOLUME{1}, 'amp', [0.3 1]);
% VOLUME{1} = setClipMode(VOLUME{1}, 'map', [0 4.1]);
% VOLUME{1} = setClipMode(VOLUME{1}, 'amp', [0 4.1]);

VOLUME{1}.ui.mapMode=setColormap(VOLUME{1}.ui.mapMode, 'hsvTbCmap');
VOLUME{1} = cmapRedgreenblue(VOLUME{1}, 'ph', 1);
VOLUME{1}.ui.displayMode='ph';
%VOLUME{1}=rotateCmap(VOLUME{1},  180);

%RHS
%VOLUME{1}=flipCmap(VOLUME{1});

VOLUME{1}=rotateCmap(VOLUME{1},  270);


%VOLUME{1} = setClipMode(VOLUME{1}, 'amp', [0 2]);
%VOLUME{1}=linearizeCmap(VOLUME{1});
VOLUME{1} = refreshScreen(VOLUME{1}, 1);