function setStimulusParameters()

M.screenSize      = [500,500];    % pixels    
M.dotSize         = [1];          % pixels
M.dotDensity      = [0.01];       % number of dots per pixel^2
M.speed           = [100];          % pixel/sec
M.direction       = [pi/2];          % radian
M.apertureDiam    = [100];
M.apertureLoc     = [350,150];
M.coherence       = [.8];          % portion of dots moving together
M.lifeTime        = [200];        % ms
M.duration        = [500];        % ms
M.framerate        = 60;

save('./StimulusParam.mat','M');


end