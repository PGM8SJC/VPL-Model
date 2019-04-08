function setStimulusParameters()

M.screenSize      = [600,600];    % pixels    
M.dotSize         = [1];          % pixels
M.dotDensity      = [0.01];       % number of dots per pixel^2
M.speed           = [100];        % pixel/sec
M.direction       = [pi];          % radian
M.apertureDiam    = [100];
M.apertureLoc     = [250,250];
M.coherence       = [1];          % portion of dots moving together
M.lifeTime        = [500];        % ms
M.duration        = [500];        % ms
M.framerate       = 60;
M.motiontype      = 'complex';
M.omega0          = 2;

save('./StimulusParam.mat','M');


end