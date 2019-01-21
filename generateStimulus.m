function [dotsPosition,motionVectors] = generateStimulus()

ScreenSize      = [500,500];    % pixels    
dotSize         = [1];          % pixels
dotDensity      = [0.01];       % number of dots per pixel^2
speed           = [100];          % pixel/sec
direction       = [pi];          % radian
apertureDiam    = [100];
apertureLoc     = [350,150];
coherence       = [1];          % portion of dots moving together
lifeTime        = [500];        % ms
duration        = [500];        % ms
framerate        = 60;

numDots = floor(apertureDiam * apertureDiam * dotDensity); % assuming that the aperture is square

dotsInitialPosition = [randperm(apertureDiam,numDots) + apertureLoc(1) - apertureDiam/2;randperm(apertureDiam,numDots)+ apertureLoc(2) - apertureDiam/2];
initialLife = rand(1,numDots) * lifeTime;
initialMotionVector = repmat([speed;direction],1,numDots);
numFramesNeeded = floor(duration * framerate / 1000);

for framecount = 1:numFramesNeeded
    
    % move the dots 
    if framecount == 1
        dotsPosition = dotsInitialPosition;
        dotsLife = initialLife;
        motionVectors = initialMotionVector;
    else
        dotsPosition(1,:,framecount) = floor(dotsPosition(1,:,framecount-1) + (1/framerate) * speed * cos(direction));
        dotsPosition(2,:,framecount) = floor(dotsPosition(2,:,framecount-1) + (1/framerate) * speed * sin(direction));
        dotsLife = dotsLife + (1/framerate) * 1000;
        motionVectors(:,:,framecount) = repmat([speed;direction],1,numDots);
    end
    
    % relocate the dead dots
    deadDotsIdx = dotsLife > lifeTime;
    numDeadDots = sum(deadDotsIdx);
    dotsPosition(:,deadDotsIdx,framecount) = [randperm(apertureDiam,numDeadDots) + apertureLoc(1) - apertureDiam/2;randperm(apertureDiam,numDeadDots)+ apertureLoc(2) - apertureDiam/2];
    dotsLife(deadDotsIdx) = 0;
    
    % set the motion vector for deadDots to zero
    if numDeadDots > 0
        motionVectors(:,deadDotsIdx,framecount) = repmat([0;direction],1,numDeadDots);
    end
    
    % replace the dots at the aperture border
    currentDotsPosition = squeeze(dotsPosition(:,:,end));
    upborderDotsIdx = currentDotsPosition(2,:) > (apertureLoc(2) + apertureDiam/2);
    lowborderDotsIdx = currentDotsPosition(2,:) < (apertureLoc(2) - apertureDiam/2);
    rightborderDotsIdx = currentDotsPosition(1,:) > (apertureLoc(1) + apertureDiam/2);
    leftborderDotsIdx = currentDotsPosition(1,:) < (apertureLoc(1) - apertureDiam/2);
    
    borderDotsIdx = upborderDotsIdx | lowborderDotsIdx | rightborderDotsIdx | leftborderDotsIdx;
    numBorderDots = sum(borderDotsIdx);
    dotsPosition(:,borderDotsIdx,framecount) = [randperm(apertureDiam,numBorderDots) + apertureLoc(1) - apertureDiam/2;randperm(apertureDiam,numBorderDots)+ apertureLoc(2) - apertureDiam/2];
    dotsLife(borderDotsIdx) = 0;
     
    
end

% show the simulated stimulus
% for i = 1:numFramesNeeded
% %     figure;
%     I = zeros(500,500);
%     idx = sub2ind(size(I), [squeeze(dotsPosition(1,:,i))], [squeeze(dotsPosition(2,:,i))]);
%     I(idx) = 1;
%     figure(1);imagesc(I);colormap(gray);pause(1/framerate)
% end

for i = 1:numFramesNeeded
    figure(1);imagesc(zeros(500,500));colormap(gray);hold on;
    figure(1);quiver(squeeze(dotsPosition(1,:,i)),squeeze(dotsPosition(2,:,i)),motionVectors(1,:,i) .* cos(motionVectors(2,:,i)),motionVectors(1,:,i) .* sin(motionVectors(2,:,i)),'Color','k');
    pause(1/framerate)
end


end