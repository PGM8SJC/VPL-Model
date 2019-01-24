function [Stimulus] = generateStimulus()

% load stimulus parameters
load ./StimulusParam.mat;

screenSize = M.screenSize;
dotSize = M.dotSize;
dotDensity = M.dotDensity;
speed = M.speed;
direction = M.direction;
apertureDiam = M.apertureDiam;
apertureLoc = M.apertureLoc;
coherence = M.coherence;
lifeTime = M.lifeTime;
duration = M.duration;
framerate = M.framerate;


numDots = floor(apertureDiam * apertureDiam * dotDensity); % assuming that the aperture is square

dotsInitialPosition = [randperm(apertureDiam,numDots) + apertureLoc(1) - apertureDiam/2;randperm(apertureDiam,numDots)+ apertureLoc(2) - apertureDiam/2];
initialLife = rand(1,numDots) * lifeTime;
initialMotionVector = repmat([speed;direction],1,numDots);
numFramesNeeded = floor(duration * framerate / 1000);

signalDotsIdx = rand(1,numDots) < coherence;
noiseDotsIdx = ~signalDotsIdx;

% Translational Motion
noiseDirection = rand(1,numDots) * 2*pi;

for framecount = 1:numFramesNeeded
    
    % move the dots 
    if framecount == 1
        dotsPosition = dotsInitialPosition;
        dotsLife = initialLife;
        motionVectors = initialMotionVector;
    else
        dotsPosition(1,signalDotsIdx,framecount) = floor(dotsPosition(1,signalDotsIdx,framecount-1) + (1/framerate) * speed * cos(direction));
        dotsPosition(2,signalDotsIdx,framecount) = floor(dotsPosition(2,signalDotsIdx,framecount-1) + (1/framerate) * speed * sin(direction));
        
        dotsPosition(1,noiseDotsIdx,framecount) = floor(dotsPosition(1,noiseDotsIdx,framecount-1) + (1/framerate) * speed * cos(noiseDirection(noiseDotsIdx)));
        dotsPosition(2,noiseDotsIdx,framecount) = floor(dotsPosition(2,noiseDotsIdx,framecount-1) + (1/framerate) * speed * sin(noiseDirection(noiseDotsIdx)));
        
        dotsLife = dotsLife + (1/framerate) * 1000;
        motionVectors(1,:,framecount) = repmat(speed,1,numDots);
        motionVectors(2,signalDotsIdx,framecount) = repmat(direction,1,sum(signalDotsIdx));
        motionVectors(2,noiseDotsIdx,framecount) = noiseDirection(noiseDotsIdx);
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


% Complex Motion


% show the simulated stimulus
% for i = 1:numFramesNeeded
% %     figure;
%     I = zeros(screenSize(1),screenSize(1));
%     idx = sub2ind(size(I), [squeeze(dotsPosition(1,:,i))], [squeeze(dotsPosition(2,:,i))]);
%     I(idx) = 1;
%     figure(1);imagesc(I);colormap(gray);pause(1/framerate)
% end

for i = 1:numFramesNeeded
    figure(1);imagesc(zeros(screenSize(1),screenSize(1)));colormap(gray);hold on;
    figure(1);quiver(squeeze(dotsPosition(2,:,i)),squeeze(dotsPosition(1,:,i)),motionVectors(1,:,i) .* sin(motionVectors(2,:,i)),motionVectors(1,:,i) .* cos(motionVectors(2,:,i)),'Color','k');
    pause(1/framerate)
end

Stimulus.dotsPosition = dotsPosition;
Stimulus.motionVectors = motionVectors;


end