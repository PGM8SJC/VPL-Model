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
motiontype = M.motiontype;
omega0 = M.omega0;

numDots = floor(apertureDiam * apertureDiam * dotDensity); % assuming that the aperture is square

% dotsInitialPosition = [randperm(apertureDiam,numDots) + apertureLoc(1) - apertureDiam/2;randperm(apertureDiam,numDots)+ apertureLoc(2) - apertureDiam/2];
dotsInitialPosition = [randi(apertureDiam,1,numDots) + apertureLoc(1) - apertureDiam/2;randi(apertureDiam,1,numDots)+ apertureLoc(2) - apertureDiam/2];
initialLife = rand(1,numDots) * lifeTime;
initialMotionVector = repmat([speed;direction],1,numDots);
numFramesNeeded = floor(duration * framerate / 1000);

signalDotsIdx = rand(1,numDots) < coherence;
noiseDotsIdx = ~signalDotsIdx;

% Translational Motion
if strcmp(motiontype,'simple')
    
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
end

% Complex Motion
if strcmp(motiontype ,'complex')
noiseDirection = rand(1,numDots) * 2*pi;
for framecount = 1:numFramesNeeded
    
    % move the dots 
    if framecount == 1
        dotsPosition = dotsInitialPosition;
        dotsLife = initialLife;
        motionVectors = initialMotionVector;
    else
        u = omega0 .* ((dotsPosition(1,signalDotsIdx,framecount-1) - apertureLoc(1)) .* cos(direction) - ...
                      (dotsPosition(2,signalDotsIdx,framecount-1) - apertureLoc(2)) .* sin(direction)); 
        v = omega0 .* ((dotsPosition(1,signalDotsIdx,framecount-1) - apertureLoc(1)) .* sin(direction) + ... 
                       (dotsPosition(2,signalDotsIdx,framecount-1) - apertureLoc(2)) .* cos(direction));
        dotsPosition(1,signalDotsIdx,framecount) = floor(dotsPosition(1,signalDotsIdx,framecount-1) + (1/framerate) * u);
        dotsPosition(2,signalDotsIdx,framecount) = floor(dotsPosition(2,signalDotsIdx,framecount-1) + (1/framerate) * v);
        
        unoise = omega0 .* ((dotsPosition(1,noiseDotsIdx,framecount-1) - apertureLoc(1)) .* cos(noiseDirection(noiseDotsIdx)) - ...
                      (dotsPosition(2,noiseDotsIdx,framecount-1) - apertureLoc(2)) .* sin(noiseDirection(noiseDotsIdx))); 
        vnoise = omega0 .* ((dotsPosition(1,noiseDotsIdx,framecount-1) - apertureLoc(1)) .* sin(noiseDirection(noiseDotsIdx)) + ... 
                       (dotsPosition(2,noiseDotsIdx,framecount-1) - apertureLoc(2)) .* cos(noiseDirection(noiseDotsIdx)));
        dotsPosition(1,noiseDotsIdx,framecount) = floor(dotsPosition(1,noiseDotsIdx,framecount-1) + (1/framerate) * unoise);
        dotsPosition(2,noiseDotsIdx,framecount) = floor(dotsPosition(2,noiseDotsIdx,framecount-1) + (1/framerate) * vnoise);
        
        dotsLife = dotsLife + (1/framerate) * 1000;
        motionVectors(1,signalDotsIdx,framecount) = sqrt(u.^2 + v.^2);
        motionVectors(1,noiseDotsIdx,framecount) = sqrt(unoise.^2 + vnoise.^2);
        motionVectors(2,signalDotsIdx,framecount) = atan2(v,u) .* (atan2(v,u) >= 0) +  (atan2(v,u) + 2*pi) .* (atan2(v,u) < 0);
        motionVectors(2,noiseDotsIdx,framecount) = atan2(vnoise,unoise) .* (atan2(vnoise,unoise) >= 0) +  (atan2(vnoise,unoise) + 2*pi) .* (atan2(vnoise,unoise) < 0);
%         motionVectors(2,signalDotsIdx,framecount) = atan(v./u);
%         motionVectors(2,noiseDotsIdx,framecount) = atan(vnoise./unoise);
%         motionVectors(1,:,framecount) = repmat(speed,1,numDots);
%         motionVectors(2,signalDotsIdx,framecount) = repmat(direction,1,sum(signalDotsIdx));
%         motionVectors(2,noiseDotsIdx,framecount) = noiseDirection(noiseDotsIdx);
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

end

% show the simulated stimulus
% for i = 1:numFramesNeeded
% %     figure;
%     I = zeros(screenSize(1),screenSize(1));
%     idx = sub2ind(size(I), [squeeze(dotsPosition(1,:,i))], [squeeze(dotsPosition(2,:,i))]);
%     I(idx) = 1;
%     figure(1);imagesc(I);colormap(gray);pause(1/framerate)
% end

% for i = 1:numFramesNeeded
%     figure(1);imagesc(zeros(screenSize(1),screenSize(1)));colormap(gray);hold on;
%     figure(1);quiver(squeeze(dotsPosition(2,:,i)),squeeze(dotsPosition(1,:,i)),motionVectors(1,:,i) .* sin(motionVectors(2,:,i)),motionVectors(1,:,i) .* cos(motionVectors(2,:,i)),'Color','k');
%     pause(1/framerate)
% end

Stimulus.dotsPosition = dotsPosition;
Stimulus.motionVectors = motionVectors;


end