function estimateMotionField(Stimulus)

load ./StimulusParam.mat;

screenSize = M.screenSize;

motionFieldResRatio = 1/5; % ratio to the screen size e.g. motionFieldResRatio = 1/5 means that every 5 pixels will have one motion vector
motionFieldSize = floor(screenSize .* motionFieldResRatio);
dotsPosition = Stimulus.dotsPosition;
motionVectors = Stimulus.motionVectors;

k = 0;
for i = 1:motionFieldSize(1)
    for j = 1:motionFieldSize(2)
        k = k + 1;
        lowLimitX = (i - 1) * (1/motionFieldResRatio) + 1;
        upLimitX = (i) * (1/motionFieldResRatio);
        motionFieldSizeLim(1,k) = lowLimitX;
        motionFieldSizeLim(2,k) = upLimitX;
        
        lowLimitY = (j - 1) * (1/motionFieldResRatio) + 1;
        upLimitY = (j) * (1/motionFieldResRatio);
        motionFieldSizeLim(3,k) = lowLimitY;
        motionFieldSizeLim(4,k) = upLimitY;
        
        whichDotsInsideX = dotsPosition(1,:,:) >= lowLimitX & dotsPosition(1,:,:) & ...
            dotsPosition(2,:,:) >= lowLimitY & dotsPosition(2,:,:) <= upLimitY;
        whichDotsInsideY = whichDotsInsideX;
        whichDotsInside = cat(1,whichDotsInsideX,whichDotsInsideY);
        
        allMotionVectorsInside = motionVectors(whichDotsInside);
        if isempty(allMotionVectorsInside)
            MotionField(i,j,1) = 0;
            MotionField(i,j,2) = 0;
        else
        end
        
        
    end
end


end