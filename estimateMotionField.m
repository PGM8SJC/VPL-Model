function MotionField = estimateMotionField(Stimulus)

% To calculate the motion vector on each grid point:
% - All the dots inside each grid point are extracted from the stimulus
% - All the dots are embedded from all the stimulus frames
% - The velocity vectors are averaged.
% - The same process is repeated for all the points on the grid of motion
% field.
% - The motion vector from different frames could be dealt with
% differently! 

load ./StimulusParam.mat;

screenSize = M.screenSize;

motionFieldResRatio = 1/10; % ratio to the screen size e.g. motionFieldResRatio = 1/5 means that every 5 pixels will have one motion vector
motionFieldSize = floor(screenSize .* motionFieldResRatio);
dotsPosition = Stimulus.dotsPosition;
motionVectors = Stimulus.motionVectors;

for i = 1:motionFieldSize(1)
    for j = 1:motionFieldSize(2)
        lowLimitX = (i - 1) * (1/motionFieldResRatio) + 1;
        upLimitX = (i) * (1/motionFieldResRatio);

        lowLimitY = (j - 1) * (1/motionFieldResRatio) + 1;
        upLimitY = (j) * (1/motionFieldResRatio);

        whichDotsInsideX = dotsPosition(1,:,:) >= lowLimitX & dotsPosition(1,:,:) <= upLimitX & ...
            dotsPosition(2,:,:) >= lowLimitY & dotsPosition(2,:,:) <= upLimitY;
        whichDotsInsideY = whichDotsInsideX;
        whichDotsInside = cat(1,whichDotsInsideX,whichDotsInsideY);
        
        
        allMotionVectorsInside = motionVectors(whichDotsInside);
        if isempty(allMotionVectorsInside)
            averageMotionInsideAmp(i,j) = 0;
            averageMotionInsideAngle(i,j) = 0;
            
        elseif ~isempty(allMotionVectorsInside)
            
            allMotionVectorsInsideTwoRows = reshape(allMotionVectorsInside,2,length(allMotionVectorsInside)./2);
            
            averageMotionInsideX = sum(allMotionVectorsInsideTwoRows(1,:) .* cos(allMotionVectorsInsideTwoRows(2,:)));
            averageMotionInsideY = sum(allMotionVectorsInsideTwoRows(1,:) .* sin(allMotionVectorsInsideTwoRows(2,:)));
            averageMotionInsideAmp(i,j) = sqrt(averageMotionInsideX^2 + averageMotionInsideY^2)./(length(allMotionVectorsInside)./2);
            averageMotionInsideAngle(i,j) = atan(averageMotionInsideY./averageMotionInsideX);
            if averageMotionInsideX == 0 && averageMotionInsideY == 0
                averageMotionInsideAngle(i,j) = 0;
            end
            if (averageMotionInsideX < 0 && averageMotionInsideY < 0) || (averageMotionInsideX < 0 && averageMotionInsideY > 0)
                averageMotionInsideAngle(i,j) = averageMotionInsideAngle(i,j) + pi;
            end
            
        end
        
       
    end
end
% figure;imagesc(zeros(motionFieldSize(2),motionFieldSize(1)));colormap(gray);hold on;
% quiver(1:motionFieldSize(2),1:motionFieldSize(1),averageMotionInsideAmp .* sin(averageMotionInsideAngle),averageMotionInsideAmp .* cos(averageMotionInsideAngle),'k');
MotionField.averageMotionInsideAmp = averageMotionInsideAmp;
MotionField.averageMotionInsideAngle = averageMotionInsideAngle;
end