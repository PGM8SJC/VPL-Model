function [MST, params] = MSTbank(MT, MTparams)

MT = X;
MTparams = params;
% This function simulates MST population from MT population. Every MST unit
% gets input from N MT units. The pattern of connectivity is simulated in
% this code. 

spatialWidth = sqrt(size(MT,3));
Thetas = unique(MTparams.thetas);
numThetas = length(Thetas);
Speeds = unique(MTparams.speeds);
numSpeeds = length(Speeds);

% MST bank parameters
spatialSigma = 5;
numMSTUnitsPerPos = 100;

% MT units at which locations send input to each MST unit? Model: a
% Gaussian centered at each MST location sets the probability of the
% location of projecting MT units. For each locati                                                                                                                           on, there are
% numMSTUnitsPerPos of MST units.
whichSptialPosition = zeros(numMSTUnitsPerPos,spatialWidth*spatialWidth);
for i = 1:spatialWidth
    for j = 1:spatialWidth
        
        posIdx = sub2ind([spatialWidth,spatialWidth],i,j);
        posTemp = round(mvnrnd([i,j],[spatialSigma,0;0,spatialSigma],numMSTUnitsPerPos));
        [ro,co] = find(posTemp <= 0 | posTemp > spatialWidth);
        posTemp(ro,1) = i;
        posTemp(ro,2) = j;
        
        whichPosIdx = sub2ind([spatialWidth,spatialWidth],posTemp(:,1),posTemp(:,2));
        whichSptialPosition(:,posIdx) = whichPosIdx;
        clear whichPosIdx posTemp;

        
    end
end





end