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
spatialSigma = 2;
numMSTUnitsPerPos = 100;
numMTsubUnitsPerMSTUnit = 10;

% MT units at which locations send input to each MST unit? 
% My model here: 
% a Gaussian centered at each MST location sets the probability of the
% location of projecting MT units (numMTsubUnitsPerMSTUnit). For each location                                                                                                                           on, there are
% numMSTUnitsPerPos of MST units are simulated.
whichSpatialPosition = zeros(numMSTUnitsPerPos,numMTsubUnitsPerMSTUnit,spatialWidth*spatialWidth);
for mstcounter = 1:numMSTUnitsPerPos
    for i = 1:spatialWidth
        for j = 1:spatialWidth
            
            posIdx = sub2ind([spatialWidth,spatialWidth],i,j);
            posTemp = round(mvnrnd([i,j],[spatialSigma,0;0,spatialSigma],numMTsubUnitsPerMSTUnit));
            [ro,co] = find(posTemp <= 0 | posTemp > spatialWidth);
            posTemp(ro,1) = i;
            posTemp(ro,2) = j;
            
            whichPosIdx = sub2ind([spatialWidth,spatialWidth],posTemp(:,1),posTemp(:,2));
            whichSpatialPosition(mstcounter,:,posIdx) = whichPosIdx;
            clear whichPosIdx posTemp;
            
            
        end
    end
end

% set the tuning parameters of the MT subunits for each MST unit;
% Note: currently, the preferred direction and speed of each subunit is
% randomly sampled from the range of thetas and speeds. But, this is
% obviously not optimal.
whichPreferredSpeed = randi(numSpeeds,numMSTUnitsPerPos,numMTsubUnitsPerMSTUnit,spatialWidth^2);
whichPreferredTheta = randi(numThetas,numMSTUnitsPerPos,numMTsubUnitsPerMSTUnit,spatialWidth^2);

% adding all the selected MT subunits in one matrix
% each element in allMTresp is an MT response to the stimulus
allMTresp = zeros(numMSTUnitsPerPos,numMTsubUnitsPerMSTUnit,spatialWidth^2);
for mstcounter = 1:numMSTUnitsPerPos
    for spatialIdx = 1:spatialWidth^2
        thisMSTsubunits = [squeeze(whichSpatialPosition(mstcounter,:,spatialIdx));squeeze(whichPreferredSpeed(mstcounter,:,1)); squeeze(whichPreferredTheta(mstcounter,:,1))];
        for subunitcounter = 1:numMTsubUnitsPerMSTUnit
            thisMT = MT(thisMSTsubunits(2,subunitcounter),thisMSTsubunits(3,subunitcounter),thisMSTsubunits(1,subunitcounter));
            allMTresp(mstcounter,subunitcounter,spatialIdx) = thisMT;
            clear thisMT
            
        end
    end
end


% assigning the synaptic weights from MT subunits to MST units 
MT2MSTweights = rand(numMSTUnitsPerPos,numMTsubUnitsPerMSTUnit,spatialWidth^2) - 0.2; % 80% of weights are positive


% summation over MT subunits to generate MST units responses
allMTrespSum = squeeze(sum(MT2MSTweights .* allMTresp,2));
postMax = max(allMTrespSum,0);





end