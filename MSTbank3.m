function [MST, params] = MSTbank3(MT, MTparams,newMST)

% Main difference compared to MSTbank: ... It involves convolution. 


global spatialSigma numSpeeds numThetas numMSTUnitsPerPos numMTsubUnitsPerMSTUnit spatialWidth
if nargin < 3
    newMST = false;
end

% MT = MT;
% MTparams = params;
% newMST = true;

% This function simulates MST population from MT population. Every MST unit
% gets input from N MT units. The pattern of connectivity is simulated in
% this code. Every MST unit has a receptive field (a kernel with specific
% connections from MT). The two-dimentsional kernel is convolved with the
% maps that are the output of MT layer.

spatialWidth = sqrt(size(MT,3));
Thetas = unique(MTparams.thetas);
numThetas = length(Thetas);
Speeds = unique(MTparams.speeds);
numSpeeds = length(Speeds);

% MST bank parameters
spatialSigma = 5;%11;%2;
numMSTUnitsPerPos = 20;%40;
numMTsubUnitsPerMSTUnit = 5;%
kernelW = 5;%11;
strideSize = 4;


% MT units at which locations send input to each MST unit? 
% My model here: 
% There is a kernel that corresponds to a single MST receptive field. The
% non-zerp elements of the kernel (chosen randomly) determine the MT-to-MST
% projections.


% set the tuning parameters of the MT subunits for each MST unit;
% Note: currently, the preferred direction and speed of each subunit is
% randomly sampled from the range of thetas and speeds. But, this is
% obviously not optimal. In order to be able to resimulate the same
% population, the random connectivity matrix (generated by
% setMTMSTConnectivity() is saved and can be reloaded.
if newMST
    Connectivity = setMTMSTConnectivity();
    whichPreferredSpeed = Connectivity.whichPreferredSpeed;
    whichPreferredTheta = Connectivity.whichPreferredTheta;
    whichSpatialPosition = Connectivity.whichSpatialPosition;
    MT2MSTweights = Connectivity.MT2MSTweights; % assigning the synaptic weights from MT subunits to MST units 

else
    load Connectivity.mat;
    whichPreferredSpeed = Connectivity.whichPreferredSpeed;
    whichPreferredTheta = Connectivity.whichPreferredTheta;
    whichSpatialPosition = Connectivity.whichSpatialPosition;
    MT2MSTweights = Connectivity.MT2MSTweights; % assigning the synaptic weights from MT subunits to MST units 
end



for mstcounter = 1:numMSTUnitsPerPos
    thisWeight = MT2MSTweights(mstcounter,:);
    thisMTsubunits = [whichPreferredSpeed(mstcounter,:);whichPreferredTheta(mstcounter,:)];
    thisSpatialPositions = whichSpatialPosition(mstcounter,:);
    [I,J] = ind2sub([spatialSigma,spatialSigma],thisSpatialPositions);
    
    % build the kernel
    Kernel = zeros(numSpeeds,numThetas,spatialSigma,spatialSigma);
    
    for subunitcounter = 1:numMTsubUnitsPerMSTUnit
        Kernel(thisMTsubunits(1,subunitcounter),thisMTsubunits(2,subunitcounter),I(subunitcounter),J(subunitcounter)) = thisWeight(subunitcounter);
    end
    
    MTreshape = reshape(MT,numSpeeds,numThetas,spatialWidth,spatialWidth);
    
    thisMST(:,:,mstcounter) = simpleNdConv(MTreshape,Kernel);
    
    
end

pooledMST = pooling(thisMST, kernelW, strideSize);
downsampledSpatialWidth = size(pooledMST,1);
% pooledMST = max(pooledMST,0).^0.1;
MST = reshape(pooledMST,downsampledSpatialWidth^2,numMSTUnitsPerPos)';
params.Connectivity = Connectivity;





end

function Connectivity = setMTMSTConnectivity()

global spatialSigma numSpeeds  numThetas numMSTUnitsPerPos numMTsubUnitsPerMSTUnit spatialWidth

for mstcounter = 1:numMSTUnitsPerPos
    whichSpatialPosition(mstcounter,:) = randperm(spatialSigma.^2,numMTsubUnitsPerMSTUnit);
end


whichPreferredSpeed = randi(numSpeeds,numMSTUnitsPerPos,numMTsubUnitsPerMSTUnit);
whichPreferredTheta = randi(numThetas,numMSTUnitsPerPos,numMTsubUnitsPerMSTUnit);
MT2MSTweights = rand(numMSTUnitsPerPos,numMTsubUnitsPerMSTUnit) - 0.2; % 80% of weights are positive
NegWeightIdx = find(MT2MSTweights<0);
MT2MSTweights(NegWeightIdx) = -rand(1,length(NegWeightIdx));
PosWeightIdx = find(MT2MSTweights>=0);
MT2MSTweights(PosWeightIdx) = rand(1,length(PosWeightIdx));

Connectivity.whichSpatialPosition = whichSpatialPosition;
Connectivity.whichPreferredSpeed = whichPreferredSpeed;
Connectivity.whichPreferredTheta = whichPreferredTheta;
Connectivity.MT2MSTweights = MT2MSTweights;

save('Connectivity.mat','Connectivity');

end


function Y = simpleNdConv(X,K)

% this function convolves the N-dimensional array X with M-dimensional
% kernel K, only over the two first dimensions. It is a 2d convolution over
% an Nd array

[d1,d2,dx,dy] = size(X);
[dk1,dk2,w,w] = size(K);

if (d1 ~= dk1) || (d2 ~= dk2)
    error('kernel and the array don not have the same speed and/or theta dimensions');
end

padsize = ceil(w./2);
padX = padarray(X,[0,0,padsize,padsize]);

for i = 1:dx
    
    for j = 1:dy
        
        x = padX(:,:,(i + padsize-floor(w/2)): (i + padsize+floor(w/2)),(j + padsize-floor(w/2)): (j + padsize+floor(w/2)));
        x_nonlinear = max(x,0).^.08;
        y = x_nonlinear .* K;
        
        Y(i,j) = sum(y(:));
        
    end
end



end


function pooledmap = pooling(maps, kernelW, strideSize)

global spatialWidth
padsize = ceil(kernelW./2);
maps = padarray(maps,[padsize,padsize,0]);

k = 0;
for i = 1:strideSize:(spatialWidth + padsize*2 - kernelW)
    k = k + 1;
    h = 0;
    for j = 1:strideSize:(spatialWidth + padsize*2 - kernelW)
        h = h + 1;
        thispart = maps(i:(i+kernelW-1),j:(j+kernelW-1),:);
%         pooledmap(k,h,:) = max(squeeze(max((thispart),[],1)),[],1);
%          pooledmap(k,h,:) = mean(squeeze(mean(thispart,1)),1);
%           pooledmap(k,h,:) = median(rsqueeze(median(thispart,1)),1);
        MAXpooled = max(squeeze(max((thispart),[],1)),[],1);
        MINpooled = min(squeeze(min((thispart),[],1)),[],1);
        MAXorMIN = abs(MAXpooled) > abs(MINpooled);
        pooledmap(k,h,MAXorMIN) = MAXpooled(MAXorMIN);
        pooledmap(k,h,~MAXorMIN) = MINpooled(~MAXorMIN);
    end
end

end