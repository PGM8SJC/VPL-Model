%% simulate MST population responses from one stimulus
[Stimulus] = generateStimulus();
MotionField = estimateMotionField(Stimulus);

Xr = reshape(log(MotionField.averageMotionInsideAmp/10),1,size(MotionField.averageMotionInsideAmp,1)*size(MotionField.averageMotionInsideAmp,2));
Xt = reshape(MotionField.averageMotionInsideAngle,1,size(MotionField.averageMotionInsideAngle,1)*size(MotionField.averageMotionInsideAngle,2));
maxspeed = single(max(Xr(:)));
[MT, params] = MTbank(Xt,Xr,maxspeed);
MT = squeeze(MT);
newMST = true;
tic;
[MST, params] = MSTbank3(MT, params,newMST);     
toc;


%% simulate responses of a MST population to a range of stimuli motions
load ./StimulusParam.mat;
i = 0;
for dir = 0:pi/4:(2*pi - pi/4)
    i = i + 1
    M.direction = dir;
    save('./StimulusParam.mat','M');
    
    [Stimulus] = generateStimulus();
    MotionField = estimateMotionField(Stimulus);
    
    Xr = reshape(log(MotionField.averageMotionInsideAmp/10),1,size(MotionField.averageMotionInsideAmp,1)*size(MotionField.averageMotionInsideAmp,2));
    Xt = reshape(MotionField.averageMotionInsideAngle,1,size(MotionField.averageMotionInsideAngle,1)*size(MotionField.averageMotionInsideAngle,2));
    maxspeed = single(max(Xr(:)));
    [MT, MTparams] = MTbank(Xt,Xr,maxspeed);
    MT = squeeze(MT);
    newMST = false;
    [MST, MSTparams] = MSTbank3(MT, MTparams,newMST);
%     [X1] = max(MST,[],2);
    X1 = MST(:,36);
    NeuroResp(i,:) = X1(:);
%     close;
end

figure;
 
for i = 1:50
  baseline = min(NeuroResp(:,i));
themax = max(NeuroResp(:,i));
rg = [baseline-(themax-baseline),themax];
    subplot(7,8,i);colormap(jet);polarmosaic(NeuroResp(:,i),rg,.35,1);box off
end

%% simulate responses of a MST population to different stimuli locations
tic;
load ./StimulusParam.mat;
dircounter = 0;
allApertureLoc = [100 100; 300 100; 500 100; 100 300; 300 300; 500 300; 100 500; 300 500; 500 500];  
for dir = 0:pi/4:(2*pi - pi/4)
    dircounter = dircounter + 1;
    M.direction = dir;
for loccounter = 1:9
    fprintf(['location ',num2str(loccounter), ' dir ',num2str(dircounter),'\n']);
    M.apertureLoc = allApertureLoc(loccounter,:);
    save('./StimulusParam.mat','M');
    
    [Stimulus] = generateStimulus();
    MotionField = estimateMotionField(Stimulus); 
    
    Xr = reshape(log(MotionField.averageMotionInsideAmp/10),1,size(MotionField.averageMotionInsideAmp,1)*size(MotionField.averageMotionInsideAmp,2));
    Xt = reshape(MotionField.averageMotionInsideAngle,1,size(MotionField.averageMotionInsideAngle,1)*size(MotionField.averageMotionInsideAngle,2));
    maxspeed = single(max(Xr(:)));
    [MT, MTparams] = MTbank(Xt,Xr,maxspeed);
    MT = squeeze(MT);
    newMST = false;
    [MST, MSTparams] = MSTbank3(MT, MTparams,newMST);
%     [X1] = max(MST,[],2);
%     X1 = max(MST(:,:),[],2);
    X1max = max(MST(:,:),[],2);
    X1min = min(MST(:,:),[],2);
    X1(abs(X1max) > abs(X1min)) = X1max(abs(X1max) > abs(X1min));
    X1(abs(X1max) <= abs(X1min)) = X1min(abs(X1max) <= abs(X1min));
    NeuroResp(dircounter,loccounter,:) = X1(:);
%     close;
end
end
toc;
figure;
% baseline = 0;%min(NeuroResp(:));
% themax = max(NeuroResp(:));
% rg = [baseline-(themax-baseline),themax];
for n = 1:50
    thisResp = squeeze(NeuroResp(:,:,n));
    baseline = min(thisResp(:));
    themax = max(thisResp(:));
    rg = [baseline-(themax-baseline),themax];
    for i = 1:9
        subplot(3,3,i);colormap(jet);polarmosaic(squeeze(NeuroResp(:,i,n)),rg,.35,1);box off;
    end
    pause;close
end

%% To decode motion type

tic;
load ./StimulusParam.mat;
numTrials = 100;
DIR1 = pi;
DIR2 = 0;
allApertureLoc = [100 100; 300 100; 500 100; 100 300; 300 300; 500 300; 100 500; 300 500; 500 500];  
for trcounter = 1:numTrials/2
    fprintf(['trial ',num2str(trcounter),'\n']);
    M.apertureLoc = allApertureLoc(1,:);

    M.direction = DIR1;
    save('./StimulusParam.mat','M');
    
    [Stimulus] = generateStimulus();
    MotionField = estimateMotionField(Stimulus); 
    
    Xr = reshape(log(MotionField.averageMotionInsideAmp/10),1,size(MotionField.averageMotionInsideAmp,1)*size(MotionField.averageMotionInsideAmp,2));
    Xt = reshape(MotionField.averageMotionInsideAngle,1,size(MotionField.averageMotionInsideAngle,1)*size(MotionField.averageMotionInsideAngle,2));
    maxspeed = single(max(Xr(:)));
    [MT, MTparams] = MTbank(Xt,Xr,maxspeed);
    MT = squeeze(MT);
    newMST = false;
    [MST, MSTparams] = MSTbank3(MT, MTparams,newMST);
    allMST(trcounter,:,1) = reshape(MST,size(MST,1)*size(MST,2),1);
    allMT(trcounter,:,1) = reshape(MT,size(MT,1)*size(MT,2)*size(MT,3),1);
end
for trcounter = 1:numTrials/2
    fprintf(['trial ',num2str(trcounter),'\n']);
    M.apertureLoc = allApertureLoc(1,:);
    
    M.direction = DIR2;
    save('./StimulusParam.mat','M');
    
    [Stimulus] = generateStimulus();
    MotionField = estimateMotionField(Stimulus); 
    
    Xr = reshape(log(MotionField.averageMotionInsideAmp/10),1,size(MotionField.averageMotionInsideAmp,1)*size(MotionField.averageMotionInsideAmp,2));
    Xt = reshape(MotionField.averageMotionInsideAngle,1,size(MotionField.averageMotionInsideAngle,1)*size(MotionField.averageMotionInsideAngle,2));
    maxspeed = single(max(Xr(:)));
    [MT, MTparams] = MTbank(Xt,Xr,maxspeed);
    MT = squeeze(MT);
    newMST = false;
    [MST, MSTparams] = MSTbank3(MT, MTparams,newMST);
    allMST(trcounter,:,2) = reshape(MST,size(MST,1)*size(MST,2),1);
    allMT(trcounter,:,2) = reshape(MT,size(MT,1)*size(MT,2)*size(MT,3),1);
end

toc


%% pca on MT/MST population

allmt = reshape(allMT,100,21600);
allmst = reshape(allMST,100,2560);
allmst = (allmst - mean(allmst(:)))./(std(allmst(:)));
allmt = (allmt - mean(allmt(:)))./(std(allmt(:)));
[coeff1,score1,latent1] = pca([allmt,allmst]);
figure;plot(coeff1(:,1),'-');


%% decoding

load ./simulated' data'/data10
allMT_normal = (allMT - mean(allMT(:)))./(std(allMT(:)));
allMST_normal = (allMST - mean(allMST(:)))./(std(allMST(:)));
epcilon = 0.2;
allMT_noisy = allMT_normal + epcilon*randn(size(allMT_normal));
allMST_noisy = allMST_normal + epcilon*randn(size(allMST_normal));

% for subsamples = 1:30
%     subsamples
%     [dprimes,selectedNeuronsIdx,L] = adaptiveDecoder2(allMT,allMST,subsamples);
%     testL(subsamples) = L(end);
% end
% figure;plot(testL,'o-k');xlabel('subspace');ylabel('AIC');box off
% [~,sIdx] = min(testL);

% [dprimes,selectedNeuronsIdx,L,Rs] = adaptiveDecoder2(allMT,allMST,2);
[dprimes,selectedNeuronsIdx] = adaptiveDecoder(allMT_noisy,allMST_noisy,2);
selectedMST = sum(selectedNeuronsIdx(:,:) > size(allMT_noisy,2));
selectedMT = sum(selectedNeuronsIdx(:,:) <= size(allMT_noisy,2));
figure;subplot(2,1,1);plot(resample(selectedMST,1,10),'-b');hold on;plot(resample(selectedMT,1,10),'-r')
subplot(2,1,2);plot(resample(double(max(dprimes,[],1)),1,10),'o-')
% proportion of MT/MST
for iter = 1:50
    iter
    [dprimes,selectedNeuronsIdx] = adaptiveDecoder(allMT_noisy,allMST_noisy,2);
    selectedMST(iter) = sum(selectedNeuronsIdx(:,end) > size(allMT_noisy,2));
    selectedMT(iter) = sum(selectedNeuronsIdx(:,end) <= size(allMT_noisy,2));
end
figure;histogram(selectedMST);hold on;histogram(selectedMT);
cohenD = (mean(selectedMST) - mean(selectedMT))./std([selectedMST,selectedMT])
    