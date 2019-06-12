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
