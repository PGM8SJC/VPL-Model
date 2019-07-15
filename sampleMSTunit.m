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
% figure;
% % baseline = 0;%min(NeuroResp(:));
% % themax = max(NeuroResp(:));
% % rg = [baseline-(themax-baseline),themax];
% for n = 1:50
%     thisResp = squeeze(NeuroResp(:,:,n));
%     baseline = min(thisResp(:));
%     themax = max(thisResp(:));
%     rg = [baseline-(themax-baseline),themax];
%     for i = 1:9
%         subplot(3,3,i);colormap(jet);polarmosaic(squeeze(NeuroResp(:,i,n)),rg,.35,1);box off;
%     end
%     pause;close
% end

%% simulate MT/MST to feed to adaptive decoder

tic;
load ./StimulusParam.mat;
numTrials = 100;
DIR1 = pi/2;
DIR2 = 3*pi/2;
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



%% adaptive decoding

% load ./simulated' data'/data32
allMT_normal = (allMT - mean(allMT(:)))./(std(allMT(:)));
allMST_normal = (allMST - mean(allMST(:)))./(std(allMST(:)));
epcilon = 0.2;
allMT_noisy = allMT_normal + epcilon*randn(size(allMT_normal));
allMST_noisy = allMST_normal + epcilon*randn(size(allMST_normal));
numReadoutNeurons = 2;

% [dprimes,selectedNeuronsIdx,L,Rs] = adaptiveDecoder2(allMT,allMST,numReadoutNeurons);
[dprimes,selectedNeuronsIdx] = adaptiveDecoder(allMT_noisy,allMST_noisy,numReadoutNeurons);
[dmax,id] = max(dprimes(:,end));
readoutNeuronIdx = selectedNeuronsIdx(id,end);
selectedMST = sum(selectedNeuronsIdx(:,:) > size(allMT_noisy,2));
selectedMT = sum(selectedNeuronsIdx(:,:) <= size(allMT_noisy,2));
figure;subplot(2,1,1);plot(resample(selectedMST,1,10),'-b');hold on;plot(resample(selectedMT,1,10),'-r')
subplot(2,1,2);plot(resample(double(max(dprimes,[],1)),1,10),'o-')
% % proportion of MT/MST
% for iter = 1:50
%     iter
%     [dprimes,selectedNeuronsIdx] = adaptiveDecoder(allMT_noisy,allMST_noisy,numReadoutNeurons);
%     selectedMST(iter) = sum(selectedNeuronsIdx(:,end) > size(allMT_noisy,2));
%     selectedMT(iter) = sum(selectedNeuronsIdx(:,end) <= size(allMT_noisy,2));
% end
% figure;histogram(selectedMST);hold on;histogram(selectedMT);
% cohenD = (mean(selectedMST) - mean(selectedMT))./std([selectedMST,selectedMT])
%     

%% tuning of the readout neuron
load ./StimulusParam.mat;
dircounter = 0;
% allApertureLoc = [100 100; 200 100; 300 100; 100 200; 200 200; 300 300; 100 300; 200 300; 300 300];
allApertureLoc = [100 100; 300 100; 500 100; 100 300; 300 300; 500 300; 100 500; 300 500; 500 500];  
for dir = 0:pi/4:(2*pi - pi/4)
    dircounter = dircounter + 1;
    M.direction = dir;
for loccounter = 1:9
    

    
    fprintf(['location ',num2str(loccounter), ' dir ',num2str(dircounter),' ']);
    M.apertureLoc = allApertureLoc(loccounter,:);
    
    % get simple motion responses
    fprintf(['simple',' ']);
    M.motiontype = 'simple';
    save('./StimulusParam.mat','M');
    
    [Stimulus] = generateStimulus();
    MotionField = estimateMotionField(Stimulus); 
    
    Xr = reshape(log(MotionField.averageMotionInsideAmp/10),1,size(MotionField.averageMotionInsideAmp,1)*size(MotionField.averageMotionInsideAmp,2));
    Xt = reshape(MotionField.averageMotionInsideAngle,1,size(MotionField.averageMotionInsideAngle,1)*size(MotionField.averageMotionInsideAngle,2));
    maxspeed = single(max(Xr(:)));
    [MT, MTparams] = MTbank(Xt,Xr,maxspeed);
    MT = squeeze(MT);
    numMT = size(MT,1)*size(MT,2)*size(MT,3);
    newMST = false;
    [MST, MSTparams] = MSTbank3(MT, MTparams,newMST);
    numMST = size(MST,1)*size(MST,2);
    allMST = reshape(MST,size(MST,1)*size(MST,2),1);
    allMT= reshape(MT,size(MT,1)*size(MT,2)*size(MT,3),1);
    if readoutNeuronIdx > numMT
        readoutNeuronResp_simple(dircounter,loccounter) = allMST(readoutNeuronIdx-numMT);
    else
        readoutNeuronResp_simple(dircounter,loccounter) = allMT(readoutNeuronIdx);
    end
    
    % get complex motion responses
    fprintf(['complex','\n']);
    M.motiontype = 'complex';
    save('./StimulusParam.mat','M');
    [Stimulus] = generateStimulus();
    MotionField = estimateMotionField(Stimulus); 
    
    Xr = reshape(log(MotionField.averageMotionInsideAmp/10),1,size(MotionField.averageMotionInsideAmp,1)*size(MotionField.averageMotionInsideAmp,2));
    Xt = reshape(MotionField.averageMotionInsideAngle,1,size(MotionField.averageMotionInsideAngle,1)*size(MotionField.averageMotionInsideAngle,2));
    maxspeed = single(max(Xr(:)));
    [MT, MTparams] = MTbank(Xt,Xr,maxspeed);
    MT = squeeze(MT);
    numMT = size(MT,1)*size(MT,2)*size(MT,3);
    newMST = false;
    [MST, MSTparams] = MSTbank3(MT, MTparams,newMST);
    numMST = size(MST,1)*size(MST,2);
    allMST = reshape(MST,size(MST,1)*size(MST,2),1);
    allMT= reshape(MT,size(MT,1)*size(MT,2)*size(MT,3),1);
    if readoutNeuronIdx > numMT
        readoutNeuronResp_complex(dircounter,loccounter) = allMST(readoutNeuronIdx-numMT);
    else
        readoutNeuronResp_complex(dircounter,loccounter) = allMT(readoutNeuronIdx);
    end


%     close;
end
end

figure;title('response to complex rdk')
thisResp = readoutNeuronResp_complex(:,:);
baseline = 0;%min(thisResp(:));
themax = max(thisResp(:));
rg = [baseline-(themax-baseline),themax];
for i = 1:9
    subplot(3,3,i);colormap(jet);polarmosaic(squeeze(thisResp(:,i)),rg,.35,1);box off;
end

figure;title('response to simple rdk')
thisResp = readoutNeuronResp_simple(:,:);
% baseline = 0;%min(thisResp(:));
% themax = max(thisResp(:));
rg = [baseline-(themax-baseline),themax];
for i = 1:9
    subplot(3,3,i);colormap(jet);polarmosaic(squeeze(thisResp(:,i)),rg,.35,1);box off;
end

%% dprime of the readout neurons at different locations

load ./StimulusParam.mat;
dircounter = 0;
% allApertureLoc = [100 100; 200 100; 300 100; 100 200; 200 200; 300 300; 100 300; 200 300; 300 300];
allApertureLoc = [100 100; 300 100; 500 100; 100 300; 300 300; 500 300; 100 500; 300 500; 500 500];
DIR1 = pi/2;
DIR2 = 3*pi/2;
tic;
for loccounter = 1:9
    fprintf(['location ',num2str(loccounter),' ']);
    for trcounter = 1:20
        
%         [dprimes,selectedNeuronsIdx] = adaptiveDecoder(allMT_noisy,allMST_noisy,numReadoutNeurons);
%         [dmax,id] = max(dprimes(:,end));
%         readoutNeuronIdx = selectedNeuronsIdx(id,end);
%         
        fprintf(['.']); 
        
        
        M.direction = DIR1;
        M.apertureLoc = allApertureLoc(loccounter,:);
        M.motiontype = 'simple';
        save('./StimulusParam.mat','M');
        
        [Stimulus] = generateStimulus();
        MotionField = estimateMotionField(Stimulus);
        
        Xr = reshape(log(MotionField.averageMotionInsideAmp/10),1,size(MotionField.averageMotionInsideAmp,1)*size(MotionField.averageMotionInsideAmp,2));
        Xt = reshape(MotionField.averageMotionInsideAngle,1,size(MotionField.averageMotionInsideAngle,1)*size(MotionField.averageMotionInsideAngle,2));
        maxspeed = single(max(Xr(:)));
        [MT, MTparams] = MTbank(Xt,Xr,maxspeed);
        MT = squeeze(MT);
        numMT = size(MT,1)*size(MT,2)*size(MT,3);
        newMST = false;
        [MST, MSTparams] = MSTbank3(MT, MTparams,newMST);
        numMST = size(MST,1)*size(MST,2);
        allMST = reshape(MST,size(MST,1)*size(MST,2),1);
        allMT= reshape(MT,size(MT,1)*size(MT,2)*size(MT,3),1);
        if readoutNeuronIdx > numMT
            readoutNeuronResp(loccounter,1,trcounter) = allMST(readoutNeuronIdx-numMT);
        else
            readoutNeuronResp(loccounter,1,trcounter) = allMT(readoutNeuronIdx);
        end
        
        M.direction = DIR2;
        M.apertureLoc = allApertureLoc(loccounter,:);
        M.motiontype = 'simple';
        save('./StimulusParam.mat','M');
        
        [Stimulus] = generateStimulus();
        MotionField = estimateMotionField(Stimulus);
        
        Xr = reshape(log(MotionField.averageMotionInsideAmp/10),1,size(MotionField.averageMotionInsideAmp,1)*size(MotionField.averageMotionInsideAmp,2));
        Xt = reshape(MotionField.averageMotionInsideAngle,1,size(MotionField.averageMotionInsideAngle,1)*size(MotionField.averageMotionInsideAngle,2));
        maxspeed = single(max(Xr(:)));
        [MT, MTparams] = MTbank(Xt,Xr,maxspeed);
        MT = squeeze(MT);
        numMT = size(MT,1)*size(MT,2)*size(MT,3);
        newMST = false;
        [MST, MSTparams] = MSTbank3(MT, MTparams,newMST);
        numMST = size(MST,1)*size(MST,2);
        allMST = reshape(MST,size(MST,1)*size(MST,2),1);
        allMT= reshape(MT,size(MT,1)*size(MT,2)*size(MT,3),1);
        if readoutNeuronIdx > numMT
            readoutNeuronResp(loccounter,2,trcounter) = allMST(readoutNeuronIdx-numMT);
        else
            readoutNeuronResp(loccounter,2,trcounter) = allMT(readoutNeuronIdx);
        end
        
        
    end
    fprintf('\n')
end
toc;
epcilon = 0.15;
readoutNeuronResp_normal = (readoutNeuronResp - mean(readoutNeuronResp(:)))./std(readoutNeuronResp(:));
readoutNeuronResp_noisy = readoutNeuronResp_normal + epcilon*randn(size(readoutNeuronResp_normal));
mean1 = mean(readoutNeuronResp_noisy(:,1,:),3);
mean2 = mean(readoutNeuronResp_noisy(:,2,:),3);
var1 = var(readoutNeuronResp_noisy(:,1,:),[],3);
var2 = var(readoutNeuronResp_noisy(:,2,:),[],3);

dprime_everywhere = abs(mean1 - mean2)./sqrt(0.5 * (var1 + var2));
figure;imagesc(reshape(log(dprime_everywhere),3,3));

%% difference of the two training conditions
diff = reshape(log(dprime_everywhere_postcomplex) - log(dprime_everywhere_postsimple),3,3);
figure;plot(0,diff(1,1),'o');hold on
plot(1,(diff(2,1)+diff(1,2))./2,'ok');
plot(2,(diff(1,3)+diff(3,1))./2,'ok');
plot(3,(diff(3,2)+diff(2,3))./2,'ok');
plot(4,diff(3,3),'ok');

