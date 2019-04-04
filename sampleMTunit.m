%% simulate MT population responses from one stimulus
[Stimulus] = generateStimulus();
MotionField = estimateMotionField(Stimulus);

Xr = reshape(log(MotionField.averageMotionInsideAmp/10),1,size(MotionField.averageMotionInsideAmp,1)*size(MotionField.averageMotionInsideAmp,2));
Xt = reshape(MotionField.averageMotionInsideAngle,1,size(MotionField.averageMotionInsideAngle,1)*size(MotionField.averageMotionInsideAngle,2));
maxspeed = single(max(Xr(:)));
[X, params] = MTbank(Xt,Xr,maxspeed);
X = squeeze(X);


%% simulate responses of a MT population to a range of stimuli 
load ./StimulusParam.mat;
i = 0;
for dir = 0:pi/6:(2*pi - pi/6)
    i = i + 1;
    M.direction = dir;
    save('./StimulusParam.mat','M');
    
    [Stimulus] = generateStimulus();
    MotionField = estimateMotionField(Stimulus);
    
    Xr = reshape(log(MotionField.averageMotionInsideAmp/10),1,size(MotionField.averageMotionInsideAmp,1)*size(MotionField.averageMotionInsideAmp,2));
    Xt = reshape(MotionField.averageMotionInsideAngle,1,size(MotionField.averageMotionInsideAngle,1)*size(MotionField.averageMotionInsideAngle,2));
    maxspeed = single(max(Xr(:)));
    [X, params] = MTbank(Xt,Xr,maxspeed);
    
    X1 = squeeze(X);
    X1 = max(X1,[],3);
    NeuroResp(i,:) = X1(3,:);
    close;
end

%% 
load ./StimulusParam.mat;
cohlevels = 0.01:0.02:.2;
for cohcount = 1:length(cohlevels)
    M.coherence = cohlevels(cohcount);
    save('./StimulusParam.mat','M');
    for trialcount = 1:20
        [Stimulus] = generateStimulus();
        MotionField = estimateMotionField(Stimulus);
        
        Xr = reshape(log(MotionField.averageMotionInsideAmp/10),1,size(MotionField.averageMotionInsideAmp,1)*size(MotionField.averageMotionInsideAmp,2));
        Xt = reshape(MotionField.averageMotionInsideAngle,1,size(MotionField.averageMotionInsideAngle,1)*size(MotionField.averageMotionInsideAngle,2));
        maxspeed = single(max(Xr(:)));
        [X, params] = MTbank(Xt,Xr,maxspeed);
        
        X1 = squeeze(X);
        X1 = max(X1,[],3);
        X1 = max(X1,[],1);
        
        NeuroResp(cohcount,trialcount) = X1(5);
        close all;
    end
end