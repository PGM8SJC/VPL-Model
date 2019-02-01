[Stimulus] = generateStimulus();
MotionField = estimateMotionField(Stimulus);

Xr = reshape(log(MotionField.averageMotionInsideAmp/10),1,size(MotionField.averageMotionInsideAmp,1)*size(MotionField.averageMotionInsideAmp,2));
Xt = reshape(MotionField.averageMotionInsideAngle,1,size(MotionField.averageMotionInsideAngle,1)*size(MotionField.averageMotionInsideAngle,2));
maxspeed = single(max(Xr(:)));
[X, params] = MTbank(Xt,Xr,maxspeed);
X = squeeze(X);

% X1 = max(X1,[],3);
% X1 = max(X1,[],1);
% figure;plot(0:(pi/8):(2*pi - pi/8),X1,'-o'); % show orientation responses
% 
% X1 = max(X1,[],1);
% X1 = max(X1,[],2);
% X1 = squeeze(X1);
% X1 = reshape(X1,25,25);
% figure;imagesc(X1)


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