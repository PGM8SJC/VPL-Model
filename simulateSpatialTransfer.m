%% measure decoding and spatial transfer repeately

numtrial = 10;
% allApertureLoc = [100 100; 200 100; 300 100; 100 200; 200 200; 300 300; 100 300; 200 300; 300 300];
allApertureLoc = [100 100; 300 100; 500 100; 100 300; 300 300; 500 300; 100 500; 300 500; 500 500];
DIR1 = pi/2;
DIR2 = 3*pi/2;
tic;
for tr = 1:numtrial
    fprintf(['Trial #',num2str(tr),'\n']);
    % adaptive decoding
    load ./simulated' data'/data37
    allMT_normal = (allMT - mean(allMT(:)))./(std(allMT(:)));
    allMST_normal = (allMST - mean(allMST(:)))./(std(allMST(:)));
    epcilon = 0.2;
    allMT_noisy = allMT_normal + epcilon*randn(size(allMT_normal));
    allMST_noisy = allMST_normal + epcilon*randn(size(allMST_normal));
    numReadoutNeurons = 2;
    
    [dprimes,selectedNeuronsIdx] = adaptiveDecoder(allMT_noisy,allMST_noisy,numReadoutNeurons);
    [dmax,id] = max(dprimes(:,end));
    readoutNeuronIdx = selectedNeuronsIdx(id,end);
    selectedMST = sum(selectedNeuronsIdx(:,:) > size(allMT_noisy,2));
    selectedMT = sum(selectedNeuronsIdx(:,:) <= size(allMT_noisy,2));
    
    % measure spatial transfer
    load ./StimulusParam.mat;
    dircounter = 0;
    
    
    for loccounter = 1:9
        
        for trcounter = 1:20

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
        fprintf('.\n')
    end
    
    epcilon = 0.2;
    readoutNeuronResp_normal = (readoutNeuronResp - mean(readoutNeuronResp(:)))./std(readoutNeuronResp(:));
    readoutNeuronResp_noisy = readoutNeuronResp_normal + epcilon*randn(size(readoutNeuronResp_normal));
    mean1 = mean(readoutNeuronResp_noisy(:,1,:),3);
    mean2 = mean(readoutNeuronResp_noisy(:,2,:),3);
    var1 = var(readoutNeuronResp_noisy(:,1,:),[],3);
    var2 = var(readoutNeuronResp_noisy(:,2,:),[],3);
    
    dprime_everywhere(:,tr) = abs(mean1 - mean2)./sqrt(0.5 * (var1 + var2));

    
    
end
toc;
