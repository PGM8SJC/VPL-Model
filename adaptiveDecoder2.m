function [dprimes,selectedNeuronsIdx,L,Rs] = adaptiveDecoder2(MTresp,MSTresp,subsamples)
warning off
motion1 = 1;
motion2 = 2;

numNeuronsMT = size(MTresp,2);
numNeuronsMST = size(MSTresp,2);
numTrials = size(MSTresp,1);

neurons = cat(2,MTresp,MSTresp);%MTresp;
MSTresp = (MSTresp - nanmean(MSTresp(:)))./nanstd(MSTresp(:));
MTresp = (MTresp - nanmean(MTresp(:)))./nanstd(MTresp(:));

numNeuronsAll = size(neurons,2);
numExtraMT = numNeuronsMT - numNeuronsMST;
numRepMST = floor(numExtraMT./numNeuronsMST);
neuronsIdx = [1:numNeuronsAll,repmat((numNeuronsMT+1):numNeuronsAll,1,numRepMST)];
numIdx = length(neuronsIdx);
whichIdx = randperm(numIdx,subsamples);
whichNeuronsIdx = neuronsIdx(whichIdx);

for trcount = 1:5000
    if trcount == 1
        selectedNeuronsIdx(:,trcount) = whichNeuronsIdx;
    else
        selectedNeuronsIdx(:,trcount) = selectedNeuronsIdx(:,trcount-1);
        selectedNeuronsIdx(dprimeminIdx,trcount) = whichNeuronsIdx;
    end
    selectedNeurons = (neurons(:,selectedNeuronsIdx(:,trcount),:));%squeeze
%     if trcount == 1
%         meanPref = mean(selectedNeurons(:,:,motion1),1);
%         meanNull = mean(selectedNeurons(:,:,motion2),1);
%         varPref = var(selectedNeurons(:,:,motion1),[],1);
%         varNull = var(selectedNeurons(:,:,motion2),[],1);
%     else
%         meanPref = mean(selectedNeurons(:,:,motion1),1);
%         meanNull = mean(selectedNeurons(:,:,motion2),1);
%         varPref = var(selectedNeurons(:,:,motion1),[],1);
%         varNull = var(selectedNeurons(:,:,motion2),[],1);    
%     end
    X = [selectedNeurons(:,:,motion1);selectedNeurons(:,:,motion2)];
    y = [zeros(numTrials,1);ones(numTrials,1)];
    
    mdl = fitglm(X,y,'Distribution','binomial','Link','logit');
    W = mdl.Coefficients.Estimate;
    W(1) = [];
    dprimetemp = abs(W)./sum(abs(W));
    dprimes(:,trcount) = dprimetemp;
%     dprimetemp = abs(meanPref - meanNull)./sqrt(.5 * (varPref + varNull));
%     if trcount == 1
%         dprimes(:,trcount) = dprimetemp;
%     else
%         dprimes(:,trcount) = dprimes(:,trcount - 1);
%         dprimes(dprimeminIdx,trcount) = dprimetemp;
%     end
        
    [dptimemin,dprimeminIdx] = min(abs(dprimes(:,trcount)));
    whichIdx = randperm(numIdx,1);
    whichNeuronsIdx = neuronsIdx(whichIdx);
    L(trcount) = mdl.ModelCriterion.AIC;
    Rs(trcount) = mdl.Rsquared.Ordinary;
    
end


end

