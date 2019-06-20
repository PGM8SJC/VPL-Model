function [dprimes,selectedNeuronsIdx] = adaptiveDecoder(MTresp,MSTresp)

motion1 = 1;
motion2 = 2;

numNeuronsMT = size(MTresp,2);
numNeuronsMST = size(MSTresp,2);
numTrials = size(MSTresp,1);

neurons = cat(2,MTresp,MSTresp);%MTresp;
numNeuronsAll = size(neurons,2);
subsamples = 10;
whichNeuronsIdx = randi(numNeuronsAll,subsamples,1);

for trcount = 1:100
    if trcount == 1
        selectedNeuronsIdx(:,trcount) = whichNeuronsIdx;
    else
        selectedNeuronsIdx(:,trcount) = selectedNeuronsIdx(:,trcount-1);
        selectedNeuronsIdx(dprimeminIdx,trcount) = whichNeuronsIdx;
    end
    selectedNeurons = (neurons(:,whichNeuronsIdx,:));%squeeze
    if trcount == 1
        meanPref = mean(selectedNeurons(:,:,motion1),1);
        meanNull = mean(selectedNeurons(:,:,motion2),1);
        varPref = var(selectedNeurons(:,:,motion1),[],1);
        varNull = var(selectedNeurons(:,:,motion2),[],1);
    else
        meanPref = mean(selectedNeurons(:,:,motion1),1);
        meanNull = mean(selectedNeurons(:,:,motion2),1);
        varPref = var(selectedNeurons(:,:,motion1),[],1);
        varNull = var(selectedNeurons(:,:,motion2),[],1);    
    end

    dprimetemp = abs(meanPref - meanNull)./sqrt(.5 * (varPref + varNull));
    if trcount == 1
        dprimes(:,trcount) = dprimetemp;
    else
        dprimes(:,trcount) = dprimes(:,trcount - 1);
        dprimes(dprimeminIdx,trcount) = dprimetemp;
    end
        
    [dptimemin,dprimeminIdx] = min(dprimes(:,trcount));
    whichNeuronsIdx = randi(numNeuronsAll,1);
    
    
    
end


end