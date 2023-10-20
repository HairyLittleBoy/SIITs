function [Pks, Lhloading, Lhunloading, creep] = Lhsplit ( obj )

if obj.unloadingNum == 0    % if the Lh curve has no loading-unloading cycles,
    % only load  and unloading, then
    % obj.unloadingNum = 0
    
    
    %-----------------------------------------------------------------------
    % This function divide the L-h curve into loading and
    % unloading parts. The L-h curve only has single loading
    % and unloading cycle
    % The Lhloading and Lhunloading are two matrixes
    %-----------------------------------------------------------------------
    
    LhCurveMatrixU = unique(obj.LhPoints,'rows','stable');
    LhCurveMatrixU(:,1) = smooth(LhCurveMatrixU(:,1),20);
    LhCurveMatrixU(:,2) = smooth(LhCurveMatrixU(:,2),20);
    
    [~,maxLoadLoc] = max(LhCurveMatrixU(:,2));
    
    Lhloading{1}(:,1) = LhCurveMatrixU(1:maxLoadLoc,1);
    Lhloading{1}(:,2) = sqrt(2*obj.indenterR.*Lhloading.part1(:,1) - Lhloading.part1(:,1).^2);
    Lhloading{1}(:,3) = LhCurveMatrixU(1:maxLoadLoc,2);
    Lhloading{1}(:,4) = Lhloading.part1(:,3)./(pi.*Lhloading.part1(:,2).^2);
    
    Lhunloading{1}(:,1) = LhCurveMatrixU(maxLoadLoc:end,1);
    Lhunloading{1}(:,2) = sqrt(2*obj.indenterR.*Lhunloading.part1(:,1) - Lhunloading.part1(:,1).^2);
    Lhunloading{1}(:,3) = LhCurveMatrixU(maxLoadLoc:end,2);
    Lhunloading{1}(:,4) = Lhunloading.part1(:,3)./(pi.*Lhunloading.part1(:,2).^2);
    
    
    Pks.Low = [Lhunloading{1}(end,1),Lhunloading{1}(end,2)];
    Pks.High1 = [LhCurveMatrixU(maxLoadLoc,1), LhCurveMatrixU(maxLoadLoc,2)];
    Pks.High2 = [LhCurveMatrixU(maxLoadLoc,1), LhCurveMatrixU(maxLoadLoc,2)];
    
elseif obj.unloadingNum > 0
    
    %----------------------------------------------------------------------------
    % This function divide the L-h curve into loading and
    % unloading parts. The L-h curve  has multiple loading
    % and unloading cycles.
    % The Lhloading and Lhunloading are two structs containing
    % several matrixes in which the penetraion is the 1st colume,
    % contact radius is the 2nd colume, load the 3rd colume and
    % contact pressure the 4th colume.
    %
    % < This LhDivideMtpl function might be useful to Lh curve with
    %    single loading-unloading cycle, this will be tested later >
    %----------------------------------------------------------------------------
    
    LhCurveMatrixU = unique(obj.LhPoints,'rows','stable') ;
    LhCurveMatrixU(:,2) = smooth(LhCurveMatrixU(:,2),20);
    LhCurveMatrixU(:,1) = smooth(LhCurveMatrixU(:,1),20);
    
    
    [maxLoad,maxLoadNum] = max(LhCurveMatrixU(:,2)) ;
    maxLoadLoc = LhCurveMatrixU(maxLoadNum,1);
    peneNum=[1:length(LhCurveMatrixU)];
    
    % pkNumLow is the number of the unloading bottom points in the whole LhCurveMatrixU data
    [pksLow, pkNumLow] = findpeaks(-LhCurveMatrixU(:,2), peneNum,  ...
        'MinPeakHeight',-10, 'SortStr','ascend', 'NPeaks',obj.unloadingNum) ;
    
        [pksLow, pkNumLow] = findpeaks(-LhCurveMatrixU(:,2), peneNum,  ...
        'MinPeakHeight',-10, 'SortStr','ascend', 'NPeaks',obj.unloadingNum) ;
    
    counter = 1;
    while length(pksLow) < obj.unloadingNum
        counter = counter+1;
%   pksNum is smaller than actual lowest point on Lh curve,change the -maxLoad
        [pksLow, pkNumLow] = findpeaks(-LhCurveMatrixU(:,2), peneNum,  ...
            'MinPeakHeight',-10*counter, 'SortStr','ascend', 'NPeaks',obj.unloadingNum) ;
    end
    
    pksLow = -pksLow;
    
    [pkNumLow, index1]= sort(pkNumLow);   % have to be sorted 
    pksLow = pksLow(index1);
    pksLocLow = LhCurveMatrixU(pkNumLow, 1);
    
    pkNumHigh1 = zeros(obj.unloadingNum,1);
    pksHigh1 = zeros(obj.unloadingNum,1);
    pksLocHigh1 = zeros(obj.unloadingNum,1);
    
    pksHigh2 = zeros(length(pkNumHigh1), 1);
    pkNumHigh2 = zeros(length(pkNumHigh1) , 1);
    
    switch obj.creepOrNot
    case 'creep'

        [pksTemp1, pkNumTemp1] = findpeaks(LhCurveMatrixU(:,2), peneNum, 'MinPeakHeight', 10) ;
        pksLocTemp1 = LhCurveMatrixU(pkNumTemp1, 1);
        
        % the maximum creep displacement is  set to 5 micrometer
        pksLocTemp1Change = find((pksLocTemp1(2:end) - pksLocTemp1(1:length(pksLocTemp1)-1)) > 5e-3);
        
        creepNum = cell(1,obj.unloadingNum);
        creepPart = cell(1,obj.unloadingNum);
        for i = 1 : obj.unloadingNum+1
            if i == 1
                creepNum{i} = find(abs(LhCurveMatrixU(:,2) - mean(pksTemp1(1:pksLocTemp1Change(i)))) < 5);
            elseif i == obj.unloadingNum+1
                creepNum{i} = find(abs(LhCurveMatrixU(:,2) - mean(pksTemp1(pksLocTemp1Change(i-1)+1:end))) < 5);
            else
                creepNum{i} = find(abs(LhCurveMatrixU(:,2) - mean(pksTemp1(pksLocTemp1Change(i-1)+1:pksLocTemp1Change(i)))) < 5);
            end
            for j = 1:length(creepNum{i})-1
                temp1 = creepNum{i}(j+1) - creepNum{i}(j);
                if temp1 ~= 1
                    creepNum{i}(j:end) = [];
                    break
                end
            end
            creepPart{i} = LhCurveMatrixU(creepNum{i},:);
            pkNumHigh1(i) = creepNum{i}(end);
            pksHigh1(i) = creepPart{i}(end,2);
            pksLocHigh1(i) = creepPart{i}(end,1);
        end
        creep.creepPart = creepPart;
        creep.creepNum = creepNum;
    case 'nocreep'
        % pkNumHigh1 is the number of the starting point of unloading part in the whole LhCurveMatrixU data
        [pksHigh1, pkNumHigh1] = findpeaks(LhCurveMatrixU(:,2), peneNum, 'MinPeakHeight', 10) ;
        pksLocHigh1 = LhCurveMatrixU(pkNumHigh1, 1);
        creep = [];
    end
    % pkNumHigh2 is the number of the end point of unloading part in the whole LhCurveMatrixU data
    
    for i = 1 : length(pkNumHigh1) - 1
        [~,tempNum] = min(abs(pksHigh1(i) - LhCurveMatrixU(pkNumLow(i) : pkNumHigh1(i+1), 2)));
        pksHigh2(i) = LhCurveMatrixU(pkNumLow(i) + (tempNum-1), 2);
        pkNumHigh2(i) = pkNumLow(i) + (tempNum-1);
    end
    pksLocHigh2 = LhCurveMatrixU(pkNumHigh2, 1);

% --------------------------------------------------    
    Lhunloading = cell(1,length(pksLocHigh2)+1);
    % Lhunloading is a cell, Lhunloading{1} is the unloading part of 1st
    % loading-unloading cycle.
    % Lhunloading{i}(:,1)  : penetration
    % Lhunloading{i}(:,2)  : contact radius a
    % Lhunloading{i}(:,3)  : load
    % Lhunloading{i}(:,4)  : contact pressure pm
    % Lhunloading{i}(:,5)  : Tangent Slope
    
    Lhloading = cell(1,length(pkNumHigh1));
    % Lhloading is also a cell, Lhloading{1} is the loading part before 1st
    % loading-unloading cycle.
    % Lhloading{i}(:,1)  : penetration
    % Lhloading{i}(:,2)  : contact radius a
    % Lhloading{i}(:,3)  : load
    % Lhloading{i}(:,4)  : contact pressure pm
    % Lhloading{i}(:,5)  : Tangent Slope
    switch obj.creepOrNot
        case 'creep'
            for i = 1 : length(pkNumHigh1)
                if i == 1
                    Lhloading{i}(:,1) = LhCurveMatrixU(1:creepNum{i}(1),1);
                    Lhloading{i}(:,2) = sqrt(2*obj.indenterR.*Lhloading{i}(:,1) - Lhloading{i}(:,1).^2);
                    Lhloading{i}(:,3) = LhCurveMatrixU(1:creepNum{i}(1),2);
                    Lhloading{i}(:,4) = Lhloading{i}(:,3)./(pi.*Lhloading{i}(:,2).^2);
                else
                    Lhloading{i}(:,1) = LhCurveMatrixU(pkNumHigh2(i-1):creepNum{i}(1),1);
                    Lhloading{i}(:,2) = sqrt(2*obj.indenterR.*Lhloading{i}(:,1) - Lhloading{i}(:,1).^2);
                    Lhloading{i}(:,3) = LhCurveMatrixU(pkNumHigh2(i-1):creepNum{i}(1),2);
                    Lhloading{i}(:,4) = Lhloading{i}(:,3)./(pi.*Lhloading{i}(:,2).^2);
                end
                Lhloading{i}(2:end,5) = smooth(diff(Lhloading{i}(:,3)) ./ diff(Lhloading{i}(:,1)),200);
                Lhloading{i}(1,5) = Lhloading{i}(2,5);
            end  
            for i = 1 : length(pksLocHigh2)+1
                if i < length(pksLocHigh2)+1
                    Lhunloading{i}(:,1) = LhCurveMatrixU(pkNumHigh1(i)+1 : pkNumLow(i)-1,1);
                    Lhunloading{i}(:,2) = sqrt(2*obj.indenterR.*Lhunloading{i}(:,1) - Lhunloading{i}(:,1).^2);
                    Lhunloading{i}(:,3) = LhCurveMatrixU(pkNumHigh1(i)+1 : pkNumLow(i)-1,2);
                    Lhunloading{i}(:,4) = Lhunloading{i}(:,3)./(pi.*Lhunloading{i}(:,2).^2);
                else
                    Lhunloading{i}(:,1) = LhCurveMatrixU(creepNum{obj.unloadingNum+1}(end):end,1);
                    Lhunloading{i}(:,2) = sqrt(2*obj.indenterR.*Lhunloading{i}(:,1) - Lhunloading{i}(:,1).^2);
                    Lhunloading{i}(:,3) = LhCurveMatrixU(creepNum{obj.unloadingNum+1}(end):end,2);
                    Lhunloading{i}(:,4) = Lhunloading{i}(:,3)./(pi.*Lhunloading{i}(:,2).^2);
                end
                Lhunloading{i}(2:end,5) = smooth(diff(Lhunloading{i}(:,3)) ./ diff(Lhunloading{i}(:,1)),200);
                Lhunloading{i}(1,5) = Lhunloading{i}(2,5);
            end

        case 'nocreep'
            for i = 1 : length(pkNumHigh1)
                if i == 1
                    Lhloading{i}(:,1) = LhCurveMatrixU(1:pkNumHigh1(i),1);
                    Lhloading{i}(:,2) = sqrt(2*obj.indenterR.*Lhloading{i}(:,1) - Lhloading{i}(:,1).^2);
                    Lhloading{i}(:,3) = LhCurveMatrixU(1:pkNumHigh1(i),2);
                    Lhloading{i}(:,4) = Lhloading{i}(:,3)./(pi.*Lhloading{i}(:,2).^2);
                else
                    Lhloading{i}(:,1) = LhCurveMatrixU(pkNumHigh2(i-1):pkNumHigh1(i),1);
                    Lhloading{i}(:,2) = sqrt(2*obj.indenterR.*Lhloading{i}(:,1) - Lhloading{i}(:,1).^2);
                    Lhloading{i}(:,3) = LhCurveMatrixU(pkNumHigh2(i-1):pkNumHigh1(i),2);
                    Lhloading{i}(:,4) = Lhloading{i}(:,3)./(pi.*Lhloading{i}(:,2).^2);
                end
                Lhloading{i}(2:end,5) = smooth(diff(Lhloading{i}(:,3)) ./ diff(Lhloading{i}(:,1)),200);
                Lhloading{i}(1,5) = Lhloading{i}(2,5);                
            end
            for i = 1 : length(pksLocHigh2)+1
                if i < length(pksLocHigh2)+1
                    Lhunloading{i}(:,1) = LhCurveMatrixU(pkNumHigh1(i)+1 : pkNumLow(i)-1,1);
                    Lhunloading{i}(:,2) = sqrt(2*obj.indenterR.*Lhunloading{i}(:,1) - Lhunloading{i}(:,1).^2);
                    Lhunloading{i}(:,3) = LhCurveMatrixU(pkNumHigh1(i)+1 : pkNumLow(i)-1,2);
                    Lhunloading{i}(:,4) = Lhunloading{i}(:,3)./(pi.*Lhunloading{i}(:,2).^2);
                else
                    Lhunloading{i}(:,1) = LhCurveMatrixU(maxLoadNum:end,1);
                    Lhunloading{i}(:,2) = sqrt(2*obj.indenterR.*Lhunloading{i}(:,1) - Lhunloading{i}(:,1).^2);
                    Lhunloading{i}(:,3) = LhCurveMatrixU(maxLoadNum:end,2);
                    Lhunloading{i}(:,4) = Lhunloading{i}(:,3)./(pi.*Lhunloading{i}(:,2).^2);
                end
                Lhunloading{i}(2:end,5) = smooth(diff(Lhunloading{i}(:,3)) ./ diff(Lhunloading{i}(:,1)),200);
                Lhunloading{i}(1,5) = Lhunloading{i}(2,5);          
            end
    end
    Pks.Low = [pksLocLow, pksLow];
    Pks.High1 = [pksLocHigh1, pksHigh1];
    Pks.High2 = [pksLocHigh2, pksHigh2];
    
else
    error('  the num of unloading should be specified  ')
end

end
