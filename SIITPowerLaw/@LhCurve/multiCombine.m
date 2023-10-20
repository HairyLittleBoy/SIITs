function LhloadingWhole = multiCombine ( obj )
% this function combine the loading part of Lh curve with mutiple
% loading-unloading cycles, for some identification method only apply to
% Lh curve with single loading-unloading
% the input Lhloading is one of the output of Lhsplit
[ ~, Lhloading, ~ ] = Lhsplit( obj );
LhloadingWholeLength = 0;
startPoint = zeros(length(Lhloading),1);
for i = 1 : length(Lhloading)
    LhloadingWholeLength = length(Lhloading{i}) + LhloadingWholeLength;
    startPoint(i) = LhloadingWholeLength + 1;
end
LhloadingWhole = zeros(LhloadingWholeLength,2);
for i = 1 : length(Lhloading)
    if i == 1
        LhloadingWhole(1: length(Lhloading{i}),1) = Lhloading{i}(:,1);
        LhloadingWhole(1: length(Lhloading{i}),2) = Lhloading{i}(:,2);
        LhloadingWhole(1: length(Lhloading{i}),3) = Lhloading{i}(:,3);
        LhloadingWhole(1: length(Lhloading{i}),4) = Lhloading{i}(:,4);
        LhloadingWhole(1: length(Lhloading{i}),5) = Lhloading{i}(:,5);
    else
        LhloadingWhole(startPoint(i-1):startPoint(i-1)+length(Lhloading{i})-1,1) = Lhloading{i}(:,1);
        LhloadingWhole(startPoint(i-1):startPoint(i-1)+length(Lhloading{i})-1,2) = Lhloading{i}(:,2);
        LhloadingWhole(startPoint(i-1):startPoint(i-1)+length(Lhloading{i})-1,3) = Lhloading{i}(:,3);
        LhloadingWhole(startPoint(i-1):startPoint(i-1)+length(Lhloading{i})-1,4) = Lhloading{i}(:,4);
        LhloadingWhole(startPoint(i-1):startPoint(i-1)+length(Lhloading{i})-1,5) = Lhloading{i}(:,5);
    end
end

end