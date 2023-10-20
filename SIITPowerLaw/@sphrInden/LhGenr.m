function Lhdata = LhGenr(obj)
%LhGenr generate a Lh curve with the method which property 'method'defines
switch obj.LhGenrMethod
    case 'FEM'
        jobName = jobNameGenr(obj); % defined in clasee abaqusSimu
        runPyMain(obj)              % defined in clasee abaqusSimu
        Lhdata= LhImport(['Lh_',jobName,'.txt']);
    case 'ANN'
        Lhdata = LhCalcAnn(obj);
    case 'interp'
        Lhdata = LhCalcNumr(obj);
end
end