classdef abaqusSimu < strsPlstrnCrv
    % generate L-h curve with abaqus
    %   
    
    properties 
        indenterR
        penetration
        niu        % poission's ratio
        miu       % fricCoef
        cpuNum
        meshSizeCoef    % meshSizeCoef can be set to 5, larger meshSizeCoef means finer mesh
    end
    
    methods
        function obj = abaqusSimu(mlPs,niu,indenterR,penetration,miu,cpuNum,meshSizeCoef)

            obj = obj@strsPlstrnCrv(mlPs);
            obj.niu = niu;
            obj.indenterR = indenterR;
            obj.penetration = penetration;
            obj.miu = miu;
            obj.cpuNum = cpuNum;
            obj.meshSizeCoef = meshSizeCoef;
        end

        jobName = jobNameGenr(obj)

        rewrtPyMain(obj)

        runPyMain(obj)
    end
end









