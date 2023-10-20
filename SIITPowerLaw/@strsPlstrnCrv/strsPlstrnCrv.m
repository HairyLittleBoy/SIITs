classdef strsPlstrnCrv
    %STRSPLSTRNCRV �˴���ʾ�йش����ժҪ
    %   �˴���ʾ��ϸ˵��
    
    properties 
        modelParas    % for power law (plStrain - stress): modelParas(1) - Elastic modulus
                                   %                      modelParas(2) - Yield stress
                                   %                      modelParas(3) - hardening index
                                   
                                   % for Voce law (plStrain - stress): modelParas(1) - Elastic modulus
                                   %                                   modelParas(2) - Yield stress
                                   %                                   modelParas(3) - A
                                   %                                   modelParas(4) - m                                  
    end
    
    methods
        function obj = strsPlstrnCrv(mlPs)
            obj.modelParas = mlPs;
        end
        
        strsPlstrnPoints = strsPlstrnPower(obj)    % generate plStrain - stress points using power law (plStrain - stress)
        
        strsPlstrnPoints = strsPlstrnVoce(obj)       % generate plStrain - stress points using Voce law (plStrain - stress)
    end
end

