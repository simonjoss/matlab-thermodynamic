classdef VdW_EoS
    %EOS_ Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant)
        % gas constant in J mol⁻1 K⁻1
        R = 8.314;
    end

    properties (SetAccess = private,GetAccess = public)
        Tc
        Pc
        a
        b
    end
    
    methods
        function obj = VdW_EoS(Tc,Pc)
            obj.Tc = Tc;
            obj.Pc = Pc;
            obj.a = getA(obj,Tc,Pc);
            obj.b = getB(obj,Tc,Pc);
        end
        
        function a = getA(obj,Tc,Pc)
            a = 27 .* obj.R^2 .* Tc.^2 ./ (64 .* Pc);
        end
        
        function b = getB(obj,Tc,Pc)
            b = obj.R .* Tc ./ (8 .* Pc);
        end
        
        function p = P(obj,T,V)            
            p = obj.R .* T ./ (V - obj.b) - obj.a ./ V.^2;
        end
    end
end
        



