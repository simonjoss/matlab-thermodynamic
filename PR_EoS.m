classdef PR_EoS
    properties (SetAccess = private, GetAccess = public)
        % critical temperature (K)
        Tc
        % critical pressure (Pa)
        Pc
        % peng robinson parameter
        b
        % peng robinson parameter, kappa
        k
        % molar wheight kg mol^-1
        m
    end
    
    properties (Constant)
        % gas constant in J mol⁻1 K⁻1
        R = 8.314;
    end
    
    properties (Access = private)
        % poly coefficients for the heat capacity
        CpCoeff@double
        % refrrence temperature
        Tref@double
        % reference pressure
        Pref@double
    end
    
    methods
        function obj = PR_EoS(Tc,Pc,omega,cpCoeff,Tref,Pref,m)
            obj.Tc = Tc;
            obj.Pc = Pc;
            obj.b = obj.getB(Tc,Pc);
            obj.k = obj.kappa(omega);
                        
            % add zero for define the a constant value of zero
            obj.CpCoeff = cpCoeff;
            obj.Tref = Tref;
            obj.Pref = Pref;
            obj.m = m;
        end
        
        function p = P(obj,T,VMolar)
            a = get_a(obj,T);            
            p = obj.R .* T ./ (VMolar - obj.b) - a ./ (VMolar .* (VMolar + obj.b) + obj.b .* (VMolar - obj.b));
        end
        
        function b = getB(obj,Tc,Pc)
            % b parameter of Peng-Robinson EoS (page 250, eq 6.7-2)
            b = 0.07780 * obj.R * Tc / Pc;
        end
        
        function Z = solveRoots(obj,T,P)
            % solve the poly roots
%             a = obj.get_a(T);
            
            % peng-robinson parameters
                       
            B = get_B(obj,T,P);
            A = get_A(obj,T,P);
            %A = a .* P ./ (obj.R .* T).^2;
            
            % define the coefficients of the poly
            p = [...
                    1, ...
                    -1 + B,...
                    A - 3 .* B.^2 - 2 .* B,...
                    -A .* B + B.^2 + B.^3];
            
            r = roots(p);
            
            if ~isreal(r)
                for ii = numel(r) :-1:1
                    if isreal(r(ii))
                        Z(ii) = real(r(ii));
                    end
                end
            else
                % the root in between is unstable. thus this root is not of
                % interest
                Z = [min(r),max(r)];
                Z(Z < 1e-12) = [];
            end
        end
        
        function showRoots(obj,T,P)
            % solves the roots and plots the results
            Z = obj.solveRoots(T,P);
            
            % calculate the molar volume for each root
            VMolarAtRoot = obj.VMolar(T,P,Z);
            
            dVMolarStep = min(abs(VMolarAtRoot))*1e-2;
            % calculate isotherm line
            VMolar = 0.5 * min(VMolarAtRoot) : dVMolarStep : 2 * max(VMolarAtRoot);
            
            p = obj.P(T,VMolar);
            
            figure()
            plot(VMolar,p.*1e-6,'b-')
            box on
            grid on
            hold on
            
            plot(VMolarAtRoot,P .* ones(numel(Z)) .* 1e-6,'x','MarkerSize',8)
            
            ylabel('Pressure (MPa)')
            xlabel('Volume (m^3 mol^-1)')
            title(['Peng Robinson EOS, Isotherm = ' num2str(T) 'K'])
            axis([min(VMolar),max(VMolar),min(p)*1e-6,max(p)*1e-6])
            
            % plot z value 
            for ii = numel(Z) : -1 : 1
                text(VMolarAtRoot(ii),P .* 1e-6,['Z = ' num2str(Z(ii))]);
            end
            
            line([0,max(VMolar)],[P,P].*1e-6)
        end
        
        function [HDepFun,SDepFun,f] = departureFun(obj,T,P,Z)
            % returns the departure functions of the enthalpy and entropy
            ac = obj.get_a(T);
            a = alpha(obj,T);
            
            % deriviate dadT (page 251)
            dadT = -0.45724 .* (obj.R .* obj.Tc).^2 ./ obj.Pc .* obj.k .* sqrt(a ./ (T .* obj.Tc));
            
            % coefficient B
            B = P .* obj.b ./ (obj.R .* T);
            
            C = log(  ( Z + (1 + sqrt(2)) .* B) ./ ( Z + (1 - sqrt(2)) .* B)  );
            
            HDepFun = obj.R .* T .* (Z - 1) + (T .* dadT - ac) ./ (2 * sqrt(2) * obj.b) .* C;
            SDepFun = obj.R .* log(Z - B) + dadT ./ (2 .* sqrt(2)) .* C;
            
            A = a .* P ./ (obj.R .* T).^2;
            
            f = P .* exp((Z - 1) - log(Z - B) - A ./ (2 .* sqrt(2) .* B) .* C);
        end
        
        function H = enthalpy(obj,T,P,Z)
            % returns the enthalpy in J mol^-1
            cp = obj.CpCoeff;

            [HDepFun,~] = departureFun(obj,T,P,Z);
            
            HIG = (cp(1) / 4) .* (T.^4 -obj.Tref.^4) +...
                    (cp(2) / 3) .* (T.^3 -obj.Tref.^3) + ...
                    (cp(3) / 2) .* (T.^2 -obj.Tref.^2)+ ...
                    (cp(4) / 1) .* (T - obj.Tref); 
            
            H = HDepFun + HIG;
        end
        
        function S = entropy(obj,T,P,Z)
            % returns the entropy in J mol^-1 K^-1
            cp = obj.CpCoeff;
            
            % get the departure function of the ideal gas
            [~,SDepFun] = departureFun(obj,T,P,Z);
            
            SIG = cp(4) .* log(T ./ obj.Tref) + ...
                    cp(3) .* (T - obj.Tref) + ...
                    cp(2) ./ 2 .* (T.^2 - obj.Tref.^2) + ...
                    cp(1) ./ 3 .* (T.^3 - obj.Tref.^3) + ...
                    - obj.R .* log(P ./ obj.Pref);
                
            S = SDepFun + SIG;
        end
        
        function V = VMolar(obj,T,P,Z)
            %returns the molar volume in mol m^-3
            V = obj.R .* Z .* T ./ (P);
        end
        
        function U = internalEnergy(obj,T,P,Z)
            V = VMolar(obj,T,P,Z);
            
            H = enthalpy(obj,T,P,Z);
            
            % internal energy
            U = H - P .* V;
        end
        
        function [f,phi] = fugacity(obj,T,P,Z)
            [~,~,f] = departureFun(obj,T,P,Z);
            
            % fugacity coefficient phi
            phi = f ./ P;
        end
        
        function G = gibbs(obj,T,P,Z)
            H = enthalpy(obj,T,P,Z);
            S = entropy(obj,T,P,Z);
            
            G = H - T .* S;
        end
    end
    
    methods (Static)
        function k = kappa(omega)
            % kappa for Peng-Robinson EoS (page 250 eq 6.7-4)
            k = 0.37464 + 1.54226 .* omega - 0.26992 .* omega.^2;
        end
    end
    
    methods (Access = private)
        function a = alpha(obj,T)
            % alpha parameter for Peng-Robinson EoS (page 250, eq 6.7-1)
            a = (1 + obj.k .* (1 - sqrt(T ./ obj.Tc))).^2;
        end
        
        function A = get_A(obj,T,P)
            a = get_a(obj,T);
            % A parameter for Peng-Robinson EoS (page 251)
            A = a .* P ./ (obj.R .* T).^2;            
        end
        
        function B = get_B(obj,T,P)
            % B parameter for Peng-Robinson EoS (page 251)
            B = P .* obj.b ./ (obj.R .* T); 
        end
        
        function a = get_a(obj,T)
            % a parameter for Peng-Robinson EoS (page 250, eq 6.7-3)
            a = obj.alpha(T);
            a = 0.45724 .* obj.R^2 .* obj.Tc^2 ./ obj.Pc .* a;            
        end
    end
end