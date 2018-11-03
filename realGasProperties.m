function result = realGasProperties(T,Tc,P,Pc,omega,cp,m,varagin)
    import PR_EoS
    
    if T <= 0 || Tc <= 0
        error('realGas:badInput',...
                'Temperature must be greater than zero!')
    end
    
    if P <= 0 || Pc <= 0
        erorr('realGas:badInput',...
                'Pressure must be greater than zero!')
    end
    
    % refrence temperature and pressure for the calculation of the
    % properties
    Tref = 25 + 273.15;
    Pref = 1e5;
    
    % load the peng robinson equations
    EoS = PR_EoS(Tc,Pc,omega,cp,Tref,Pref,m);
    
    
    Z = solveRoots(EoS,T,P);
            
    H = enthalpy(EoS,T,P,Z);
    
    S = entropy(EoS,T,P,Z);
    
    G = gibbs(EoS,T,P,Z);
    
    U = internalEnergy(EoS,T,P,Z);
    
    V = VMolar(EoS,T,P,Z);
    
    [f,phi] = fugacity(EoS,T,P,Z);
    
    fprintf('-+-+-+-+-+-+-+-+-Peng Robinson-+-+-+-+-+-+-+-\n')
    fprintf('Presure P (MPa): \t\t %2.2f \n',P.*1e-6)
    fprintf('Temperature T (K): \t\t %2.2f \n',T)
    fprintf('Molar volume V (m^3/mol): \t %2.4f \n',V)
    fprintf('Entropy S (J/(mol K)): \t %2.2f \n',S)
    fprintf('Enthalpy H (J/mol): \t\t %2.2f \n',H);
    fprintf('Gibbs G (J/mol): \t\t %2.2f \n',G);
    fprintf('Internal Energy U (J/mol): \t %2.2f \n',U);
    fprintf('Fugacity f (MPa): \t\t %2.4f \n',f .* 1e-6);
    fprintf('Fugacity coefficient phi: \t %2.2f \n',phi)
    fprintf('-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-\n')
    
    % fill out the return struct
    result.Z = Z;
    result.H = H;
    result.S = S;
    result.G = G;
    result.U = U;
    result.V = V;
    result.f = f;
end