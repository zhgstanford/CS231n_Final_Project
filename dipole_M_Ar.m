function d = dipole_M_Ar(Px, Py, Pz, config)
%d = dipole_M_Ar(Px, Py, Pz)
%This function calculates the dipole matrix element for atomic Ar. In the
%vicinity of the 3s->4p resoance. 
    
    E_AU = 27.2114; % eV to a.u.

    %Fano Parameters
    q = -0.389;
    GAMMA = 0.0874/E_AU;
    Ip = 15.7596/E_AU;
    E_R = 26.585/E_AU;
    
    %mesh of P^2 values
    P2 = Px.^2+Py.^2+Pz.^2;
    costh2 = Px.^2 ./ (P2);
        
    %Reduced Energy   
    epsilonMesh = (P2./2 - (E_R-Ip))/(GAMMA/2);

    %For 3p->nd
    AR_d = 10;
    AD_d = 0.75;
    A_d = ((q+epsilonMesh)./(1i+epsilonMesh))*AR_d + AD_d;

    %For 3p->ns
    AR_s = -1.7;
    AD_s = 4.4;
    A_s = ((q+epsilonMesh)./(1i+epsilonMesh))*AR_s + AD_s;

    %Angular Distributions
    Y_0 = 1/sqrt(4*pi);
    Y_2 = (sqrt(5/pi)/4) .* (3.* costh2 + ones(size(P2)));
    
    d = A_s .* Y_0 + A_d .* Y_2;
        

end