function d = dipole_M_H(Px, Py, Pz, config)
%out = dipole_M_H(Px, Py, Pz)
%This function calculates the dipole moment for ionization of a hydrogen
%atom. d = P ./ (P^2 + 2*Ip).^3;

    if isfield(config,'Ip')
        Ip = config.Ip / 27.2114; %Ip in a.u.
    else
        Ip = 1/2; % in au.
    end

    %mesh of P^2 values
    P2 = Px.^2+Py.^2+Pz.^2;

    d = abs(Px) ./ (P2 + 2*Ip).^3;
end