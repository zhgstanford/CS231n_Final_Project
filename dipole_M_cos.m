function d = dipole_M_cos(Px, Py, Pz, config)
%d = dipole_M_cos(Px, Py, Pz)

    d = Px./sqrt(Px.^2+Py.^2+Pz.^2);

end