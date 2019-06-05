function [P2D,P3D,bp]=streak_au(E_X,A,t_X,config)
%[P2D,P3D,bp] = streak_au(E_X,A,t_X,config)
%This function calculates the probability amplitude, b(p) for observing an 
%electorn with momentum p after the system is ionized by an XUV pulse (E_X) 
%in the presence of an IR laser field (A). The calculation is done in the 
%strong-field approximation (SFA).
%
%%%%%%%% inputs:
%E_X: input electric field, with the carrier phase, ie. includes the term
%       exp(-1i*omega_xray*t)
%A: vector potential, A(1,:) = Ax, A(2,:) = Ay 
%t_X: array of time points in fs. Assume time points are equally spaced.
%config - struct():
%config.Ip: ionization potential in eV [default = 13.6 eV]
%config.Np: number of momentum points to use for calculation. [default = 64]
%config.dipole_matrix: function name for computing dipole moment.
%                       [default is dipole_M_H] 
%config.Kmax: The maximum kinetic energy (in a.u.). The largest momentum is then 
%                       (2*K_max/m_e)^0.5. [Default = 1 a.u.]
%config.Kmin: The minimum kinetic energy (in a.u.). Probability for
%                       energies below Kmin are set to zero. [Default = 0]
%config.Pz_Slice: Allows for the calculation of a single P_z slice of the
%                       photoelectron momentum distribution. When using this 
%                       option, function returns P3D = 0, and P2D returns 
%                       the desired slice. To remove this field and 
%                       calcuate the full distribution, use: config = rmfield(config,'Pz_Slice'). 
%%%%%%%% outputs:
%P2D, P3D: 2D and 3D momentum distribution, i.e. |bp|^2
%size(P3D) = [config.Np, (config.Np)^2]
%size(P2D) = [1, (config.Np)^2]
%bp: transition amplitude, size(bp) = [config.Np, (config.Np)^2]

%Check that config fields exist. 
if (~isfield(config,'dipole_matrix'))
    config.dipole_matrix = 'dipole_M_H';
end
if (~isfield(config,'Ip'))
    config.Ip = 13.6;
end
if (~isfield(config,'Np'))
    config.Np = 64;
end
if(~isfield(config,'Kmax'))
    config.Kmax = 1;
end
if(~isfield(config,'Kmin'))
    config.Kmin = 0;
end
    
%E_X, A, and t_x should have the same length
if ( (numel(t_X) == numel(E_X)) && ( numel(t_X) == size(A,2) ) )
    
    E_X = reshape(E_X, numel(E_X), 1);
    bp = Photoelectron_Momentum_Distribution(E_X,A,t_X,config);

    if( isfield(config,'Pz_Slice') )
        P2D = abs(bp).^2; %dimension [1,px*py]
        P3D = 0;
    else
        P3D = abs(bp).^2; %dimension [pz,px*py]
        P2D = trapz(P3D,1);%dimension [px*py,1]
    end
else
    sprintf('Error: Matrix size is incosistent')
    P3D = 0;
    P2D = 0;
    bp = 0;
end


end

function  b_p = Photoelectron_Momentum_Distribution(E_X, A, t_X, config )
%
% 
%Constants used in calculations
T_AU = 0.024189; %fs to a.u.
E_AU = 27.2114; %eV to a.u.

%Apply Configurations and convert units
I_p = config.Ip / E_AU; %Change to a.u.
K_max = config.Kmax;
K_min = config.Kmin;
N_p = config.Np/2; 
dipoleM = str2func(config.dipole_matrix);
t_X = t_X / T_AU; % Change fs to a.u.

N_t = numel(t_X); %Number of time points
T=t_X(end)-t_X(1);%time window in a.u.
dt = T / ( numel(t_X) - 1 ); %Time step in a.u.

%Setup momentum grid
P = sqrt(2*K_max);  % The largest momentum.
dp = (2*P)/(2*N_p);  % Momentum step

P_xy = -P+dp/2:dp:P-dp/2; % vector for x,y momentum
%vector for z momentum
if ( isfield(config,'Pz_Slice') )
    P_z = config.Pz_Slice;
else
    P_z = -P+dp/2:dp:P-dp/2;
end

[P_xm,P_ym] = meshgrid(P_xy); %Mesh Grid for momentum
%P_xm = reshape(P_xm,1,(2*N_p)^2); %Reshape momentum to a vector.
%P_ym = reshape(P_ym,1,(2*N_p)^2);

N_pz = numel(P_z); % number of z-points
N_pi = (2*N_p)^2; % number of x,y points

%Prealocate memory for b_p 
b_p = zeros(N_pz,N_pi); % size should be [N_pz,(N_px*N_py)]

%V(t)=p-A(t); %Canonical Momentum, used for calculating the Action
V_x = ones(N_t,1) * P_xm(:)' - A(1,:)' * ones(1,N_pi);
V_y = ones(N_t,1) * P_ym(:)' - A(2,:)' * ones(1,N_pi);

I_xy = (1/2) * ( abs(V_x).^2 + abs(V_y).^2 ); % Integrand for Action

%S_xy(t) = int(t,T,I_xy)
%S_xy(t)  =       int(0,T,I_xy)      -    int(0,t,I_xy)
%S_xy = ones(N_t,1) * sum ( I_xy , 1) - cumsum( I_xy , 1);
%Phase_xy = exp(-1i * (dt) .* S_xy); %Phase part of time integral

%S_xy(t) = int(T,t,I_xy)
S_xy = cumsum( I_xy , 1, 'reverse') - I_xy;
Phase_xy = exp(-1i * (dt) .* S_xy ); %Phase part of time integral 

for ind_z = 1:N_pz
    
    p = P_z(ind_z);
    K_z = (1/2) * p^2;
    
    if (K_min>0) % Only do indexing if neccessary, speeds up calc.
        ind_p = find( (P_xm(:).^2 + P_ym(:).^2 + p.^2) > K_min );
        N_i = numel(ind_p);
    
        Phase_z = exp(-1i * (K_z + I_p) .* (t_X(end) - t_X') ) * ones(1,N_i);
        
        %Calculate dipole d(p_x-A_x,p_y-A_y,p_z)
        dipole = dipoleM(V_x(:,ind_p),V_y(:,ind_p),p,config);
        
        %Int[E_X*d * Exp(phi_xy) * Exp(phi_z)]
        b_p(ind_z,ind_p) = dt * sum( ( E_X * ones(1,N_i) ) .* dipole .* Phase_xy(:,ind_p) .* Phase_z, 1 );
            
    else   
        Phase_z = exp(-1i * (K_z + I_p) .* (t_X(end) - t_X') ) * ones(1,N_pi);

        %Calculate dipole d(p_x-A_x,p_y-A_y,p_z)
        dipole = dipoleM(V_x,V_y,p,config);

        %Int[E_X*d * Exp(phi_xy) * Exp(phi_z)]
        b_p(ind_z,:) = dt * sum( ( E_X * ones(1,N_pi) ) .* dipole .* Phase_xy .* Phase_z, 1 );
    end
    
end
    
end

function d = dipole_M_cos(Px, Py, Pz, config)
%d = dipole_M_cos(Px, Py, Pz)

    d = Px./sqrt(Px.^2+Py.^2+Pz.^2);

end