function [P2D, Q] = generate_train(N_p, tau_x, Ip, phi_init, xi, tsep, dE, ratio, dphi, U_p)
% N_p is the image size
% tau_x is the FWHM pulse duration of the XFEL pulse
% Ip is the ponderomotive energy of the streaking laser


% ratio: a small number between 0 and 0.3
% tsep: The separation between two pulses, on the order of 0.5 fs
% dE: The jitter in the pulse energy of the second Gaussian pulse

%% SET UP

N_w=8;N_t=8;%vN lattice points in time and frequency
hbar=6.6e-16;%in ev*s
T=4.33e-15;%time range in s
energy_x=80;%xray energy in eV

% Ponderomotive energy
%I_A=sqrt(4*4/27.2);%streaking laser vector potential strength, in a.u.
I_A = sqrt(U_p);


tau_l=4.33e-15;%streaking laser period in s
[alpha,axis,~]=vNbasis(T, N_w, N_t, 1500);
% alpha_t=alpha.t;
% alpha_w=alpha.w;
alpha_tdata=alpha.t_sample;
% alpha_wdata=alpha.w_sample;
% tb=axis.t;
% wb=axis.w;
tdata=axis.t_sample;
alpha_tdata=alpha_tdata./max(abs(alpha_tdata(:)));
%% constructed pulse(s)
chirp=1*(1+1i*xi);

time_jitter = phi_init*tau_l;

% Take the time jittering into our consideration
dt = mean(diff(tdata));
t_shift_jitter = floor(time_jitter/dt);

E_Xq1 = 2.^(-4*chirp*((tdata(1,:)+tsep/2)/sqrt(1+xi^2)/tau_x).^2);
E_Xq2 = ratio*2.^(-4*(1+1*1i*0)*((tdata(1,:)-tsep/2)/sqrt(1+1*0^2)/(tau_x*0.5*rand)).^2)*exp(-1i*dphi).*exp(-1i*dE/hbar.*tdata);

E_Xq = circshift(E_Xq1+E_Xq2, t_shift_jitter);

Q=(inv(alpha_tdata*alpha_tdata'))*(E_Xq*alpha_tdata')';%get Q
%set up streaking laser
%t_p=tau_l*0.5+phi_init*tau_l;%initial streaking phase, in s
%streaking laser vector potential. first row x, second row y.
t_p=tau_l*0.5;
t_X=tdata;%in s
A_L(1,:)=I_A.*cos(2*pi*(t_X(1,:)-t_p)/tau_l);
A_L(2,:)=I_A.*cos(2*pi*(t_X(1,:)-t_p)/tau_l+pi/2);

%setup config for streaking
config.Np = N_p;
config.Ip = Ip;
config.dipole_matrix = 'dipole_M_cos';
config.Kmax = 128/27.2;%Max kinetic energy in a.u.

%Run streaking calculation

energy_jitter = 1+0.03*randn;
[ P2D, ~, ~]=streak_au((E_Xq.*exp(-1i*energy_x*energy_jitter/hbar.*tdata)),A_L,tdata.*1e15,config);
end

