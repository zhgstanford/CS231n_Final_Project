function [alpha, axis, vNgrid] = vNbasis(T, N_w, N_t, ndata)
%%%inputs
%T: total time range in femtoseconds
%N_w: number of vN omega lattice points
%N_t: number of vN time lattice points
%ndata: number of points in resampled axis, or time axis to be sampled.
%
%%%outputs
%alpha - struct():
%   alpha.t - time domain representation of basis functions, evaluated at
%   axis.t
%   alpha.w - frequency domain representation of basis functions, evaluated
%   at axis.w
%   alpha.t_sample - time domain representation of basis functions,
%   evaluated at axis.t_sample
%   alpha.w_sample - frequency domain representation of basis functions,
%   evaluated at axis.w_sample
%axis - struct():
%   axis.t - time domain representation time axis, N_t*N_w points spanning
%   (-T/2:T/2).
%   axis.w - frequency domain representation frequency axis, N_t*N_w points
%   spanning (-Omega/2:Omega/2); Omega = 2*pi * (N_t / T);.   
%   axis.t_sample - Resampled time axis, ndata points spanning (-T/2:T/2).
%   axis.w_sample - Resampled frequency axis, ndata points spanning 
%   (-Omega/2:Omega/2).
%vNgrid - vN lattice points.

N_basis=N_t*N_w; %Total number of vN lattice points

%vN basis functions
T_Basis = @(t,tn,wm,alpha) (2*alpha*pi).^(-1/4) .* exp( -(t-tn).^2 ./ (4 * alpha) ) .* exp(-1i * wm .* t);
W_Basis = @(w,tn,wm,alpha) (2*alpha/pi).^(1/4) .* exp( -alpha .* (w + wm).^2 ) .* exp( -1i * tn .* (w + wm) );

%Time Axis
Dt = T / N_basis; %Time step in time rep. 
t = (-T/2 + Dt/2): Dt: (T/2 - Dt/2); %Time axis in time domain

%Frequency Axis
%Omega = 2*pi * (N_t / T); %frequency range
Omega = 2*pi * (N_basis / T); %frequency range
DOmega = Omega / N_basis; %Frequency step in freq rep.

wmin=-Omega/2;
w = (wmin + DOmega/2): DOmega: (wmin + Omega - DOmega/2); %Frequency axis in frequency domain

%Resampling Axes
if numel(ndata) == 1 %If the input was the number of points to be sampled
    %Dts = T / ndata; %Time step for resampled data
    %t_sample = (-T/2 + Dts/2): Dts: (T/2 - Dts/2); %Resampled time axis
    t_sample = linspace(t(1),t(end),ndata); %Resampled time axis]
    %DOmegas = Omega / ndata; %Time step for resampled data
    %w_sample = (wmin + DOmegas/2): DOmegas: (wmin + Omega - DOmegas/2); %Resampled Omega axis
    w_sample = linspace(w(1), w(end), ndata); %Resampled Omega axis
else %If the input was a time vector to be sampled
    t_sample = ndata;
    T = max(ndata) - min(ndata);
    Os = 2*pi * (numel(ndata) / T);
    DOs = Os / numel(ndata) ;
    Oms = -Os /2;
    w_sample = (Oms + DOs/2): DOs: (Oms + Os - DOs/2);
end


%vN alpha
alpha_vN = T / (2*Omega); %in t^2

%vN lattice
dt = T / N_t; %Time step for vN grid
dOmega = Omega / N_w; %Freq step for vN grid
tn = (-T/2 + dt/2): dt: (T/2 - dt/2); %vN lattice points in time
wm = (wmin+dOmega/2): dOmega: (wmin+Omega-DOmega/2);%vN lattice points in frequency

%Basis Functions
alpha_t = zeros(N_basis); %time domain basis function
alpha_w = zeros(N_basis); %Frequency domain basis function
alpha_t_sample = zeros(N_basis, ndata); %Resampled time domain basis function
alpha_w_sample = zeros(N_basis, ndata); %Resampled frequency domain basis function

ind = 0;
for ind_t = 1:N_t
    for ind_w = 1:N_w
        ind = ind + 1;
        alpha_t(ind,:) = T_Basis(t, tn(ind_t), wm(ind_w), alpha_vN)';
        alpha_w(ind,:) = W_Basis(w, tn(ind_t), wm(ind_w), alpha_vN)';
        alpha_t_sample(ind,:) = T_Basis(t_sample, tn(ind_t), wm(ind_w), alpha_vN)';
        alpha_w_sample(ind,:) = W_Basis(w_sample, tn(ind_t), wm(ind_w), alpha_vN)';
    end
end

%%%%Build Output Structures
%Basis Fns
alpha.t = alpha_t;
alpha.w = alpha_w;
alpha.t_sample = alpha_t_sample;
alpha.w_sample = alpha_w_sample;
%Axes
axis.t = t;
axis.w = w;
axis.t_sample = t_sample;
axis.w_sample = w_sample;
%Grid
vNgrid = [tn;wm];

end