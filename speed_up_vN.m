%% SET UP
N_w=8;N_t=8;%vN lattice points in time and frequency
N_basis=N_w*N_t;
hbar=6.6e-16;%in ev*s
T=4.33e-15;%time range in s
energy_x=80;%xray energy in eV
N_p=32;%image size
I_A=sqrt(4*4/27.2);%streaking laser vector potential strength, in a.u.
Ip=50;%ionization potential in eV
tau_l=4.33e-15;%streaking laser period in s
[alpha,axis,vNgrid]=vNbasis(T, N_w, N_t, 1500);
alpha_t=alpha.t;
alpha_w=alpha.w;
alpha_tdata=alpha.t_sample;
alpha_wdata=alpha.w_sample;
tb=axis.t;
wb=axis.w;
tdata=axis.t_sample;
alpha_tdata=alpha_tdata./max(abs(alpha_tdata(:)));

%% Set up the streaking laser

t_p=tau_l*0.5;
t_X=tdata;%in s
A_L(1,:)=I_A.*cos(2*pi*(t_X(1,:)-t_p)/tau_l);
A_L(2,:)=I_A.*cos(2*pi*(t_X(1,:)-t_p)/tau_l+pi/2);

%% Run the nonlinear fitting algorithm
%setup config for streaking
config.Np = N_p;
config.Ip = Ip;
config.dipole_matrix = 'dipole_M_cos';
config.Kmax = 100/27.2;%Max kinetic energy in a.u.

global Bp_basis; Bp_basis=[];
global bp_basis; bp_basis=[];

M = reshape(P2D_train(1,:),[N_p N_p]);

[Qre,cost]=vNreconstruction_nlfit_au(M,N_w,N_t,A_L,energy_x,alpha_tdata,tdata,300,config);

%% We already have the Bp and bp. Now we can just do the nonlinear optimization directly

N_check = 1000;

grad_no_init = zeros(N_check,1);
grad_init_guess = zeros(N_check,1);
KL_no_init = zeros(N_check,1);
KL_init_guess = zeros(N_check,1);
cost_no_init = zeros(N_check,1);
cost_init_guess = zeros(N_check,1);

for ind = 1:N_check
    ind
    M = reshape(P2D_train(ind, :),[N_p N_p]);
    ground_truth_Q = Q_train(ind, :);
    [Qre1,cost1,grad1]=vNreconstruction_nlfit_au(M,N_w,N_t,A_L,energy_x,alpha_tdata,tdata,400,config); 
    [Qre2,cost2,grad2]=vNreconstruction_nlfit_CNN(M,N_w,N_t,A_L,energy_x,alpha_tdata,tdata,400,config, x0_initial(ind,:));
    
    grad_no_init(ind,1) = norm(grad1);
    grad_init_guess(ind,1) = norm(grad2);
    
    cost_no_init(ind, 1) = cost1;
    cost_init_guess(ind, 1) = cost2;
    
    true_Q_vN = abs(ground_truth_Q).^2/sum(abs(ground_truth_Q).^2);
    Qre1_vN = abs(Qre1).^2/sum(abs(Qre1).^2);
    Qre2_vN = abs(Qre2).^2/sum(abs(Qre2).^2);
    
    KL_no_init(ind,1) = KLDiv(Qre1_vN,true_Q_vN);
    KL_init_guess(ind,1) = KLDiv(Qre2_vN,true_Q_vN);
end

% %% Compare the results with and without the good initialization.
% [counts_1,centers_1] = hist(grad_no_init(:,1), 20);
% [counts_2,centers_2] = hist(grad_init_guess(:,1), 20);
% 
% grad_1 = sort(grad_no_init(:,1));
% grad_2 = sort(grad_init_guess(:,1));
% 
% subplot(2,1,1);
% hist(grad_1(1:135), 15);
% ylabel('Counts')
% xlabel('||\nabla_{c_n} loss||_2')
% title('Without Initialization')
% 
% subplot(2,1,2);
% hist(grad_2(1:135), 15);
% ylabel('Counts')
% xlabel('||\nabla_{c_n} loss||_2')
% title('With Initialization')
