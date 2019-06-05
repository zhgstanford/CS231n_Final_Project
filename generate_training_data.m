%% Set up the script to generate the training data for the reconstruction of the ponderomtive energy

N_training = 1000; % The number of training data to be generated
% N_save = 100; % Save the training data after every N_save runs

Ip = 50; % eV
central_tau_x = 0.6; % fs
central_Up = 4*4/27.2;

N_p = 32; % Image size

P2D_train = zeros(N_training, N_p*N_p);
%Ip_train = zeros(N_training, 1);
Up_train = zeros(N_training, 1);
tau_x_train = zeros(N_training, 1);
phi_init_train = zeros(N_training, 1);

xi_train = zeros(N_training, 1);
tsep_train = zeros(N_training, 1);
dE_train = zeros(N_training, 1);
ratio_train = zeros(N_training, 1);
dphi_train = zeros(N_training, 1);
Q_train = zeros(N_training, 8*8);


%% Generate the training set

for n = 1:N_training 
    tau_x_fs = central_tau_x+0.2*randn;
    phi_init = 0.3*rand;
    
    xi = 1.5*rand;
    tsep = 1*randn*1e-15;
    dE = 0.8*randn; % eV
    ratio = 0.2*rand;
    %dphi = 2*pi*rand;
    dphi = 0;
    U_p = central_Up;
    
    if tau_x_fs < 0.1
        tau_x_fs = 0.1;
    end
    
    tau_x =tau_x_fs*1e-15;
    
    tau_x_train(n, 1) = tau_x;
    phi_init_train(n, 1) = phi_init;
    xi_train(n, 1) = xi;
    tsep_train(n ,1) = tsep;
    dE_train(n ,1) = dE;
    ratio_train(n, 1) = ratio;
    Up_train(n, 1) = U_p;
    
    [P2D, Q] = generate_train(N_p, tau_x, Ip, phi_init, xi, tsep, dE, ratio, dphi, U_p);
    P2D_train(n, :) = P2D(1,:);
    Q_train(n, :) = Q(:,1);
    fprintf("The %d-th shot has been finished! \n",n);
end

%% Save the data

savefile_time = floor(now);

P2D_name = 'P2D_train_'+string(savefile_time)+'.mat';
Up_name = 'Up_train_'+string(savefile_time)+'.mat';
tau_x_name = 'tau_x_train_'+string(savefile_time)+'.mat';
phi_init_name = 'phi_init_train_'+string(savefile_time)+'.mat';
xi_name = 'xi_train_'+string(savefile_time)+'.mat';
tsep_name = 'tsep_train_'+string(savefile_time)+'.mat';
dE_name = 'dE_train_'+string(savefile_time)+'.mat';
ratio_name = 'ratio_train_'+string(savefile_time)+'.mat';
Q_name = 'Q_train_'+string(savefile_time)+'.mat';
dphi_name = 'dphi_train_'+string(savefile_time)+'.mat';


save(P2D_name, 'P2D_train');
save(Up_name, 'Up_train');
save(tau_x_name, 'tau_x_train');
save(phi_init_name, 'phi_init_train');
save(xi_name, 'xi_train');
save(tsep_name, 'tsep_train');
save(dE_name, 'dE_train');
save(dphi_name, 'dphi_train');
save(ratio_name, 'ratio_train');
save(Q_name, 'Q_train');