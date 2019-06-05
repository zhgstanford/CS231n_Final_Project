%% SET UP
N_w=8;N_t=8;%vN lattice points in time and frequency
N_basis=N_w*N_t;
hbar=6.6e-16;%in ev*s
T=4.33e-15;%time range in s
energy_x=80;%xray energy in eV
N_p=32;%image size
I_A=sqrt(4*4/27.2);%streaking laser vector potential strength, in a.u.
Ip=40;%ionization potential in eV
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
%% constructed pulse(s)
xi=sqrt(rand);
chirp=1*(1+1i*xi);
tau_x=.40e-15;%xray pulse duration in s, fwhm
tsep=2.165e-15;%separation between pulses in s
E_Xq=2.^(-4*chirp*((tdata(1,:)+0*tsep/2)/sqrt(1+xi^2)/tau_x).^2)+...
    0.3*2.^(-4*(1+1*1i*0)*((tdata(1,:)-tsep/2)/sqrt(1+xi^2)/(tau_x*rand)).^2);
Q=(inv(alpha_tdata*alpha_tdata'))*(E_Xq*alpha_tdata')';%get Q
Et=Q'*alpha_tdata.*exp(-1i*energy_x/hbar.*tdata);%efield in time domain
Et2=Q'*alpha_t;%efiled in time domain at vN pts
Ef=Q'*alpha_w;%efield in frequency domain at vN pts
Ef=Ef./max(abs(Ef));
phi_init=0;

%set up streaking laser
t_p=tau_l*0.5+phi_init*tau_l/2/pi;%initial streaking phase, in s
%streaking laser vector potential. first row x, second row y.
t_X=tdata;%in s
A_L(1,:)=I_A.*cos(2*pi*(t_X(1,:)-t_p)/tau_l);
A_L(2,:)=I_A.*cos(2*pi*(t_X(1,:)-t_p)/tau_l+pi/2);

%setup config for streaking
config.Np = N_p;
config.Ip = Ip;
config.dipole_matrix = 'dipole_M_cos';
config.Kmax = 100/27.2;%Max kinetic energy in a.u.

%Run streaking calculation
tic
[P2D,P3D,bp]=streak_au((E_Xq.*exp(-1i*energy_x/hbar.*tdata)),A_L,tdata.*1e15,config);
toc
M=reshape(abs(P2D),[N_p N_p]);

%% Make plots for the report

% Plot the input pulses

subplot(1,3,1);
yyaxis left;
plot(tdata(1,500:1200)*1e15, (abs(E_Xq(1,500:1200))/max(abs(E_Xq))).^2);
hold on;
yyaxis left;
scatter(tb(1,22:50)*1e15, (abs(Et2(1,22:50))/max(abs(Et2))).^2, '*');
hold on;
ylabel('Intensity (a. u.)')

yyaxis right;
plot(tdata(1,500:1200)*1e15, unwrap(angle(E_Xq(1,500:1200))));
hold on;
yyaxis right;
scatter(tb(1,22:50)*1e15, unwrap(angle(Et2(1,22:50)))+0*pi, '*');
ylabel('Phase (rad.)')

xlabel('t (fs)');

title('Temporal XFEL Pulse');

subplot(1,3,2);
imagesc(vNgrid(1,:)*1e15, hbar*vNgrid(2,:), abs(reshape(Q, [8,8])).^2/(max(abs(Q)))^2);
xlabel('t (fs)')
ylabel('\Delta \omega (eV)')
colorbar;

title('von Neuman Coefficients')

subplot(1,3,3);
imagesc(M/max(max(M)));
colorbar;
xlabel('p_x');
ylabel('p_y');
title('VMI Image')

%% Reconstruction
global Bp_basis; Bp_basis=[];
global bp_basis; bp_basis=[];
%or load saved basis functions if already created
%[Qre,cost]=vNreconstruction_nlfit(M,N_w,N_t,tau_l,Ip,I_A,energy_x,alpha_tdata,tdata,300);
[Qre,cost]=vNreconstruction_nlfit_au(M,N_w,N_t,A_L,energy_x,alpha_tdata,tdata,100,config);
%% plot reconstruction
reE=Qre*alpha_t;
reEf=Qre*alpha_w;



Et2=Et2./max(abs(Et2));
origI=Et2.*conj(Et2);
reE=reE./max(abs(reE));
reI=reE.*conj(reE);
phi1=unwrap(angle(Et2));
phi2=unwrap(angle(reE));
phi2 = phi2+ phi1(1,32)-phi2(1,32);
idx0=abs(Et2).^2>max(abs(Et2).^2)*0.2;
Eterr=sqrt(sum(((abs(abs(Et2(idx0))-abs(reE(idx0))))).^2)./sum(idx0));
disp(['cost=' num2str(cost) ', Eterr=' num2str(Eterr)]);
figure;yyaxis left;plot(tb*1e15,origI,'LineWidth',2);
yyaxis right;plot(tb*1e15,phi1,'LineWidth',2);
yyaxis left;hold on;plot(tb*1e15,reI,'--*','LineWidth',1.5,'MarkerSize',8);xlim([-1.5 1.5]);
yyaxis right;hold on;plot(tb*1e15,phi2,'--*','LineWidth',1.5,'MarkerSize',8);
xlim([-1.5 1.5]);
set(gcf,'color','w');set(gca,'fontsize', 30);
xlabel('t (fs)');
yyaxis left;ylabel('intensity (a.u.)');
yyaxis right;ylabel('phase (rad.)');
legend('input intensity','input phase','reconstructed intensity','reconstructed phase');


origS=Ef.*conj(Ef);
reEf=reEf./max(abs(reEf));
reS=reEf.*conj(reEf);
phi1=unwrap(angle(Ef));
phi2=unwrap(angle(reEf));
phi2 = phi2+ phi1(1,32)-phi2(1,32);
figure;yyaxis left;plot(wb*hbar+energy_x,origS,'LineWidth',2);
yyaxis right;plot(wb*hbar+energy_x,phi1,'LineWidth',2);
yyaxis left;hold on;plot(wb*hbar+energy_x,reS,'--*','LineWidth',1.5,'MarkerSize',8);
yyaxis right;hold on;plot(wb*hbar+energy_x,phi2,'--*','LineWidth',1.5,'MarkerSize',8);
xlim([60 100]);
set(gcf,'color','w');set(gca,'fontsize', 30);
xlabel('\omega (eV)');
yyaxis left;ylabel('intensity (a.u.)');
yyaxis right;ylabel('phase (rad.)');
legend('input intensity','input phase','reconstructed intensity','reconstructed phase');