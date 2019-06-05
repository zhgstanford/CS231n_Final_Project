function [Q,cost,grad]=vNreconstruction_nlfit_au(M,N_w,N_t,A_L,energy_x,alpha_tdata,tdata,iterMax,config)
hbar=6.6e-16;%in ev*s
N_basis=N_w*N_t;
N_p=config.Np;
%%%calculate basis:bp
global bp_basis
if size(bp_basis)==0
    %%%setup vN basis and axes
    alpha_tdata=alpha_tdata.*exp(-1i.*energy_x/hbar.*tdata);
    bp_basis=make_bp(alpha_tdata,tdata,A_L,config,N_basis);%dimension: [n*m, Px*Py, Pz]
end
%%%calculate basis: Bp
global Bp_basis
if size(Bp_basis)==0
    Bp_basis=make_Bp(bp_basis,N_basis,N_p);%dimension: [(n*m)^2, Px*Py]
end

%%%nonlinear fitting
global P2D
P2D=M;
x_initial=rand(1,2*N_basis-1);
tic
fun=@Complex_Residual_Error_Nonlinear_Slope;
options = optimoptions(@fminunc, 'MaxIterations', iterMax, 'Algorithm','quasi-newton','SpecifyObjectiveGradient',true,'CheckGradients',false);
[xout,cost,~,~,grad,~] = fminunc(fun, x_initial, options);
toc

ReCo=xout(1,1:N_basis);
ImCo=[xout(1,(N_basis+1):(ceil(1.5*N_basis)-1)),0,xout(1,ceil(1.5*N_basis):(2*N_basis-1))];
Q=ReCo+1i*ImCo;
end

function [bp_basis]=make_bp(E,t,A_L,config,N_basis)
N_p=config.Np;
bp_basis=zeros(N_basis,N_p*N_p,N_p);
for i = 1:N_basis
    [~,~,bp]=streak_au((E(i,:)),A_L,t*1e15,config);
    bp_basis(i,:,:)=transpose(bp);%bp_basis dimension: [n*m, Px*Py, Pz]. edited on 5/3/2018
end
end

function Bp_basis=make_Bp(bp_basis,N_basis,N_p)
Bp_basis=zeros(N_basis^2,N_p*N_p);
for i=1:N_basis
    for j=1:N_basis
        index=(i-1)*N_basis+j;
        Bp=bp_basis(j,:,:).*conj(bp_basis(i,:,:));%dimension: [1,Px*Py,Pz]
        Bp_basis(index,:)=trapz(Bp,3);%dimension: [(n*m)^2, Px*Py]
        if i==j
            figure(77);imagesc(reshape(abs(Bp_basis(index,:)),[N_p N_p]));drawnow;
        end
    end
end
end

function [ error, gradient ] = Complex_Residual_Error_Nonlinear_Slope(X)
% For the complex coefficients.

global Bp_basis;
global P2D
N_basis=(numel(X)+1)/2;
N_p=size(P2D,1);

middle = ceil(1.5*N_basis);
ReCo = X(1,1:N_basis);
ImCo = [X(1,(N_basis+1):(middle-1)),0,X(1,middle:(2*N_basis-1))];

Qguess=ReCo+1i*ImCo;
Qguess=Qguess(:);
diff=real(reshape(Qguess*Qguess',[1 N_basis^2])*Bp_basis)-transpose(P2D(:));
error=sum(abs(diff).^2);

%The following part calculates the gradiant array
if nargout>1
    g=zeros(N_basis,N_p^2);%partial derivative of Q
    for j=1:N_basis
        index=((j-1)*N_basis+1):j*N_basis;
        g(j,:)=2*(transpose(Qguess))*Bp_basis(index,:);%factor of two because \partial Q/\partial Q*=1 (NOT 0!)
    end
    
    Slope_Re=(2*sum(diff.*real(g),2))';
    Slope_Im=(2*sum(diff.*imag(g),2))';
    
    gradient = [Slope_Re, Slope_Im(1,1:(middle-1-N_basis)), Slope_Im(1,(middle+1-N_basis):N_basis)];
end

end