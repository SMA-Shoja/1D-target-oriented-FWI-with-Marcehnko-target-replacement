% 1D target-oriented FWI with Marchenko target-replacement
% Author: Aydin Shoja Email: se.mo.ay.sh@gmail.com; s.m.a.shoja@tudelft.nl
clear; close all; clc;

%% Initialization

full_band = 0;                                                              % full_band = 1 means using the full frequency spectrum and 0 means using the wavelet spectrum and multi-scale method.
optim = 1;                                                                  % Optimization method, Gradient decent = 0 and Gauss-Newton = 1. Gausse-Newton only works with full_band = 0.

nt = 2000;                                                                  % Number of time samples
dt = 1e-2;                                                                  % Time sampling in seconds
t = nt*dt;                                                                  % recording time

k = 1:nt;                                                                   % Time samples vector
time = dt.*(k);                                                             % Time vector

fs = 1/dt;                                                                  % Maximum allowed frequency
df = 1/(t);                                                                 % Frequency sampling in Hz
Nyq_fq = 1/(2*dt);                                                          % Nyquist frequency

if mod(nt,2) == 0                                                           % Nt is even
    fvec = -Nyq_fq : df : Nyq_fq;
else                                                                        % Nt is odd
    fvec = [sort(-1*(df:df:Nyq_fq)) (0:df:Nyq_fq)];                         % Frequency vector
end
fvec(end) = [];
fvec(fvec==0) = eps;
w = 2*pi*fftshift(fvec)';                                                   % Angular frequency vector
dw = 2*pi*df;


%% model definition

% Model with two velocity anomalies

realv = [2500*ones(1,100),1500*ones(1,100),2500*ones(1,300)...
    ,2700*ones(1,100),2500*ones(1,50)]';                                    % true velocity model
refv = realv;
refv(210:end,1) = 2500;
refv = velsmooth(refv,5,5,150);                                             % initial velocity
  a = refv(1:450,1);                                                        % overburden smooth velocity for Marchenko method
  b = refv(451:end,1);                                                      % target initial velocity for FWI

A = realv(1:450,1);                                                         % true overburden velocity
B = realv(451:end,1);                                                       % true target velocity

ro_b = ones(size(b));                                                       % target region density vector
ro_a = ones(size(a));                                                       % overburden density vector
ro_ref = ones(size(realv));                                                 % initial model density

dz=10;                                                                      % Depth sampling in meters
nzf = length(a);                                                            % Focusing depth index
zf = dz*nzf;                                                                % Focusing depth in m

z = dz*length(realv);                                                       % Depth vector in m
nz = length(realv);                                                         % Number of depth samples


z_mod = dz*ones(nz,1);                                                      % Vector containing depth samples
z_modb = dz*ones(length(b),1);                                              % Vector containing depth samples of target
zvecr = dz*(1:nz);
zvec = dz*(nzf+1:nz);                                                       % Depth vector
W = repmat(w,1,nz);                                                         % Angular frequency matrix
%% source definition
f0 =1;                                                                      % peak frequency of the wavelet
[wav]	=	ricker(f0,nt,dt);                                               % Ricker wavelete


%% Computing the data
np	=	1;
dp	=	0.0002;
p0	=	0;
norm = 0;
mul = 1;
nprim = 0;
kzf = w./refv(nzf);                                                         % wavenumber vector at focusing depth
kz0 = w./refv(1);                                                           % wavenumber vector at surface
[~,data]=layercac(realv,ro_ref,dz,nt,dt,np,dp,p0,norm,mul,nprim);           % making data with true model

dataw = real(ifft(fft(data) .* fft(wav)));                                  % convolving data with wavelete

%% Marchenko method
[Td,~,F1p0,~] = layercac(A,ro_a,dz,nt,dt,np,dp,p0,norm,0,0);
Tdinv = F1p0;
[win]=tdeps(Tdinv,nt);
f1m	=	flipud(win).*real(ifft(fft(data).*fft(Tdinv)));
Mp	=	win.*real(ifft(conj(fft(data)).*fft(f1m)));
for iter = 1:20
    f1m	=	flipud(win).*real(ifft(fft(data).*fft(Tdinv+Mp)));
    Mp	=	win.*real(ifft(conj(fft(data)).*fft(f1m)));
end
f1p=Tdinv+Mp;
f2m=f1p;
f2p=-real(ifft(conj(fft(f1m))));
TA	=	real(ifft((fft(f1p)).^(-1)));
RAdown = real(ifft(((fft(f2m)).^(-1)).*fft(f2p)));
RAup	=	real(ifft(((fft(f1p)).^(-1)).*fft(f1m)));
%% Computing target response
%inserting the response of b into the response of A
[Tbnew,Rbnewup,~,~]=layercac(b,ro_b,dz,nt,dt,np,dp,p0,norm,mul,0);
% inverting equation 42
GBbp 	=	real(ifft(((ones(nt,1)-fft(RAdown).*fft(Rbnewup)).^(-1)).*fft(TA)));
% equation 40
calc	=	RAup + real(ifft(fft(TA).*fft(Rbnewup).*fft(GBbp)));

%% making the green's functions of initial target region

[Gbb]=vsplayrc(b,ro_b,dz,nt,dt,p0,norm);                                    % isolated target region Green's function
Gbb = real(ifft((fft(Gbb))));
Gs = real(ifft(fft(Gbb).*repmat(fft(GBbp),1,length(b))));                   % convolving Gbb with overburden interactions
Gsw = repmat(fft(wav),1,nz-nzf).*fft(Gs);
Gswnew = real(ifft(Gsw));                                                   % Green's function convolved with wavelet



%% FWI loop
iteration = 20;
Gbbm = zeros(nt,length(b));
Gbbp = Gbbm;
Gsnew = Gs;
newv = refv(nzf+1:nz);
green = zeros(nt,nz-nzf,iteration);
B = (-2*1i*kz0);
Ha = zeros(nz-nzf,nz-nzf,iter);
    
for iter = 1:iteration
    
    iter
    
    
    green(:,:,iter) = Gsnew;
    newv2 = repmat(newv,1,nt)';
    Kz = W(:,nzf+1:end)./newv2;
    D = (-2*1i*Kz);

  if full_band == 1  

    res = calc - data;    
    cost(1,iter) = res'*res;  
    cost_norm = cost./cost(1,1);
    pred(:,iter) = calc;
    residual(:,iter) = res;

    Grf = fft(Gsnew)./D;                                                    % Green's functions of reciever's location
    Gsf = fft(Gsnew)./D;                                                    % Green's functions at source's location
    resf = fft(res)./B;

    % data visualization
    figure(1)
    subplot(3,1,1)
    plot(time,data);
    title('Data');
    xlabel(' Time(s) ')
    ylim([min(data) max(data)])
    pause(0.1)
    subplot(3,1,2)
    plot(time,calc);
    title('Predicted Data');
    xlabel(' Time(s) ')
    ylim([min(data) max(data)])
    pause(0.1)
    subplot(3,1,3)
    plot(time,res);
    title('residual');
    xlabel(' Time(s) ')
    ylim([min(data) max(data)])
    pause(0.1)

  elseif full_band == 0  
    calcw = real(ifft(fft(calc) .* fft(wav))); 
    res = calcw - dataw;
    pred(:,iter) = calcw;
    residual(:,iter) = res;
    cost(1,iter) = res'*res;  
    cost_norm = cost./cost(1,1);

    Grf = fft(Gsnew)./D;                                                      
    Gsf = fft(Gswnew)./D;                                                   
    resf = fft(res)./B;
    % Multi-scale condition
    if (cost_norm(1,iter) < 0.15) && (f0 < 50)
         f0 = f0+5;                                                                    
         [wav]	=	ricker(f0,nt,dt);
         dataw = real(ifft(fft(data) .* fft(wav))); 
         calcw = real(ifft(fft(calc).*fft(wav)));
         res = calcw - dataw;
         cost(1,iter) = res'*res;  
    end
   % data visualization
    figure(1)
    subplot(3,1,1)
    plot(time,dataw);
    title('Data');
    xlabel(' Time(s) ')
    ylim([min(dataw) max(dataw)])
    pause(0.1)
    subplot(3,1,2)
    plot(time,calcw);
    title('Predicted Data');
    xlabel(' Time(s) ')
    ylim([min(dataw) max(dataw)])
    pause(0.1)
    subplot(3,1,3)
    plot(time,res);
    title('residual');
    xlabel(' Time(s) ')
    ylim([min(dataw) max(dataw)])
    pause(0.1)
  end   
    
    


    
    ds =  (2*W(:,nzf+1:nz).^2./(newv2.^3));     
    J = ds.*(Gsf.*Grf);                                                                 % Jacobian
    gradient = real((J'*resf));                                                         % Gradient
    
    if optim == 0
       dm = -gradient;
    elseif optim == 1 && full_band == 0              
        Ha(:,:,iter) = J'*J;                                                            % Ha = Hessian approximate
        dm = -real(((Ha(:,:,iter) + 1e-2*(max(max(diag(Ha(:,:,iter)))))...              % Applying inverse of Ha to gradient
            *eye(size(Ha(:,:,iter))))\gradient));   
        figure(4)
        imagesc(real(Ha(:,:,iter)));
        title( ' Hessian Approximate' );
        xlabel('Depth(m)')
        ylabel('Depth(m)')
        pause(1)
    end
    
    
    %% step length
    alpha_t = 5e30;
    while max(max(abs(alpha_t*(dm)))) > max(max(b))/100
        alpha_t = 0.9 * alpha_t;
    end
    
    vv = newv + alpha_t*(dm);
    
    [Tvv,Rvvup,~,~]=layercac(vv,ro_b,dz,nt,dt,np,dp,p0,norm,mul,nprim);

    GBbp 	=	real(ifft(((ones(nt,1)-fft(RAdown).*fft(Rvvup)).^(-1)).*fft(TA)));

    calc_new	=	RAup + real(ifft(fft(TA).*fft(Rvvup).*fft(GBbp)));

    if full_band == 1
    r_calc = calc_new - calc;
    elseif full_band == 0
    calc_neww = real(ifft(fft(calc_new).*fft(wav)));
    r_calc = calc_neww - calcw;
    end
    alpha = -((r_calc'*res)/(r_calc'*r_calc))*alpha_t;
    
    
    
    %% Model update
    newv = newv + alpha*(dm);
    ModelMatrix(:,iter) = newv;
    updatedirection(:,iter) = alpha*dm;
    
    %% Making the Green's functions and predicted data for the next iteration
   
    
    [Gbb]=vsplayrc(newv,ro_b,dz,nt,dt,p0,norm);
    [~,Rbnewup,~,~]=layercac(newv,ro_b,dz,nt,dt,np,dp,p0,norm,mul,0);

    GBbp 	=	real(ifft(((ones(nt,1)-fft(RAdown).*fft(Rbnewup)).^(-1)).*fft(TA)));
    % New predicted data
    calc	=	RAup + real(ifft(fft(TA).*fft(Rbnewup).*fft(GBbp)));
    calcw = real(ifft(fft(calc).*fft(wav)));
    refv = [a;newv];
    
    % New Green's functions inside the target zone
    % with source at the surface and ricevers inside the target zone
    Gsnew = real(ifft(fft(Gbb).*repmat(fft(GBbp),1,length(b))));
    Gswnew = real(ifft(repmat(fft(wav),1,nz-nzf).*fft(Gsnew)));

    
    
    figure(2);
    subplot(2,1,1)
    plot(zvec,-gradient);
    title('Gradient');
    xlabel(' Depth(m) ')
    pause(1)
    
    subplot(2,1,2)
    plot(zvec,alpha*dm);
    title('update direction (alpha*Gradient)');
    xlabel(' Depth(m) ')
    pause(1)
    

    
    figure(3)
    subplot(1,2,1)
    plot(zvecr,realv)
    hold on
    plot(zvec,newv)
    title('updated model')
    ylabel( ' velocity (m/s) ' )
    xlabel( ' depth(m) ')
    hold off
    pause(1)
    subplot(1,2,2)
    plot(cost_norm)
    title('cost function')
    xlabel('iteration')
    ylabel('L_2 norm of residuals')
    pause(1)
    
    
    
    
    %     figure(5)
    %     plot(zvecr,realv);
    %     hold on
    %     plot(zvecr,refv(1:end))
    %     ylabel( ' velocity (m/s) ' )
    %     xlabel(' Depth(m) ')
    %
    
    
end
    


