function [T,R,F1p,F1m]=layercac(cp,rho,dz,nt,dt,np,dp,p0,norm,mul,nprim)
% LAYERCODE 	Compute the  response of
%		an 1D acoustic medium
%
% syntax: [T,R]=layercode(cp,rho,dz,nt,dt,np,dp,p0,norm,mul)
%
%
% R	= Reflection response (t,p)
% T	= Transmission response (t,p)
%
% Cp	= Velocity log
% Rho	= Density log
% dz	= depth step
% nt	= number of time samples
% dt	= time sampling
% np	= number op rayparameters
% dp	= rayparameter sampling
% p0	= first rayparameter
% norm	= 0: flux normalization;  1: pressure normalization
% mul = 0: no multiples;  1: multiples
% nprim = 0: standard scheme, otherwise generate only the primary response of layer nprim


%number of layers
%------------------
if length(cp)~=length(rho)
    disp('WARNING: discripance between density and velocity log');
    disp(' Smallest is chosen!');
end
N	= min([length(cp) length(rho)]);


N	=	N+1;
cp(N)	=	cp(N-1);
rho(N)	=	rho(N-1);


%frequencies
%------------------
om      = (0:(nt/2)).*(2*pi/(nt*dt));
om      = om.';
nf      = (nt/2)+1;
p     	= (0:np-1).*dp + p0;


%initialise the GLOBAL quantities
%-------------------------
Rd    =    zeros(nf,np);
Ru    =    zeros(nf,np);
Td    =    ones(nf,np);
Tu    =    ones(nf,np);
T2    =    ones(nf,np);

%start the recursion loop over the N-1 layers
for n=1:N-1
    
    if 10*round(n/10) == n
        n;
    end
    
    %calculate the local operators
    %-vertical slowness- (number of layers (N), number of p(np))
    
    q1     =    sqrt(cp(n).^(-2)-p.^2);
    q2     =    sqrt(cp(n+1).^(-2)-p.^2);
    
    %-reflection coefficients- and
    %-flux normalised transmission coefficients-
    
    r     =     (rho(n+1).*q1-rho(n).*q2)./(rho(n+1).*q1+rho(n).*q2);
    if norm == 0
        td    =     sqrt(rho(n+1).*q1.*rho(n).*q2)./(0.5.*(rho(n+1).*q1+rho(n).*q2));
        tu    =     td;
    else
        td    =	    1 + r;
        tu    =     1 - r;
    end
    %if mul==0
    %	td=1;
    %	tu=1;
    %end
    t2    =     tu.*td;
    
    %calculate the phase shift operator
    
    q     =    ones(nf,1)*q1;
    q     =    real(q) + 1i*(sign(real(om))*ones(1,np)).*imag(q);
    r     =    ones(nf,1)*r;
    td    =    ones(nf,1)*td;
    tu    =    ones(nf,1)*tu;
    t2    =    ones(nf,1)*t2;
    om1   =    om*ones(1,np);
    w     =    exp(1i*om1.*q*dz);
    w2    =    w.^2;
    
    %calculate the multiple operator
    if mul == 1
        M     =    (1-w2.*r.*Ru).^(-1);
    else
        M	=	1;
    end
    
    if nprim==0
        
        %calculate the R downgoing
        Rd    =    Rd + T2.*w2.*r.*M;
        
    else
        if n==nprim
            
            Rd    =  Rd + T2.*w2.*r.*M;
        end
    end
    
    %calculate the R upgoing
    Ru    =    -r + t2.*w2.*Ru.*M;
    %calculate the T downgoing
    Td    =    td.*w.*Td.*M;
    %calculate the T upgoing
    Tu    =    tu.*w.*Tu.*M;
    %calculate the T square
    T2    =    Td.*Tu;
    
    
end

Fp	=	Td.^(-1);
Fm	=	Rd.*Fp;

%calculate the inverse fft's

T = [Td(1:nt/2,:) ;real(Td(nt/2+1,:)) ;conj(Td(nt/2:-1:2,:))];
T         =    real(ifft(conj(T)));


R = [Rd(1:nt/2,:) ;real(Rd(nt/2+1,:)) ;conj(Rd(nt/2:-1:2,:))];
R         =    real(ifft(conj(R)));

F1p = [Fp(1:nt/2,:) ;real(Fp(nt/2+1,:)) ;conj(Fp(nt/2:-1:2,:))];
F1p         =    real(ifft(conj(F1p)));

F1m = [Fm(1:nt/2,:) ;real(Fm(nt/2+1,:)) ;conj(Fm(nt/2:-1:2,:))];
F1m         =    real(ifft(conj(F1m)));





