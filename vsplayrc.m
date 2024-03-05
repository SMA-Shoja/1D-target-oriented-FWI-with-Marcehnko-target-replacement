function [VSP]=vsplayrc(cp,rho,dz,nt,dt,p,norm)
% LAYERCODE 	Compute the VSP response of 
%		an 1D acoustic medium for one rayparameter
% 
% syntax: [VSP]=vsplayrc(Cp,Rho,dz,nt,dt,p,norm)
%
%
% VSP	= VSP response (t,z)
%
% Cp	= Velocity log
% Rho	= Density log
% dz	= depth step
% nt	= number of time samples
% dt	= time sampling
% p	= rayparameter 
% norm	= 0: flux normalization;  1: pressure normalization
%

%number of layers
%------------------
if length(cp)~=length(rho)
	disp('WARNING: discripance between density and velocity log');
	disp(' Smallest is chosen!');
end;
N	= min([length(cp) length(rho)]);

%frequencies
%------------------
om      = (0:(nt/2)).*(2*pi/(nt*dt));
om      = om.';
nf      = (nt/2)+1; 
np	= 1;


%initialise the GLOBAL quantities
%-------------------------
Rd    =    zeros(nf,np);
Ru    =    zeros(nf,np);
Td    =    ones(nf,np);
Tu    =    ones(nf,np);
T2    =    ones(nf,np);

TT1   =    zeros(nf,N-1);
RR1   =    zeros(nf,N-1);

%start the first (downward) recursion loop over the N-1 layers
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
	td    =	    sqrt(rho(n+1).*q1.*rho(n).*q2)./(0.5.*(rho(n+1).*q1+rho(n).*q2));
	tu    =	    td;
else
	td    =	    1 + r;
	tu    =     1 - r;
end
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
	M     =    (1-w2.*r.*Ru).^(-1);
	%calculate the R downgoing
	Rd    =    Rd + T2.*w2.*r.*M;
	%calculate the R upgoing
	Ru    =    -r + t2.*w2.*Ru.*M;
	%calculate the T downgoing
	Td    =    td.*w.*Td.*M;
	%calculate the T upgoing
	Tu    =    tu.*w.*Tu.*M;
	%calculate the T square
	T2    =    Td.*Tu;

	TT1(:,n) =   Td;
	RR1(:,n) =   Ru;
end

%reset the GLOBAL quantities
%-------------------------

Rd    =    zeros(nf,np);

RR2   =    zeros(nf,N-1);

%start the second (upward) recursion loop over the N-1 layers
for m=1:N-2


	n      =    N-m;

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
	td    =	    sqrt(rho(n+1).*q1.*rho(n).*q2)./(0.5.*(rho(n+1).*q1+rho(n).*q2));
	tu    =     td;
else
	td    =	    1 + r;
	tu    =     1 - r;
end
	t2    =     tu.*td;

	%calculate the phase shift operator

	q     =    ones(nf,1)*q1;
	q     =    real(q) + 1i*(sign(real(om))*ones(1,np)).*imag(q);
	r     =    ones(nf,1)*r;
	t2    =    ones(nf,1)*t2;
	om1   =    om*ones(1,np);
	w     =    exp(1i*om1.*q*dz);
	w2    =    w.^2;

	%calculate the multiple operator
	M     =    (1+r.*Rd).^(-1);
	%calculate the R downgoing
	Rd    =    w2.*(r+t2.*Rd.*M);

	RR2(:,n-1) =   Rd;
end


 VSP   =    ((1-RR2.*RR1).^(-1)).* (TT1 + RR2.*TT1);


%calculate the inverse fft's 




VSP = [VSP(1:nt/2,:) ;real(VSP(nt/2+1,:)) ;conj(VSP(nt/2:-1:2,:))];
VSP         =    real(ifft(conj(VSP)));

[T,R]=layerc(cp,rho,dz,nt,dt,1,0,p,norm,1);
VSP=[R VSP];
VSP(1,1)=VSP(1,1)+1;