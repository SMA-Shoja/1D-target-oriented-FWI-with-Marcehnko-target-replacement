function[wav]=ricker(f0,nt,dt);	

t	=	[-nt/2*dt:dt:(nt/2-1)*dt]';
wavsym	=	(1-2*pi^2*f0^2*t.^2).*exp(-pi^2*f0^2*t.^2);
wav	=	[wavsym(nt/2+1:nt);wavsym(1:nt/2)];