	
function[win]=tdeps(h,nt);	
win=ones(nt,1);
count=0;
for i=nt:-1:nt/2+2

if count==1
win(i)=0;
end

if abs(h(i)) > 0.005

if count==0
ii=i;
count=1;
end

end

end

if count==1
win(ii)=0;

end