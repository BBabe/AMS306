function z=Vex(t,x,y)
ox=ones(size(x));
oy=ones(size(y));
r1=sqrt((x-1).^2*oy' + ox*(y.^2)'); 
r2=sqrt((x+1).^2*oy' + ox*(y.^2)'); 
ro=0.5;
z=min(1,min(max(r1-ro-t,-1),max(r2-ro-t,-1)));

