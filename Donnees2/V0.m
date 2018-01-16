function z=V0(x,y)
ox=ones(size(x));
oy=ones(size(y));
r1=sqrt((x-1).^2*oy' + ox*(y.^2)');  %- r1=||X-A||
r2=sqrt((x+1).^2*oy' + ox*(y.^2)');  %- r2=||X-B||
ro=0.5;
z=min(1,min(r1-ro,r2-ro));


