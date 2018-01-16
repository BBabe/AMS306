function z=V0(x,y)
ox=ones(size(x));
oy=ones(size(y));
r=sqrt(x.^2*oy' + ox*(y.^2)');
z=min(1,r-1);



