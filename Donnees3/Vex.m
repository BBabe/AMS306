function z=Vex(t,x,y)
Rmt=[cos(-t) -sin(-t); sin(-t) cos(-t)];
z=zeros(length(x),length(y));
for i=1:length(x)
for j=1:length(y)
  X=[x(i); y(j)];
  Y=Rmt*X; 
  z(i,j)=V0(Y(1),Y(2));
end
end
