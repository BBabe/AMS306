function ploot(t,V,titre)
global x y


X=x; 
Y=y; 
Z=V;

hold on
colormap cool
%surf(X,Y,Z,'EdgeColor',[.8 .8 .8]);  %- surface

%- Autres modes:
%surf(X,Y,Z,'EdgeColor',[.8 .8 .8],'FaceColor','none');
%surf(X,Y,Z);
%surfl(X,Y,Z);
surfc(X,Y,Z');  

contour3(X,Y,Z',[0 0],'r'); %- contour
grid on

view(-15,25);
xlabel('x');
ylabel('y');
title(strcat(titre,'; t=',num2str(t)));


