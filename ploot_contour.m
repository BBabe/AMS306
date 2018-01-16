function ploot_contour(t,V,titre)
global x y 

Z=Vex(t,x,y);
contour(x,y,V',[0 0],'r'); %- contour
hold on
contour(x,y,Z',[0 0],'b'); %- contour
legend('front numerique','front exact');

grid on

%view(-15,25);
xlabel('x');
ylabel('y');
title(strcat(titre,'; t=',num2str(t)));



