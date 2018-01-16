%- tp2: Methode Semi-Lagrangienne
%- Janvier 2017
clear all;

global x y

TEST   = 1;		%- TEST=1,2, ou 3
SCHEMA = 'ENO2';	%- SCHEMA='SL' (4 points) - 'SL2' (3points) - 'SL-RK2' - 'LLF'
                %- 'ENO2'
fprintf('TEST=%i, SCHEMA=%s\n', TEST, SCHEMA)
if strcmp(SCHEMA,'LLF')
    cx = 1 ; cy = 1;
end

%!!---------------------------------------------------------------------!!
%- Donnees du probleme: vitesse, Condition initiale / solution exacte
%       (+dynamique modifiee eventuellement)
%!!---------------------------------------------------------------------!!

switch TEST
    case 1
        %--------------------------------------
        %- Exemple avec cdi v_0(x)=min(1,|x|-1).
        %--------------------------------------
        
        Nu=6; ulist=(2*pi) * (0:Nu-1)/Nu;
        A=2;  T=1.; Mx=21; My=Mx; N=10;
        
        Vbord=1;
        
        %- Vitesses (f1,f2), fonctions V0, Vex;
        addpath('Donnees1')
        
    case 2
        %--------------------------------------
        %- Exemple des deux trous en expansion: v_0=min(1,|x-A|-0.5,|x-B|-0.5)
        %--------------------------------------
        Nu=6; ulist=(2*pi) * (0:Nu-1)/Nu;
        A=2;  T=1; Mx=20; My=Mx; N=50;
        %A=2;  T=1; Mx=30; My=Mx; N=30;
        
        Vbord=1;
        
        %- Vitesses (f1,f2), fonctions V0, Vex;
        addpath('Donnees2')
        
    case 3
        %-------------------------------------------------
        %- Exemple de la rotation (+ expansion lente eventuelle)
        %-    Mx=My=20: on observe qqsoit la CFL que le front numerique "disparait"
        %-    au bout d'un demi-tour environ
        %-    [C'est un pb de diffusion numerique]
        %----------------------------------------------------
        Vbord=1;
        Nu=1; ulist=[0];
        A=2;  T=2*pi; Mx=30; My=Mx; N=10;
        
        %- Vitesses (f1,f2), fonctions V0, Vex;
        addpath('Donnees3')
        
end

%!!-------------------------------------------!!
%- Maillage. CFL. INITIALISATION DES GRAPHIQUES
%!!-------------------------------------------!!

dt=T/N;
Xmin=-A; Xmax=A; Ymin=-A; Ymax=A;

%- Maillage avec Mx * My points INCLUANT le bord du domaine
hx=(Xmax-Xmin)/(Mx-1);
hy=(Ymax-Ymin)/(My-1);
x=Xmin+(0:Mx-1)'*hx;  %- Maillage x spatial =  x_0 ... x_{Mx-1}
y=Ymin+(0:My-1)'*hy;  %- Maillage y spatial =  y_0 ... y_{My-1}


%- Initialisation
V=V0(x,y);
%- CFL
switch TEST
    case 1
        CFL=0; fprintf('CFL A COMPLETER!')
    case 2
        CFL=0; fprintf('CFL A COMPLETER!')
    case 3
        CFL=0; fprintf('CFL A COMPLETER!')
end
fprintf('CFL= %f\n',CFL);

%- Initialisations pour Graphique
Zmin=-1; Zmax=1;

%- Graphique initial
%- Affichage de la solution numerique - et de la courbe V=0
t=0;
clf;
subplot(121); ploot(t,V,'Sol. Numerique');
subplot(122); ploot_contour(t,V,'Frontiere V=0');

fprintf('appuyer sur la touche Return');
input('');

%!!-------------------------------------------------------!!
%- Boucle principale: Calcul de V^{n+1} en fonction de V^n
%!!-------------------------------------------------------!!
for n=0:N-1
    
    t=(n+1)*dt;
    
    switch SCHEMA
        
        case 'SL'
            %- Methode SL avec interpolation sur 4 points.
            %-Calcul d'un nouveau V en fonction du precedent.
            %-On notera Vu les calculs intermediaires pour chaque u.
            %-On notera Vnew le min des vecteurs Vu precedement calcule's.
            
            for iu=1:length(ulist);
                u=ulist(iu);
                for i=1:Mx
                    for j=1:My
                        %- On remonte la caracteristique sur un temps -dt  par Euler
                        xb=x(i)- f1(u,x(i),y(j))*dt;
                        yb=y(j)- f2(u,x(i),y(j))*dt;
                        
                        if (xb<=Xmin | xb>=Xmax | yb<=Ymin | yb>=Ymax);
                            Vu(i,j)=Vbord;
                        else
                            
                            %- Determination du carre contenant Xb=[xb,yb], puis interpolation
                            m=floor((xb-Xmin)/hx)+1;
                            l=floor((yb-Ymin)/hy)+1;
                            p=(xb-x(m))/hx;
                            q=(yb-y(l))/hy;
                            Vu(i,j) = (1-p)*(1-q)*V(m,l) ...
                                +p*(1-q)*V(m+1,l) +(1-p)*q*V(m,l+1) +p*q*V(m+1,l+1);
                        end
                    end
                end
                if iu==1;  Vnew=Vu; else; Vnew=min(Vnew,Vu); end;
            end
            V=Vnew;
            
        case 'SL2'
            %- Methode SL avec interpolation sur 3 points (une solution).
            
            fprintf('A vous de programmer cette methode !')
            
        case 'SL-RK2'
            %- Methode SL interpolation sur 4 points + Runge Kutta ordre 2
            %Calcul d'un nouveau V en fonction du precedent.
            %On notera Vu les calculs intermediaires pour chaque u.
            %On notera Vnew le min des vecteurs Vu precedement calcule's.
            
            fprintf('A vous de programmer cette methode!')
            
        case 'LLF'
            % Local Lax Friedrichs
            Vold = zeros(Mx+2,My+2);
            Vold(2:Mx+1,2:My+1) = V;
            Vold(1,:) = Vbord; Vold(Mx+2,:   ) = Vbord;
            Vold(:,1) = Vbord; Vold(:   ,My+2) = Vbord;
            uxp = (Vold(3:Mx+2,2:My+1) - Vold(2:Mx+1,2:My+1)) / hx;
            uxm = (Vold(2:Mx+1,2:My+1) - Vold(1:Mx  ,2:My+1)) / hx;
            uyp = (Vold(2:Mx+1,3:Mx+2) - Vold(2:Mx+1,2:Mx+1)) / hy;
            uym = (Vold(2:Mx+1,2:Mx+1) - Vold(2:Mx+1,1:Mx  )) / hy;
            g = sqrt(((uxp+uxm)/2).^2+((uyp+uym)/2).^2) -cx/2*(uxp-uxm) -cy/2*(uyp-uym);
            V = V -dt*g;
            
        case 'ENO2'
            Vold = zeros(Mx+4,My+4); C = 2;
            Vold(C+(1:Mx),C+(1:My)) = V;
            Vold(C+0,:  ) = Vbord; Vold(C+Mx+1,:     ) = Vbord;
            Vold(:  ,C+0) = Vbord; Vold(:     ,C+My+1) = Vbord;
            Vold(C-1,:  ) = Vbord; Vold(C+Mx+2,:     ) = Vbord;
            Vold(:  ,C-1) = Vbord; Vold(:     ,C+My+2) = Vbord;
            uxp = (Vold(C+(2:Mx+1),C+(1:My)) - Vold(C+(1:Mx  ),C+(1:My))) / hx;
            uxm = (Vold(C+(1:Mx  ),C+(1:My)) - Vold(C+(0:Mx-1),C+(1:My))) / hx;
            uyp = (Vold(2:Mx+1,3:Mx+2) - Vold(2:Mx+1,2:Mx+1)) / hy;
            uym = (Vold(2:Mx+1,2:Mx+1) - Vold(2:Mx+1,1:Mx  )) / hy;
            g = sqrt(((uxp+uxm)/2).^2+((uyp+uym)/2).^2) -cx/2*(uxp-uxm) -cy/2*(uyp-uym);
            V = V -dt*g;
            
        otherwise
            fprintf('SCHEMA not programmed\n'); abort
            
    end
    
    
    %- Affichages: solution numerique - courbes niveaux V=0
    clf
    subplot(121); ploot(t,V,'Sol. Numerique');
    subplot(122); ploot_contour(t,V,'Frontiere V=0');
    
    %- Affichage de l'erreur Linfty ou de l'erreur L1:
    %printf('t=%5.2f, ErrLinfty:%f \n',t,norm(V-Vex(t,x,y),'infty'));
    fprintf('t=%5.2f, ErrL1    :%f\n',t,norm(V-Vex(t,x,y),1)*hx*hy);
    
    %fprintf(' (Taper q pour quit)');
    %out=input('','s'); if out=='q'; abort; end;
    
    %input('','s')
    %fprintf('\n');
    
    pause(0.2)
end

%- tp2-3: avec Mx=My=20, N=10: cfl~1,  eL1=0.13
%-        avec Mx=My=20, N=1:  cfl~10, eL1=0.23



