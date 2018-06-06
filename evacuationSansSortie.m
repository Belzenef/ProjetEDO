% MODELE INDIVIDU-CENTRE
% Evacuation d'une salle avec des individus calmes, paniqués et appeurés
% Calmes (X)
% Apeurés (Y)
% Paniqués (Z) 
clear all; close all;

%% Paramètres des équations :

% Déplacement des individus
D = 0.2;    % coeficient de diffusion
Dmax = 3;   % densité max acceptable pour individus immobiles dans une foule (RATP)
vp=1.5;     % vitesse d'un individu paniqué (m/s)
vc=1;       % vitesse d'un individu calme (m/s)

% Diffusion de la panique
a1 =0.2;       % proportion d'individus calmes qui deviennent apeurés
a2 =0.01;      % proportion d'individus apeurés qui se calment
b1 =0.1;       % taux de contamination des gens apeurés par les gens paniqués
b2 =0.2;       % proportion de personnes calmes qui deviennent apeurés
c1 =0.2;       % proportion de personnes paniqués qui deviennent calmes

% Distribution des individus :
N = 160;       % population totale
Z0 = 10;        % nombre initial d'individus paniqués
Y0 = 20;        % nombre initial d'individus apeurés
X0 = N-Z0-Y0;   % les autres individus sont calmes
etat(N,3)=false;
etat(1:X0,1)=true;
etat((X0+1):(X0+Y0),2)=true;
etat((X0+Y0+1):(X0+Y0+Z0),3)=true;

%% Discretisation temporelle :
t0 = 0;          % temps initial       
tf = 100;         % temps final
dt = 0.1;       % pas de temps
X = zeros(tf/dt,1);  
Y = zeros(tf/dt,1);   
Z = zeros(tf/dt,1);

%% Discretisation spatiale
h=1;        % pas spatial
l=250;      % longueur salle (m)
L=80;       % Largeur salle (m)
xPorte=l/2; % Position de la porte
yPorte=0;
lPorte=4;   % longueur de la porte (m)
position = [l*rand([N,1]),L*rand([N,1])]; % positionnement des individus dans la salle
Densite = N/(l*L);

%% Initialisation
i=1;
ti=dt;  % temps courant
Xi=X0;
Yi=Y0;
Zi=Z0;
X(1,1)=X0;
Y(1,1)=Y0;
Z(1,1)=Z0;

%% Affichage Initial
figure(2);
%individus
x=position(:,1);
y=position(:,2);
calme=etat(:,1);
peur=etat(:,2);
panique=etat(:,3);
plot(x(calme),y(calme),'b.',x(peur),y(peur),'g.',x(panique),y(panique),'r.');
% murs de la salle
line([0 l],[L L],'Color','black');
line([0 l],[0 0],'Color','black');
line([l l],[0 L],'Color','black');
line([0 0],[0 L],'Color','black');
% porte
line([xPorte-lPorte/2 xPorte+lPorte/2],[yPorte yPorte],'Color',[0.5,0,0.5])
line([xPorte-lPorte/2 xPorte+lPorte/2],[yPorte-0.5 yPorte-0.5],'Color',[0.5,0,0.5])
axis([-5,l+5,-5,L+5])
    
%% Boucle principale
tic
while ti<tf
    
    n=length(x);
    drawnow;
    
    % Diffusion de la panique ---------------------------------------------
    
    % Schéma explicite pour les différents groupes
    X(i+1)= X(i) + dt*(-a1*X(i) + a2*Y(i) + c1*Z(i));
    Y(i+1)= Y(i) + dt*(a1*X(i) - a2*Y(i) + b2*Z(i) -b1*Y(i)*Z(i)*Densite);
    Z(i+1)= Z(i) + dt*(b1*Y(i)*Z(i)*Densite - b2*Z(i) - c1*Z(i));% Diffusion du carractère de panique
    
    % Changement d'état des individus
    Yi=round(Y(i));
    Xi=round(X(i));
    Zi=round(Z(i));
    etat(:,:)=false;
    etat(1:Xi,1)=true;
    etat((Xi+1):(Xi+Yi),2)=true;
    etat((Xi+Yi+1):N,3)=true;
    
    % Déplacement des individus -------------------------------------------
    
    % MAJ densité (autour de la porte)
    Densite=(N)/((max(y)-min(y))*(max(x)-min(x)))
    calme=etat(:,1);
    peur=etat(:,2);
    panique=etat(:,3);
    
    % Individus paniqués
    position(:,1) = x + (xPorte-x)*dt.*panique*(dt*(vp-Densite*vp/Dmax)) + randn(size(x))*sqrt(2*D*dt);
    position(:,2) = y + (yPorte-y)*dt.*panique*(dt*(vp-Densite*vp/Dmax)) + randn(size(y))*sqrt(2*D*dt);
    x=abs(position(:,1));
    y=abs(position(:,2));
    % Individus appeurés
    position(:,1) = x + (xPorte-x)*dt.*peur*(dt*(vp-Densite*vp/Dmax)) + randn(size(x))*sqrt(D*dt);
    position(:,2) = y + (yPorte-y)*dt.*peur*(dt*(vp-Densite*vp/Dmax)) + randn(size(y))*sqrt(D*dt);
    x=abs(position(:,1));
    y=abs(position(:,2));
    % Individus calmes
    position(:,1) = x + (xPorte-x)*dt.*calme*(dt*(vc-Densite*vc/Dmax));
    position(:,2) = y + (yPorte-y)*dt.*calme*(dt*(vc-Densite*vc/Dmax));
    x=abs(position(:,1));
    y=abs(position(:,2));
    
    % Affichage -----------------------------------------------------------
    figure(4);
    %individus
    plot(x(calme),y(calme),'b.',x(peur),y(peur),'g.',x(panique),y(panique),'r.');
    % murs de la salle
    line([0 l],[L L],'Color','black');
    line([0 l],[0 0],'Color','black');
    line([l l],[0 L],'Color','black');
    line([0 0],[0 L],'Color','black');
    % porte
    line([xPorte-lPorte/2 xPorte+lPorte/2],[yPorte yPorte],'Color',[0.5,0,0.5])
    line([xPorte-lPorte/2 xPorte+lPorte/2],[yPorte-0.5 yPorte-0.5],'Color',[0.5,0,0.5])
    axis([-5,l+5,-5,L+5]);
    
    
    ti=ti+dt;
    i=i+1;
end
toc

%% Affichage portrait de phase
temps=transpose((t0):dt:tf);

figure(1)
plot(temps,X,'blue');
hold on;
plot(temps,Y,'green');
hold on;
plot(temps,Z,'red');
hl = legend(['Calmes ';'Peur   ';'Panique']);
title('Evolution temporelle de la panique dans une population');    %% titre du graphe
