% MODELE INDIVIDU-CENTRE
% Evacuation d'une salle avec des individus calmes, paniqués et appeurés
clear all; close all;

%% Paramètres des équations :
D = 0.2;    % coeficient de diffusion
cp= 0.5;
vp=1.5;     % vitesse d'un individu paniqué (m/s)
vc=1;       % vitesse d'un individu calme (m/s)

%% Paramètre spatiaux
l=210;      % longueur salle (m)
L=80;       % Largeur salle (m)
xPorte=l/2; % Position de la porte
yPorte=0;
lPorte=4;   % longueur de la porte (m)

%% Discretisation spatial
h=1;            % pas spatial

%% Discretisation temporelle
t0=0;      % temps initial (s)
tf=100;     % temps final (s)
t=t0;      % temps courant (s)
dt=0.1;    % pas de temps (s)

%% Distribution des individus
% p=0 si calme ; p=1 si paniqué
N=170;         % population totale
P0=20;         % nb d'individus paniqués
p=false(1,N);   % distribution du carractère panique dans la population
p(1:P0)=true;
x=l*rand([N,1]);% répartition spatiale des individus
y=L*rand([N,1]);
w=zeros(1,N);   % probabilité de paniquer
Densite = N/(l*L);

%% Boucle Principale
% Affichage
figure(1);
plot(x(p),y(p),'r.',x(~p),y(~p),'b.');
line([xPorte-lPorte/2 xPorte+lPorte/2],[yPorte yPorte],'Color','green')
axis([0,l,0,L]);
    

tic
while t<tf
    drawnow;
    N=length(x);
    % Diffusion du carractère de panique
    % influencée par densité avoisinante
    w=dt*cp*(~p)*Densite;
    irep=find(rand(1,N)<w);
    for i=1:length(irep)
        [~,rempl] = min((x - x(irep(i))).^2 ...
           + (y - y(irep(i))).^2 ... 
           + 10*(x == x(irep(i))).*(y == y(irep(i))));
        p(rempl) = true;
    end
    
    % Déplacement des individus
    % vitesse diminuée en fonction de la densité
    T=transpose(p);
    x= x + (xPorte-x)*dt.*(T*(dt*(vp-Densite*vp/3)-dt*(vc-Densite*vc/3))+dt*(vc-Densite*vp/3)) ...
        + randn(size(x))*sqrt(2*D*dt);
   
    y= y + (yPorte-y)*dt.*(T*(dt*(vp-Densite*vp/3)-dt*(vc-Densite*vc/3))+dt*(vc-Densite*vp/3)) ...
        + randn(size(y))*sqrt(2*D*dt);
    
    x=abs(x);
    y=abs(y);
    
    % Sortie des individus
    j=1;
    while j<length(x)
        if (y(j)<=(yPorte+0.5) && x(j)<=(xPorte+lPorte/2) && x(j)>=(xPorte-lPorte/2))
            x(j)=[];
            p(j)=[];
            y(j)=[];
        end
        j=j+1;
    end
    
    % MAJ densité
    Densite=length(p)/((max(y)-min(y))*(max(x)-min(x)));
    
    % Affichage
    figure(1);
    plot(x(p),y(p),'r.',x(~p),y(~p),'b.');
    line([xPorte-lPorte/2 xPorte+lPorte/2],[yPorte yPorte],'Color','green')
    axis([0,l,0,L]);
    
    t=t+dt;
end
toc
