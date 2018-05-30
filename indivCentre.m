% MODELE INDIVIDU-CENTRE
% Dans une salle, individus Calmes vs individus Paniqués
clear all; close all;
%% Paramètres des équations :
D = 0.2;    % coeficient de diffusion
cp = 0.01; % coeficient de diffusion de la panique
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
Nl=round(l/h);
NL=round(L/h);

%% Discretisation temporelle
t0=0;      % temps initial (s)
tf=50;     % temps final (s)
t=t0;      % temps courant (s)
dt=0.1;    % pas de temps (s)

%% Distribution des individus
% p=0 si calme ; p=1 si paniqué
N=500;         % population totale
P0=100;         % nb d'individus paniqués
p=false(1,N);   % distribution du carractère panique dans la population
p(1:P0)=true;
x=l*rand([N,1]);% répartition spatiale des individus
y=L*rand([N,1]);
w=zeros(1,N);   % probabilité de paniquer
Densite = N/(l*L);

%% Boucle Principale

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
    x=x + (xPorte-x)*dt.*(transpose(p)*(dt*(vp-Densite*vp/3)-dt*(vc-Densite*vc/3))+dt*(vc-Densite*vp/3)) ...
        + randn(size(x))*sqrt(2*D*dt);
    x=abs(x);
    y=y + (yPorte-y)*dt.*(transpose(p)*(dt*(vp-Densite*vp/3)-dt*(vc-Densite*vc/3))+dt*(vc-Densite*vp/3)) ...
        + randn(size(y))*sqrt(2*D*dt);
    y=abs(y);
    
    % Sortie des individus
    j=1;
    while j<length(x)
        if (y(j)<=(yPorte+1) && x(j)<=(xPorte+lPorte/2) && x(j)>=(xPorte-lPorte/2))
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
    N
end
toc