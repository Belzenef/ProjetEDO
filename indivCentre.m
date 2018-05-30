% MODELE INDIVIDU-CENTRE
% Dans une salle, individus Calmes vs individus Paniqués

%% Paramètres des équations :
D = 0.5; % coeficient de diffusion
cp = 0.1; % coeficient de panique

%% Paramètre spatiaux
l=210;  % longueur salle
L=80;   % Largeur salle 

%% Discretisation spatial
h=1;    % pas spatial
Nl=round(l/h);
NL=round(L/h);
Xi=0:h:l;
Yi=0:h:L;


%% Discretisation temporelle
t0=0;
tf=50;
t=t0;
dt=0.1;

%% Distribution des individus
% p=0 si calme ; p=1 si paniqué
N=1700; % population totale
P0=100; % nb d'individus paniqués
p=false(1,N); % distribution du carractère panique dans la population
p(1:P0)=true;
x=l*rand([N,1]); % répartition spatiale des individus
y=L*rand([N,1]);
w=zeros(1,N) % probabilité de paniquer

%% Boucle Principale
figure(1); clf;
plot(x(p),y(p),'r.',x(~p),y(~p),'b.')

tic
while t<tf
    drawnow;
    
    % Diffusion du carractère de panique
    w=dt*cp*(1-p);
    irep=find(rand(1,N)<w);
    for i=1:length(irep)
        [~,rempl] = min((x - x(irep(i))).^2 ...
           + (y - y(irep(i))).^2 ... 
           + 10*(x == x(irep(i))).*(y == y(irep(i))));
        p(rempl) = true;
    end
    
    % Déplacement aléatoire des individus
    x=x+sqrt(dt)*sqrt(2*D)*randn(size(x));
    x=abs(x);
    y=y+sqrt(dt)*sqrt(2*D)*randn(size(y));
    y=abs(y);
    
    % Affichage
    plot(x(p),y(p),'r.',x(~p),y(~p),'b.')    
    axis([0,l,0,L]);
    t=t+dt;
end
toc
