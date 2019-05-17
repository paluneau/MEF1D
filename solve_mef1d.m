
function [u_nodal,u_graph,err_L2] = solve_mef1d(p,q,R,lim_dom,n_elem,cond_lim,deg,sol_analytique)
% SOLVE_MEF1D
% Calcule une solution approchée de u(x) t.q.
% p(x)*u(x)-d_x[q(x)*d_x[u(x)]]=R(x) avec conditions aux limites de type Dirichlet.
%
% !!! Paramètres !!!
%
% - p,q,R: les paramètres du problème (function handles).
%
% - lim_dom: un vecteur [a b] qui représente le domaine sur lequel on veut
% approximer u, l'intervalle [a,b].
%
% - n_elem: le nombre d'éléments voulu pour le maillage uniforme sur le domaine.
%
% - cond_lim: un vecteur de la forme [x y], tel que u(a)=x et u(b)=y, les
% conditions aux limites du domaine. Si on ne veut pas imposer un bord, on
% doit mettre NaN.
%
% - deg: le degré des fonctions de base (1 ou 2)
%
% - sol_analytique (facultatif): Solution analytique du problème, pour
% calculer l'erreur sur la solution numérique.
%
% !!! Sortie !!!
%
% - u_nodal: le vecteur des valeurs nodales de de la solution approchée.
%
% - u_graph: un tableau 2 par n_elem qui contient la fonction u_i et l'élément K_i
% sur lequel elle est définie.
%
% - err_L2: l'erreur en norme L^2 par rapport à la solution analytique
% fournie (si pas fournie, NaN).

%% Construction des tableaux

% Calcul du nombre de noeuds nécessaires
n_total_noeuds_elem = 3;
if deg==1
    n_total_noeuds_elem = 2;
end
n_noeuds = n_total_noeuds_elem*n_elem-n_elem+1;

% Prétraitement des conditions initiales et des noeuds imposés
noeuds_imposes = [];
cond_lim_valides = [];
if ~isnan(cond_lim(1))
    noeuds_imposes = [noeuds_imposes 1];
    cond_lim_valides = [cond_lim_valides cond_lim(1)];
end
if ~isnan(cond_lim(2))
    noeuds_imposes = [noeuds_imposes n_noeuds];
    cond_lim_valides = [cond_lim_valides cond_lim(2)];
end
[~,n_noeuds_imposes] = size(noeuds_imposes);

% Coordonnées des noeuds 
coord = linspace(lim_dom(1),lim_dom(2),n_noeuds); % Maillage uniforme sur le domaine

% Connectivité des noeuds
connec = zeros(n_elem,n_total_noeuds_elem);
for i = 1:n_elem
    k = 1;
    init = 1+(i-1)*(n_total_noeuds_elem-1);
    j = init;
    saut = n_total_noeuds_elem - 1; 
    while j<=init+(n_total_noeuds_elem-1)
       connec(i,k)=j;
       j=j+saut;
       k=k+1;
    end
    j = init + 1;
    while k<=n_total_noeuds_elem
       connec(i,k)=j;
       j=j+1;
       k=k+1;
    end
end

% Numérotation des DDLs
numer = zeros(1,n_noeuds); 
numer(noeuds_imposes(:))=-1;
impose = n_noeuds - n_noeuds_imposes+1;
j=1;
for i=1:n_noeuds
    if(numer(i)==-1)
        numer(i) = impose;
        impose = impose + 1;
    else
        numer(i) = j;
        j=j+1;
    end
end

% Relèvement des conditions
U_g=zeros(n_noeuds,1); 
U_g(numer(noeuds_imposes))= cond_lim_valides;

%% Passage à l'élément de référence
T=@(x1,x2,y)(x1+x2+(abs(x2-x1).*y))./2;
T_i=@(x1,x2,x)(2.*x-(x1+x2))./abs(x2-x1);

%% Définition des fonctions d'interpolation de Lagrange sur l'élément de référence
if deg==2 % Quadratique
    f1_hat = @(y)y.*(y-1)./2;
    f2_hat = @(y)y.*(y+1)./2;
    f3_hat= @(y)1-y.^2;
    fi_hat = {f1_hat f2_hat f3_hat};

    df1_hat = @(y)(2.*y-1)./2;
    df2_hat = @(y)(2.*y+1)./2;
    df3_hat = @(y)-2.*y;
    dfi_hat = {df1_hat df2_hat df3_hat};
elseif deg==1 % Linéaire
    f1_hat = @(y)(1-y)./2;
    f2_hat = @(y)(y+1)./2;
    fi_hat = {f1_hat f2_hat};

    df1_hat = @(y)-1/2;
    df2_hat = @(y)1/2;
    dfi_hat = {df1_hat df2_hat};
end
    

%% Points d'intégration (Gauss-Legendre 3 points)
ptsGauss = [-sqrt(3./5.) 0 sqrt(3./5.)];
wGauss = [5/9 8/9 5/9];

%% Assemblage
A = zeros(n_noeuds,n_noeuds); % matrice de rigidité
M = zeros(n_noeuds,n_noeuds); % matrice de masse (pour calculer l'erreur)
F = zeros(1,n_noeuds);
for i=1:n_elem % Boucle sur les éléments (i)
    x1 = coord(connec(i,1));
    x2 = coord(connec(i,2));
    h=abs(x2-x1);
    % Construction de la fonction de relèvement locale
        rel_loc=@(x)0;
        drel_loc=@(x)0;
        for c=1:n_total_noeuds_elem
           rel_loc=@(x)rel_loc(x)+U_g(numer(connec(i,c)))*fi_hat{c}(x);
           drel_loc=@(x)drel_loc(x)+U_g(numer(connec(i,c)))*dfi_hat{c}(x);
        end
    for j=1:n_total_noeuds_elem % Boucle sur les DDLs (j), 1 par noeud
        adr1=numer(connec(i,j)); 
        Fj=@(y)((h/2).*R(T(x1,x2,y)).*fi_hat{j}(y))-((h/2).*p(T(x1,x2,y)).*fi_hat{j}(y).*rel_loc(y))-((2/h).*q(T(x1,x2,y)).*dfi_hat{j}(y).*drel_loc(y));
        F(adr1)=F(adr1)+(Fj(ptsGauss) * wGauss');
        for k=1:n_total_noeuds_elem % Boucle sur les DDLs (k), 1 par noeud
            adr2 = numer(connec(i,k));
            Ajk = @(y)(h/2).*p(T(x1,x2,y)).*fi_hat{k}(y).*fi_hat{j}(y)+(2/h).*q(T(x1,x2,y)).*dfi_hat{k}(y).*dfi_hat{j}(y);
            Mjk = @(y)(h/2).*fi_hat{k}(y).*fi_hat{j}(y);
            A(adr1,adr2)=A(adr1,adr2)+(Ajk(ptsGauss) * wGauss');
            M(adr1,adr2)=M(adr1,adr2)+(Mjk(ptsGauss) * wGauss'); % Pour calculer l'erreur
        end  
    end  
end
%% Résolution du système global (approche en correction)
[~,dim22]=size(noeuds_imposes);
dim11 = n_noeuds - dim22;
del_U_I = A(1:dim11,1:dim11)\(F(1:dim11)');
del_U_C = zeros(dim22,1);
del_U = [del_U_I;del_U_C];
U = del_U + U_g;
u_nodal = U(numer); % Renuméroter les DDLs selon l'ordre des noeuds
%% Calcul de l'erreur en norme L^2 (si sol_analytique donnée)
if nargin >7
    ecart = u_nodal-sol_analytique(coord');
    err_L2 = sqrt(ecart'*M*ecart);
else
    err_L2 = NaN;
end
%% Solution par élément
% Sur chaque élément, on regarde les coefficients qu'on a trouvé et on
% exprime les points par combinaison linéaire des fonctions de base.
u_graph = cell(2,n_elem);
for i = 1:n_elem
    x1 = coord(connec(i,1));
    x2 = coord(connec(i,2));
    K = linspace(x1,x2,(x2-x1)*20);
    u_i = u_nodal(connec(i,:));
    u = @(x)0;
    for j=1:n_total_noeuds_elem
        u = @(x)u(x)+u_i(j)*fi_hat{j}(x);
    end
    u_graph{1,i}=@(x)u(T_i(x1,x2,x));
    u_graph{2,i}=K;
end
end

