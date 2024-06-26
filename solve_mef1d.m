
function [u_nodal,u_graph,err_L2,err_H10,err_H1] = solve_mef1d(p,q,R,lim_dom,n_elem,cond_dirichlet,cond_neumann,deg,sol_analytique,sol_analytique_grad)
% SOLVE_MEF1D
% Calcule une solution approch�e de u(x) t.q.
% p(x)*u(x)-d_x[q(x)*d_x[u(x)]]=R(x) avec conditions aux limites de type Dirichlet ou Neumann.
%
% !!! Param�tres !!!
%
% - p,q,R: les param�tres du probl�me (function handles).
%
% - lim_dom: un vecteur [a b] qui repr�sente le domaine sur lequel on veut
% approximer u, l'intervalle [a,b].
%
% - n_elem: le nombre d'�l�ments voulu pour le maillage uniforme sur le domaine.
%
% - cond_dirichlet: un vecteur de la forme [x y], tel que u(a)=x et u(b)=y, les
% conditions aux limites du domaine. Si on ne veut pas imposer un bord, on
% doit mettre NaN.
%
% - cond_neumann: un vecteur de la forme [x y], tel que d_x[u](a)=x et d_x[u](b)=y, les
% conditions aux limites du domaine. Si on ne veut pas imposer un bord, on
% doit mettre NaN.
%
% - deg: le degr� des fonctions de base (1 ou 2)
%
% - sol_analytique (facultatif): Solution analytique du probl�me, pour
% calculer l'erreur sur la solution num�rique.
%
% !!! Sortie !!!
%
% - u_nodal: le vecteur des valeurs nodales de de la solution approch�e.
%
% - u_graph: un tableau 2 par n_elem qui contient la fonction u_i et l'�l�ment K_i
% sur lequel elle est d�finie.
%
% - err_L2: l'erreur en norme L^2 par rapport � la solution analytique
% fournie (si pas fournie, NaN).

%% Construction des tableaux

% Calcul du nombre de noeuds n�cessaires
n_total_noeuds_elem = 3;
if deg==1
    n_total_noeuds_elem = 2;
end
n_noeuds = n_total_noeuds_elem*n_elem-n_elem+1;

% Pr�traitement des conditions initiales dirichlet et des noeuds impos�s
noeuds_imposes = [];
cond_lim_valides = [];
if ~isnan(cond_dirichlet(1))
    noeuds_imposes = [noeuds_imposes 1];
    cond_lim_valides = [cond_lim_valides cond_dirichlet(1)];
end
if ~isnan(cond_dirichlet(2))
    noeuds_imposes = [noeuds_imposes n_noeuds];
    cond_lim_valides = [cond_lim_valides cond_dirichlet(2)];
end
[~,n_noeuds_imposes] = size(noeuds_imposes);


% Coordonn�es des noeuds 
coord = linspace(lim_dom(1),lim_dom(2),n_noeuds); % Maillage uniforme sur le domaine

% Connectivit� des noeuds
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

% Num�rotation des DDLs
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

% Traitement des conditions de Neumann
S = zeros(1,n_noeuds);
if ~isnan(cond_neumann(1))
    S(numer(connec(1,1)))=cond_neumann(1);
end
if ~isnan(cond_neumann(2))
    S(numer(connec(n_elem,2)))=cond_neumann(2);
end

% Rel�vement des conditions
U_g=zeros(n_noeuds,1); 
U_g(numer(noeuds_imposes))= cond_lim_valides;

%% Passage � l'�l�ment de r�f�rence
T=@(x1,x2,y)(x1+x2+(abs(x2-x1).*y))./2;
T_i=@(x1,x2,x)(2.*x-(x1+x2))./abs(x2-x1);

%% D�finition des fonctions d'interpolation de Lagrange sur l'�l�ment de r�f�rence
if deg==2 % Quadratique
    f1_hat = @(y)y.*(y-1)./2;
    f2_hat = @(y)y.*(y+1)./2;
    f3_hat= @(y)1-y.^2;
    fi_hat = {f1_hat f2_hat f3_hat};

    df1_hat = @(y)(2.*y-1)./2;
    df2_hat = @(y)(2.*y+1)./2;
    df3_hat = @(y)-2.*y;
    dfi_hat = {df1_hat df2_hat df3_hat};
elseif deg==1 % Lin�aire
    f1_hat = @(y)(1-y)./2;
    f2_hat = @(y)(y+1)./2;
    fi_hat = {f1_hat f2_hat};

    df1_hat = @(y)-1/2;
    df2_hat = @(y)1/2;
    dfi_hat = {df1_hat df2_hat};
end
    

%% Points d'int�gration (Gauss-Legendre 3 points)
ptsGauss = [-sqrt(3./5.) 0 sqrt(3./5.)];
wGauss = [5/9 8/9 5/9];

%% Assemblage (approche en correction)
A = sparse(n_noeuds,n_noeuds); % matrice de rigidit�
M = sparse(n_noeuds,n_noeuds); % matrice de masse (pour calculer l'erreur)
F = zeros(1,n_noeuds);
for i=1:n_elem % Boucle sur les �l�ments (i)
    x1 = coord(connec(i,1));
    x2 = coord(connec(i,2));
    h=abs(x2-x1);
    % Construction de la fonction de rel�vement locale
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
%% R�solution du syst�me global (approche en correction)
n_inconnus = n_noeuds - n_noeuds_imposes;
del_U_I = A(1:n_inconnus,1:n_inconnus)\((F(1:n_inconnus)+S(1:n_inconnus))');
del_U_C = zeros(n_noeuds_imposes,1);
del_U = [del_U_I;del_U_C];
U = del_U + U_g;
u_nodal = U(numer); % Renum�roter les DDLs selon l'ordre des noeuds
%% Solution par �l�ment + Calcul d'erreur
% Sur chaque �l�ment, on regarde les coefficients qu'on a trouv� et on
% exprime les points par combinaison lin�aire des fonctions de base.
if nargin >8
    err_L2 = 0;
    err_H10 = 0;
    err_H1 = 0;
else
    err_L2 = NaN;
    err_H10 = NaN;
    err_H1 = NaN;
end

u_graph = cell(3,n_elem);
for i = 1:n_elem
    x1 = coord(connec(i,1));
    x2 = coord(connec(i,2));
    h = abs(x2 - x1);
    K = linspace(x1,x2,(x2-x1)*100);
    u_i = u_nodal(connec(i,:));
    u = @(x)0;
    du = @(x)0;
    for j=1:n_total_noeuds_elem
        u = @(x)u(x)+u_i(j)*fi_hat{j}(x);
        du = @(x)du(x)+u_i(j)*dfi_hat{j}(x)*2/h;
    end
    u_graph{1,i}=@(x)u(T_i(x1,x2,x));
    u_graph{2,i}=@(x)du(T_i(x1,x2,x));
    if nargin>8
        err_L2 = err_L2 + integral(@(x)(u_graph{1,i}(x)-sol_analytique(x)).^2,x1,x2,'AbsTol',1e-15);
        err_H10 = err_H10 + integral(@(x)(u_graph{2,i}(x)-sol_analytique_grad(x)).^2,x1,x2,'AbsTol',1e-15);
    end
    u_graph{3,i}=K;
end
err_L2 = sqrt(err_L2);
err_H10 = sqrt(err_H10);
err_H1 = sqrt(err_L2^2 + err_H10^2);
end



