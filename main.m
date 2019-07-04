close all;clear all;

%% Solution manufactur�e (MMS): probl�me d'�lasticit� lin�aire (solution exacte P2)
% Champ de d�placement solution: u(x)=3*x^2.
% Mat�riau: Acier (Loi de comportement: Hooke, lambda=96.95 et mu=76.17)
% Probl�me: -div(sigma)=R(x) <=> -d_x[(lambda+2*mu)*d_x[u(x)]]=R(x)
% Terme source calcul�: R(x)=-1495.74
% Conditions aux limites: u(0)=0 et u(1)=3
lambda = 96.95;
mu = 76.17;
p = @(x)0;
q = @(x)(lambda+2*mu);
R = @(x)-1495.74;
sol_analytique = @(x)3*x.^2;
[~,u_graph,err] = solve_mef1d(p,q,R,[0 1],5,[0 3],[NaN NaN],2,sol_analytique);
plot_mef1d(u_graph,sol_analytique,[0 1],err);

%% Probl�me de dirichlet (solution approch�e)
% Trouver u(x) t.q. u(x)-d^2_x[u(x)]=x avec conditions aux
% limites de Dirichlet homog�nes u(0)=u(1)=0
p = @(x)1;
q = @(x)1;
R = @(x)x;
sol_analytique = @(x)x-(sinh(x)./sinh(1));
[~,u_graph,err] = solve_mef1d(p,q,R,[0 1],10,[0 0],[NaN NaN],1,sol_analytique);
plot_mef1d(u_graph,sol_analytique,[0 1],err);

%% Probl�me de laplacien dirichlet (avec solution exacte en P2)
% Trouver u(x) t.q. -d^2_x[u(x)]=5 avec conditions aux
% limites de Dirichlet u(0)=0 et u(5)=2
p = @(x)0;
q = @(x)1;
R = @(x)5; % En mettant une constante non-nulle, on s'assure d'avoir un polynome de degr� 2 comme solution
sol_analytique = @(x)-5/2.*x.^2+129/10.*x;
[~, uk1, err_L2] = solve_mef1d(p,q,R,[0 5],10,[0 2],[NaN NaN],2,sol_analytique);
plot_mef1d(uk1,sol_analytique,[0 5],err_L2);

%% Probl�me de laplacien dirichlet-neumann (solution approch�e)
% Trouver u(x) t.q. -d^2_x[u(x)]=x avec conditions aux
% limites u(0)=0 et d_x[u](1)=0
p = @(x)0;
q = @(x)1;
R = @(x)x;
sol_analytique = @(x)-(x.^3)./6+(1/2).*x;
[~,u_graph,err] = solve_mef1d(p,q,R,[0 1],50,[0 NaN],[NaN 0],1,sol_analytique);
plot_mef1d(u_graph,sol_analytique,[0 1],err);

%% Probl�me de tension dans un c�ble
% Trouver u(x) t.q. -d_x[400*d_x[u(x)]]=r(x) avec conditions aux
% limites de Dirichlet homog�nes u(0)=u(5)=0
p = @(x)0;
q = @(x)400; % Tension dans le cable
R = @(x)r(x); % Poids lin�aire du cable + Poids ajout�, d�fini plus bas 

[~,uk2,~] = solve_mef1d(p,q,R,[0 5],10,[0 0],[NaN NaN],2);
plot_mef1d(uk2);

% Repr�sente le poids lin�aire du c�ble et une masse non constante r�partie
% sur la partie gauche du c�ble.
function y = r(x)
[~,dim] = size(x);
y = zeros(1,dim);
for i = 1:dim
    if 0 <= x(i) && x(i)<=2
        y(i) = -6-40*x(i);
    elseif 2 <= x(i) && x(i)<=5
        y(i) = -6;
    end
end
end
