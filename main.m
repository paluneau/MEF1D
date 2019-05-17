close all;clear all;
%% Probl�me analytique (solution approch�e)
% Trouver u(x) t.q. u(x)-d^2_x[u(x)]=x avec conditions aux
% limites de Dirichlet homog�nes u(0)=u(1)=0

p = @(x)1;
q = @(x)1;
R = @(x)x;
sol_analytique=@(x)x-(sinh(x)./sinh(1));
[~,u_graph,err]=solve_mef1d(p,q,R,[0 1],20,[0 0],2,sol_analytique);
plot_mef1d(u_graph,sol_analytique,[0 1],err);

%% Probl�me de laplacien (avec solution exacte en P2)
% Trouver u(x) t.q. -d^2_x[u(x)]=5 avec conditions aux
% limites de Dirichlet u(0)=0 et u(5)=2
p1=@(x)0;
q1=@(x)1;
R1=@(x)5; % En mettant une constante non-nulle, on s'assure d'avoir un polynome de degr� 2 comme solution
Omega1 = [0 5];
dirichlet1 = [0 2];

sol_analytique_laplacien = @(x)-5/2.*x.^2+129/10.*x;

[~, uk1, err_L2] = solve_mef1d(p1,q1,R1,Omega1,10,dirichlet1,2,sol_analytique_laplacien);
plot_mef1d(uk1,sol_analytique_laplacien,Omega1,err_L2);


%% Probl�me de tension dans un c�ble
% Trouver u(x) t.q. -d_x[400*d_x[u(x)]]=r(x) avec conditions aux
% limites de Dirichlet homog�nes u(0)=u(5)=0
p2=@(x)0;
q2=@(x)400; % Tension dans le cable
R2=@(x)r(x); % Poids lin�aire du cable + Poids ajout�, d�fini plus bas.
dirichlet2=[0 0];

[~,uk2,~] = solve_mef1d(p2,q2,R2,Omega1,10,dirichlet2,2);
plot_mef1d(uk2);

% Repr�sente le poids lin�aire du c�ble et une masse non constante r�partie
% sur la partie gauche du c�ble.
function y = r(x)
[~,dim] = size(x);
y = zeros(1,dim);
for i=1:dim
    if 0<=x(i) && x(i)<=2
        y(i)=-6-40*x(i);
    elseif 2<=x(i) && x(i)<=5
        y(i)=-6;
    end
end
end
