function [uhx, varargout] = eval_mef1d(U,x,grad)
% EVAL_MEF1D
% Évalue la solution EF (et son gradient) en x.
%
% !!! Paramètres !!!
%
% - U: la sortie de la fonction SOLVE_MEF1D
%
% - x : le point à évaluer
%
% - grad (facultatif) : booleen qui indique si on veut sortir le gradient
% en premier
%

uhx = NaN;
if nargout >1 
    varargout = {NaN};
end

fir = 1;
sec = 2;

if nargin > 2
    if grad
        fir = 2;
        sec = 1;
    end
end

[~,m] = size(U);
for j=1:length(x)
    for i=1:m
        if x(j)<=U{3,i}(end) && x(j)>=U{3,i}(1)
            uhx(j) = U{fir,i}(x(j));
            if nargout>1
                varargout(j) = {U{sec,i}(x(j))};
            end
            continue
        end
    end
end


end