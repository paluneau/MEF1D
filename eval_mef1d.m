function [uhx, varargout] = eval_mef1d(U,x)
% EVAL_MEF1D
% Évalue la solution EF (et son gradient) en x.
%
% !!! Paramètres !!!
%
% - U: la sortie de la fonction SOLVE_MEF1D
%
% - x : le point à évaluer
%

[~,m] = size(U);
for i=1:m
    if x<=U{3,i}(end) && x>=U{3,i}(1)
        uhx = U{1,i}(x);
        if nargout>1
            varargout = {U{2,i}(x)};
        end
        return 
    end
end

uhx = NaN;
if nargout >1 
    varargout = {NaN};
end

end