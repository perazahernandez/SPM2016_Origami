% CODE FOR SIMULATION OF ORIGAMI SHEETS WITH SMOOTH ELASTIC FOLDS
% AUTHOR: Edwin A. Peraza Hernandez
% 05-16-2016
% PAPER DOI: http://dx.doi.org/10.1016/j.cad.2016.05.010

function Ui = Ui_FUN(I, EL, tet, FN)
% ELASTIC STRAIN ENERGY OF A FOLD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUTS
%
% tet: Array with all fold angles at the current increment (N_F x 1)
% FN: Fold index
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% OUTPUTS
%
% Ui: Elastic strain energy of the fold with index FN
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   


kapi = KAP(tet(FN,1), I.w(FN), I.ord);  % DETERMINING MAXIMUM CURVATURE

if I.ord == 2 % 2nd ORDER OF GEOMETRIC CONTINUITY
    Ui = (4/15)*( EL.E(FN,1) * I.Ii(FN) * I.w(FN) ) * kapi^2 ...
        / ( 1 - EL.nu(FN,1)^2);
elseif I.ord == 1 % 1st ORDER OF GEOMETRIC CONTINUITY
    Ui = (1/2)*( EL.E(FN,1) * I.Ii(FN) * I.w(FN) ) * kapi^2 ...
        / ( 1 - EL.nu(FN,1)^2);
end

