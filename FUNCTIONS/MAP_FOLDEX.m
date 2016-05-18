% CODE FOR SIMULATION OF ORIGAMI SHEETS WITH SMOOTH ELASTIC FOLDS
% AUTHOR: Edwin A. Peraza Hernandez
% 05-16-2016
% PAPER DOI: http://dx.doi.org/10.1016/j.cad.2016.05.010

function MFE = MAP_FOLDEX(I, tet, FN)
% RETURNS 4 x 4 MATRIX WITH FURTHER TRANSFORMATIONS FOR FOLD MAPPING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% INPUTS
% 
% tet: Array with all fold angles at the current increment (N_F x 1)
% FN: Fold number
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% OUTPUT
%
% MFE: 4 x 4 MATRIX WITH FURTHER TRANSFORMATION FOR FOLD MAPPING
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R3 = @(t) [cos(t)  -sin(t)    0  0
        sin(t)     cos(t)     0  0
        0   0   1  0
        0   0   0  1]; % MATRIX FOR ROTATION WITH RESPECT TO 3-AXIS IN HOMOGENEOUS COORDINATES
    
R1 = @(t) [1   0   0   0
           0  cos(t)  -sin(t)  0
           0  sin(t)  cos(t)   0
           0   0   0   1]; % MATRIX FOR ROTATION WITH RESPECT TO 1-AXIS IN HOMOGENEOUS COORDINATES
       
B = @(b) [1  0   0  b(1,1)
    0   1   0   b(2,1)
    0   0   1   b(3,1)
    0   0   0   1]; % MATRIX FOR TRANSLATION BY VECTOR b IN HOMOGENEOUS COORDINATES

wstr = WSTR(tet(FN), I.w(FN), I.ord); % w PARAMETER

MFE = B( [I.p1(FN, 1:2)'; 0] ) ...
    *R3(I.halp(FN)) ...
    *R1(0.5 * tet(FN)) ...
    *B( [0; wstr/2; 0] ); 
    