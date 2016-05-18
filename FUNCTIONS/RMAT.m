% CODE FOR SIMULATION OF ORIGAMI SHEETS WITH SMOOTH ELASTIC FOLDS
% AUTHOR: Edwin A. Peraza Hernandez
% 05-16-2016
% PAPER DOI: http://dx.doi.org/10.1016/j.cad.2016.05.010

function R = RMAT(alpv, tetv)
% ROTATION CONSTRAINT MATRIX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUTS
%
% alpv: Array with face angles for the considered fold intersection
% tetv: Array with fold angles for the considered fold intersection
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% OUTPUTS
% 
% R: Rotation constraint matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R3 = @(t) [cos(t),  -sin(t),    0;
        sin(t),     cos(t),     0;
        0,      0,      1]; % MATRIX FOR ROTATION WITH RESPECT TO 3-AXIS
    
R1 = @(t) [1,   0,      0;
           0,   cos(t), -sin(t);
           0,   sin(t), cos(t)]; % MATRIX FOR ROTATION WITH RESPECT TO 1-AXIS

R = eye(3); % INITIALIZING MATRIX

for i = 1:size(alpv,2)
    R = R * R1(tetv(1,i))* R3(alpv(1,i));
end


