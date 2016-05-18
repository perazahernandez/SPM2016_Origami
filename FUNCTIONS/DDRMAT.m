% CODE FOR SIMULATION OF ORIGAMI SHEETS WITH SMOOTH ELASTIC FOLDS
% AUTHOR: Edwin A. Peraza Hernandez
% 05-16-2016
% PAPER DOI: http://dx.doi.org/10.1016/j.cad.2016.05.010

function DDR = DDRMAT(alpv, tetv, NF1, NF2)
% 2nd DERIVATIVE OF CONSTRAINT ROTATION MATRIX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUTS
%
% alpv: Array with face angles for the considered fold intersection
% tetv: Array with fold angles for the considered fold intersection
% NF1, NF2: Derivatives are taken with respect to the fold angles with these indices 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% OUTPUTS
% 
% DDR: 2nd DERIVATIVE OF CONSTRAINT ROTATION MATRIX
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R3 = @(t) [cos(t),  -sin(t),    0;
        sin(t),     cos(t),     0;
        0,      0,      1]; % MATRIX FOR ROTATION WITH RESPECT TO 3-AXIS
    
R1 = @(t) [1,   0,      0;
           0,   cos(t), -sin(t);
           0,   sin(t), cos(t)]; % MATRIX FOR ROTATION WITH RESPECT TO 1-AXIS

DR1 = @(t) [0,   0,      0;
           0,   -sin(t), -cos(t);
           0,   cos(t), -sin(t)]; % MATRIX FOR 1st DERIVATIVE OF ROTATION WITH RESPECT TO 1-AXIS
       
DDR1 = @(t) [0,   0,      0;
           0,   -cos(t), sin(t);
           0,   -sin(t), -cos(t)]; % MATRIX FOR 2nd DERIVATIVE OF ROTATION WITH RESPECT TO 1-AXIS
       
DDR = eye(3); % INITIALIZING MATRIX

for i = 1:size(alpv,2)
    
    if NF1 == NF2
        if i == NF1 
            DDR = DDR * DDR1(tetv(1,i))* R3(alpv(1,i));
        else
            DDR = DDR * R1(tetv(1,i))* R3(alpv(1,i));
        end        
    else
        if i == NF1 || i == NF2
            DDR = DDR * DR1(tetv(1,i))* R3(alpv(1,i));
        else
            DDR = DDR * R1(tetv(1,i))* R3(alpv(1,i));
        end
    end
end