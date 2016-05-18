% CODE FOR SIMULATION OF ORIGAMI SHEETS WITH SMOOTH ELASTIC FOLDS
% AUTHOR: Edwin A. Peraza Hernandez
% 05-16-2016
% PAPER DOI: http://dx.doi.org/10.1016/j.cad.2016.05.010

function d = dVEC(alpv, tetv, twv, lvecs, ord)
% TRANSLATION CONSTRAINT VECTOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUTS
%
% alpv: Array with face angles for the considered fold intersection
% tetv: Array with fold angles for the considered fold intersection
% twv: Array with fold widths for the considered fold intersection
% lvecs: order of geometric continuity
% ord: Order of geometric continuity
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% OUTPUTS
% 
% d: Translation constraint vector
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R3 = @(t) [cos(t),  -sin(t),    0;
    sin(t),     cos(t),     0;
    0,      0,      1]; % MATRIX FOR ROTATION WITH RESPECT TO 3-AXIS

R1 = @(t) [1,   0,      0;
    0,   cos(t), -sin(t);
    0,   sin(t), cos(t)]; % MATRIX FOR ROTATION WITH RESPECT TO 1-AXIS


d = zeros(3,1); % INITIALIZING VECTOR 

for i = 1:size(alpv,2)
    
    ljk = [lvecs(i,1), lvecs(i,2), 0]'; % OBTAINING ljk VECTOR

    wsjk = [0
        WSTR(tetv(1,i), twv(1,i), ord)
        0]; % OBTAINING wjk VECTOR
    
    R = eye(3); % INITIALIZING ROTATION MATRIX
	
    for l = 1:i-1
        R =  R * R1(tetv(1,l))* R3(alpv(1,l));
    end
    
    d = d + (R * R1(0.5*tetv(1,i)) * wsjk) + ...
        (R * R1(tetv(1,i)) * ljk); 
end

