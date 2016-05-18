% CODE FOR SIMULATION OF ORIGAMI SHEETS WITH SMOOTH ELASTIC FOLDS
% AUTHOR: Edwin A. Peraza Hernandez
% 05-16-2016
% PAPER DOI: http://dx.doi.org/10.1016/j.cad.2016.05.010

function MF = MAP_FACE(I, tet, Fs)
% RETURNS 4 x 4 FACE MAPPING MATRIX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUTS
%
% tet: Array with all fold angles at the current increment (N_F x 1)
% Fs: List of fold numbers (signed) crossed by path connecting fixed
%        face to considered face (column array)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% OUTPUTS
%
% MF: 4 x 4 FACE MAPPING MATRIX
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     

TR_MAT = @(b1,b2,b3,a,t) [-2*cos(t)^2*cos(a)^2+2*cos(t)^2+2*cos(a)^2-1,...  % 1,1
    2*sin(t)^2*sin(a)*cos(a),...
    2*sin(a)*sin(t)*cos(t),...
    2*cos(t)^2*b3(2)*sin(a)*cos(a)+2*cos(t)^2*b3(1)*cos(a)^2-2*sin(t)*cos(t)*b3(3)*sin(a)-cos(t)*b2(1)*cos(a)^2-sin(a)*cos(a)*b2(2)*cos(t)+sin(a)*cos(a)*b1(2)*cos(t)+cos(t)*b1(1)*cos(a)^2+sin(a)*sin(t)*b2(3)-sin(a)*sin(t)*b1(3)-2*cos(t)^2*b3(1)+b2(1)*cos(a)^2+sin(a)*cos(a)*b2(2)-2*sin(a)*cos(a)*b3(2)-2*b3(1)*cos(a)^2-sin(a)*cos(a)*b1(2)-b1(1)*cos(a)^2+cos(t)*b2(1)-cos(t)*b1(1)+b3(1)+b1(1);
    2*sin(t)^2*sin(a)*cos(a), ... % 2,1
    2*cos(t)^2*cos(a)^2-2*cos(a)^2+1,...
    -2*cos(a)*sin(t)*cos(t), ...
    -2*cos(t)^2*b3(2)*cos(a)^2+2*cos(t)^2*b3(1)*sin(a)*cos(a)+2*sin(t)*cos(t)*b3(3)*cos(a)-cos(t)*b2(1)*sin(a)*cos(a)+cos(t)*b2(2)*cos(a)^2-cos(t)*b1(2)*cos(a)^2+cos(t)*b1(1)*sin(a)*cos(a)-cos(a)*sin(t)*b2(3)+cos(a)*sin(t)*b1(3)+cos(a)*sin(a)*b2(1)-b2(2)*cos(a)^2+2*cos(a)^2*b3(2)-2*cos(a)*sin(a)*b3(1)+b1(2)*cos(a)^2-cos(a)*sin(a)*b1(1)+b2(2)-b3(2);
    -2*sin(a)*sin(t)*cos(t),... % 3,1
    2*cos(a)*sin(t)*cos(t), ...
    2*cos(t)^2-1,...
    -2*sin(t)*cos(t)*b3(2)*cos(a)+2*sin(t)*cos(t)*b3(1)*sin(a)-sin(a)*sin(t)*b2(1)+cos(a)*sin(t)*b2(2)-cos(a)*sin(t)*b1(2)+sin(a)*sin(t)*b1(1)-2*cos(t)^2*b3(3)+cos(t)*b2(3)-cos(t)*b1(3)+b3(3)+b1(3);
    0,0,0,1]; % FOLD TRANSFORMATION MATRIX

MF = eye(4); % INITIALIZING FACE MAPPING MATRIX

for i = 1:size(Fs,1)
    k = abs(Fs(i));
    if Fs(i) == 0
        break
    elseif Fs(i) > 0 % POSITIVELY CROSSED FOLD
        
        MF = MF * TR_MAT( [I.p1(k, 1:2)'; 0] , ...
            [I.p4(k, 1:2)'; 0] - (I.w(k) - WSTR(tet(k), I.w(k), I.ord)) *[-I.unitmi(k,2); I.unitmi(k,1); 0] , ...
            [I.p4(k, 1:2)'; 0],... 
            I.halp(k), 0.5 * tet(k) );
        
    elseif Fs(i) < 0 % NEGATIVELY CROSSED FOLD

        MF = MF * TR_MAT( [I.p4(k, 1:2)'; 0] ,  ...
            [I.p1(k, 1:2)'; 0] + (I.w(k) - WSTR(tet(k), I.w(k), I.ord)) *[-I.unitmi(k,2); I.unitmi(k,1); 0], ...
            [I.p1(k, 1:2)'; 0], ...
            I.halp(k), -0.5 * tet(k) );        
        
    end

end