% CODE FOR SIMULATION OF ORIGAMI SHEETS WITH SMOOTH ELASTIC FOLDS
% AUTHOR: Edwin A. Peraza Hernandez
% 05-16-2016
% PAPER DOI: http://dx.doi.org/10.1016/j.cad.2016.05.010

function FSP = FOLD_SURF_POINTS(I, tet, FN, MF)
% POINTS TO PLOT FOLDS USING THE surf FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUTS
%
% tet: Array with all fold angles at the current increment (N_F x 1)
% FN: Fold number
% MF: Linear mapping from points in local coordinate system to global coords
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% OUTPTUS
%
% FSP: Fold surface points
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

npoints1 = 20; % Number of points through curved cross-section direction
npoints2 = 2; % Number of points through straight direction 

FSP = zeros(npoints1, npoints2, 3); % INITIALIZING

if I.ord == 2 % 2nd ORDER OF GEOMETRIC CONTINUITY
    % BASIS FUNCTIONS
    h55 = @(t) ((1/16)*t^5+(1/16)*t^4-(1/8)*t^3-(1/8)*t^2+(1/16)*t+1/16);
    h54 = @(t) (-(1/16)*t^5+(1/16)*t^4+(1/8)*t^3-(1/8)*t^2-(1/16)*t+1/16);
    h53 = @(t) (-(3/16)*t^5-(1/16)*t^4+(5/8)*t^3+(3/8)*t^2-(7/16)*t-5/16);
    h52 = @(t) (-(3/16)*t^5+(1/16)*t^4+(5/8)*t^3-(3/8)*t^2-(7/16)*t+5/16);
    h51 = @(t) ((3/16)*t^5-(5/8)*t^3+(15/16)*t+1/2);
    h50 = @(t) (-(3/16)*t^5+(5/8)*t^3-(15/16)*t+1/2);
    
    % OBTAINING FOLD SHAPE PARAMETERS
    wstr = WSTR(tet(FN), I.w(FN), I.ord);
    c1 = C1(tet(FN), I.w(FN), I.ord);
    c2 = C2(tet(FN), I.w(FN));
    
    cm1 = [0;   -wstr/2;    0];
    cp1 = [0;   wstr/2;    0];
    dcm1 = c1*[0;   cos(tet(FN)/2);    -sin(tet(FN)/2)];
    dcp1 = c1*[0;   cos(tet(FN)/2);    sin(tet(FN)/2)];
    ddcm1 = c2*[0;   cos(tet(FN)/2);    -sin(tet(FN)/2)];
    ddcp1 = c2*[0;   -cos(tet(FN)/2);    -sin(tet(FN)/2)];
    h = [norm(I.m(FN,1:2)) - I.d(FN,1) - I.d(FN,2); 0 ; 0];
    
    % APPLYING TRANSFORMATIONS
    temp = MF*[cm1;1];    cm1 = temp(1:3);
    temp = MF*[cp1;1];    cp1 = temp(1:3);
    dcm1 = MF(1:3,1:3)*dcm1;
    dcp1 = MF(1:3,1:3)*dcp1;
    ddcm1 = MF(1:3,1:3)*ddcm1;
    ddcp1 = MF(1:3,1:3)*ddcp1;
    h = MF(1:3,1:3)*h; 
    
    z1v = linspace(-1,1,npoints1)';
    z2v = linspace(0,1,npoints2)';
    
    for i = 1:size(z1v,1);
        for j = 1:size(z2v,1);
            z1 = z1v(i);   z2 = z2v(j);
            FSP(i,j,1:3) = ...
                h50(z1)*cm1 + h51(z1)*cp1 ...
                + h52(z1)*dcm1 + h53(z1)*dcp1 ...
                + h54(z1)*ddcm1 + h55(z1)*ddcp1 ...
                + z2*h;
            
        end
    end
    
elseif I.ord == 1 % 1st ORDER OF GEOMETRIC CONTINUITY
    % BASIS FUNCTIONS
    h33 = @(t) ((1/4)*t^3+(1/4)*t^2-(1/4)*t-1/4);
    h32 = @(t) ((1/4)*t^3-(1/4)*t^2-(1/4)*t+1/4);
    h31 = @(t) (-(1/4)*t^3+(3/4)*t+1/2);
    h30 = @(t) ((1/4)*t^3-(3/4)*t+1/2);
    
    % OBTAINING FOLD SHAPE PARAMETERS
    wstr = WSTR(tet(FN), I.w(FN), I.ord);
    c1 = C1(tet(FN), I.w(FN), I.ord);
    
    cm1 = [0;   -wstr/2;    0];
    cp1 = [0;   wstr/2;    0];
    dcm1 = c1*[0;   cos(tet(FN)/2);    -sin(tet(FN)/2)];
    dcp1 = c1*[0;   cos(tet(FN)/2);    sin(tet(FN)/2)];
    h = [norm(I.m(FN,1:2)) - I.d(FN,1) - I.d(FN,2); 0 ; 0];
    
    % APPLYING TRANSFORMATIONS
    temp = MF*[cm1;1];    cm1 = temp(1:3);
    temp = MF*[cp1;1];    cp1 = temp(1:3);
    dcm1 = MF(1:3,1:3)*dcm1;
    dcp1 = MF(1:3,1:3)*dcp1;
    h = MF(1:3,1:3)*h;
    
    z1v = linspace(-1,1,npoints1)';
    z2v = linspace(0,1,npoints2)';
    
    for i = 1:size(z1v,1);
        for j = 1:size(z2v,1);
            z1 = z1v(i);   z2 = z2v(j);
            FSP(i,j,1:3) = ...
                h30(z1)*cm1 + h31(z1)*cp1 ...
                + h32(z1)*dcm1 + h33(z1)*dcp1 ...
                + z2*h;
            
        end
    end
        
end




