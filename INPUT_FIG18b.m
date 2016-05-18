% CODE FOR SIMULATION OF ORIGAMI SHEETS WITH SMOOTH ELASTIC FOLDS
% AUTHOR: Edwin A. Peraza Hernandez
% 05-16-2016
% PAPER DOI: http://dx.doi.org/10.1016/j.cad.2016.05.010

clear all; close all; clc % INITIALIZING ENVIRONMENT
% CODE FOR SIMULATION OF ORIGAMI SHEETS WITH SMOOTH ELASTIC FOLDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% MAIN STRUCTURES
% I: Structure with general inputs
% O: Structure with general outputs
% EL: Structure with elastic simulation inputs


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REQUIRED INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I.v: Vertex coordinates (I.N_V x 2)   | INPUT
% I.CM: Matrix with vertex indices of start/end points of fold vectors 
%       (I.N_F x 2) | INPUT
% I.w: Fold widths (I.N_F x 1) | INPUT
% I.d: d-parameters (perpendicular distance between vertex and fold boundary, 
%      1st column d1, 2nd column d2 (I.N_F x 2) | INPUT
% I.CF: Matrix with fold indices counterclockwise (row) associated with 
%       each interior fold intersection. Negative if fold vector points 
%       towards intersection (I.N_I x max(I.n_j)) | INPUT
% I.tet0: Initial fold angles in radians (I.N_F x 1)   | INPUT
% I.Dtet: Delta theta for derivatives using finite differences   | INPUT
% I.plotAng: Fold angles to plot (column array)  | INPUT
% I.plotInc: Increments to plot configuration (column array)   | INPUT
% I.FixedFace: Index of the face which is fixed | INPUT
% I.P_FoldVert: Folds neighboring faces (ccw order) (signed depending on
%           side neighboring sign
%           If extra vertices are added, the should be called as 
%           I.N_F + vertex #  (I.N_P x max(fold/vertex neighbors))   | INPUT
% I.saveGIF: Whether it should save a gif of the requested
%            outputs, 1: YES  Otherwise: NO   | INPUT
% I.ord : Order of geometric continuity | INPUT
% I.parallel : 1: Parallel computing YES, Otherwise: Parallel computing NO | INPUT

% EL.incs: # of load increments   | INPUT
% EL.E: Young's modulus of the folds (I.N_F x 1)   | INPUT
% EL.nu: Poisson's ratio of the folds (I.N_F x 1)   | INPUT
% EL.h: Thickness of the folds (I.N_F x 1)  | INPUT
% EL.f: Applied force vectors (3 x O.EL.N_f)    | INPUT
% EL.X0: Initial location of points where load is applied (3 x O.EL.N_f)    | INPUT
% EL.FACE0: Face index of points where load is applied (1 x O.EL.N_f)    | INPUT
% EL.lamR: Weight for rotation constraint matrix | INPUT
% EL.lamd: Weight for translation constraint vector | INPUT
% EL.tolr: Tolerance for the residual | INPUT
% EL.totdtet: Tolerance in change of fold angle | INPUT
% EL.maxiter: Maximum number of correction iterations | INPUT
% EL.maxdtet: Maximum change in fold angle per increment  | INPUT
% EL.XUP: Locations of points where displacement is requested (3 x #of requests)  | INPUT
% EL.FACEUP: Face index of location of points where displacement is requested (1 x #of requests)  | INPUT
% EL.gravity: Gravity vector (3 x 1), 
%             if equal to [0; 0; 0] then it is intended as no gravity  | INPUT
% EL.Pdensity: Density of the faces, can be left empty if gravity is not
%              considered (I.N_P x 1)    | INPUT
% EL.Fdensity: Density of the folds, can be left empty if gravity is not
%              considered (I.N_F x 1)    | INPUT
% EL.Ph: Thickness of the faces, can be left empty if gravity is not
%              considered (I.N_P x 1)  | INPUT


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS AND CALCULATED INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I.N_V: Number of vertices | CALCULATED
% I.N_F: Number of folds | CALCULATED
% I.N_I: Number of interior fold intersections | CALCULATED
% I.n_j: Number of folds associated with each interior fold intersection 
%        (I.N_I x 1) | CALCULATED
% I.m: Fold vectors (I.N_F x 2) | CALCULATED
% I.halp: Angle between fold vector and 1-axis (I.N_F x 1) | CALCULATED
% I.mI: Fold vectors per intersection (max(I.n_j) x 2*I.N_I) | CALCULATED
% I.alp: Face corner angles for intersections (I.N_I x max(I.n_j))  | CALCULATED
% I.N_P: Number of faces | CALCULATED 
% I.PFolds: Indices of folds crossed by the paths connecting the fixed face to each other
%           face (signed depending on the crossing direction)  (I.N_P x I.N_F)   | CALCULATED
% I.p1, I.p2, I.p3, I.p4: Corner points of folds | CALCULATED
% I.Pcorners: Face corner points (I.N_P x 2*2*size(I.P_FoldVert,2)) | CALCULATED
% I.Pnpoints: Number of face corner poitns (I.N_P x 1) | CALCULATED
% I.PconF: Index of face neighboring the p1-p2 side of the fold (I.N_F x
%            1)   | CALCULATED
% I.Parea: Surface area of the faces (I.N_P x 1)   | CALCULATED
% I.Pcentroid: Centroid point of the faces (I.N_P x 2)   | CALCULATED
% I.Farea: Surface area of the folds (I.N_F x 1)   | CALCULATED
% I.Fcentroid: Centroid point of the folds (I.N_F x 2)   | CALCULATED
% I.overlap: 1 folds overlap, 0 folds do not overlap | CALCULATED
% I.overlap_Pairs: Fold pairs that overlap (ignore first row of 2
%                  zeros) | CALCULATED
% I.overlap_nPairs: Number of fold pairs that overlap | CALCULATED
% I.overlap_stimated_l: Stimated overlap length (minimum of non-zero) | CALCULATED 
% I.Li: Length of fold  (I.N_F x 1)   | CALCULATED
% I.Ii: Area moment of inertia of folds  (I.N_F x 1)   | CALCULATED
% I.unitmi: Unit vector in the direction of mi  (I.N_F x 2)   | CALCULATED

% O.ninc: Total number of increments (for any process)
% O.tet: Fold angles (I.N_F x O.ninc+1)
% O.totaltime: Total time of running code

% O.EL.N_f: Number of applied loads
% O.EL.rnormal: 2-norm of residual vector at the last iteration
% O.EL.dtetnormal: 2-norm of fold angle increment vector last iteration
% O.EL.exitFlag: exitFlag (0: max#iterations), (1: r below tol), 
%                (2: dtet below tol)
% O.EL.niter: Number of iterations per increment
% O.EL.res: Norm of translation constraint vector (O.ninc + 1, I.N_I)
% O.EL.fo: Force vector for each increment (3*(O.ninc+1) x O.EL.N_f)
% O.EL.xf: Location of the points where point loads applied (3*(O.ninc+1) x O.EL.N_f)
% O.EL.UP: Displacements of points of interest (3*(# of requested points) x (O.ninc+1))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEGINNING OF INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I.plotAng = []';

I.v = 2*0.01*[-20 -10 0 10 20, -20 -10 10 20 -20, 20 -20 -10 10 20, -20 -10 0 10 20
    -20 -20 -20 -20 -20, -10 -10 -10 -10 0, 0 10 10 10 10, 20 20 20 20 20]';

I.CM = [7 7 7 8 8, 8 7 7 8 7, 7 8 8 13 14, 13 13 14 13 13, 13 14 14 14
    1 2 3 3 4, 5 6 8 9 10, 13 14 11 10 11, 12 14 15 16 17, 18 18 19 20]';  

I.w = 2*0.01*2*[1 1 1 1 1, 1 1 1 1 1, 1 1 1 1 1, 1 1 1 1 1, 1 1 1 1]';

I.d = 2*0.01*2*[2 2 2 2 2, 2 2 2 2 2, 2 2 2 2 2, 2 2 2 2 2, 2 2 2 2
    2 1 2 2 1, 2 1 2 1 2, 2 2 2 2 2, 1 2 1 2 1, 2 2 1 2]';
    
I.CF = [8 11 10 7 1 2 3
    9 13 12 -8 4 5 6
    17 21 20 19 16 14 -11
    18 24 23 22 -17 -12 15];

I.tet0 = zeros(1,24)';

I.Dtet = 4*pi/180;

I.FixedFace = 9;

I.P_FoldVert = [ 1 -2 0 0
    2 -3 0 0
    4 -5 0 0
    5 -6 0 0
    7 -1 0 0;
    3 -4 -8 0
    6 -9 0 0
    10 -7 0 0
    -11 8 12 -17
    -13 9 0 0;
    14 -10 11 0
    -12 13 -15 0
    16 -14 0 0
    -18 15 0 0
    19 -16 0 0;
    -21 17 22 0
    -24 18 0 0
    20 -19 0 0
    21 -20 0 0
    23 -22 0 0;
    24 -23 0 0];
    
I.ord = 1;

I.saveGIF = 0;

I.parallel = 0;


EL.incs = 80;
I.plotInc = [0 40 80]';


EL.tolr = 1e-6;
EL.totdtet = 1e-4;
EL.maxiter = 40;

EL.maxdtet = 15*pi/180;

EL.lamR = 100;
EL.lamd = 1e0;

EL.E = [15.2e6    15.2e6    15.2e6    15.2e6    15.2e6, ...
    15.2e6    15.2e6    15.2e6    15.2e6    15.2e6,...
    15.2e6    15.2e6    15.2e6    15.2e6    15.2e6,...
    15.2e6    15.2e6    15.2e6    15.2e6    15.2e6,...
    15.2e6    15.2e6    15.2e6    15.2e6    ]';

EL.nu = [0.45   0.45    0.45    0.45    0.45, ...
    0.45   0.45    0.45    0.45    0.45,...
    0.45   0.45    0.45    0.45    0.45,...
    0.45   0.45    0.45    0.45    0.45,...
    0.45   0.45    0.45    0.45    ]';

EL.h = [1.8e-3    1.8e-3    1.8e-3    1.8e-3    1.8e-3,...
    1.8e-3    1.8e-3    1.8e-3    1.8e-3    1.8e-3,...
    1.8e-3    1.8e-3    1.8e-3    1.8e-3    1.8e-3,...
    1.8e-3    1.8e-3    1.8e-3    1.8e-3    1.8e-3,...
    1.8e-3    1.8e-3    1.8e-3    1.8e-3    ]';

EL.f = [-2  0 2  0  
    0  -2   0  2    
    -1.2  -1.2   -1.2  -1.2 ];

EL.X0 = 2*[0.1646    0   -0.1646   0
    0   0.1646  0  -0.1646
    0   0  0   0];

EL.FACE0 = [12    16   11   6];

EL.XUP = 2*[ 0
     -0.1646
     0];

EL.FACEUP = [6];

EL.gravity = [0;    0;  0];

EL.Pdensity = [1040     1040    1040    1040    1040 ...
    1040     1040    1040    1040    1040 ...
    1040     1040    1040    1040    1040 ...
    1040     1040    1040    1040    1040 ...
    1040     ]';
	
EL.Ph = [2e-3   2e-3   2e-3    2e-3     2e-3 ...
    2e-3   2e-3   2e-3    2e-3     2e-3 ...
    2e-3   2e-3   2e-3    2e-3     2e-3 ...
    2e-3   2e-3   2e-3    2e-3     2e-3 ...
    2e-3   ]'; 

EL.Fdensity = [1200     1200    1200	1200	1200	1200	1200	1200 ...
		1200     1200    1200	1200	1200	1200	1200	1200 ...
        1200     1200    1200	1200	1200	1200	1200	1200 ]';
    
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALLING MAIN_MODULE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[I, O] = MAIN_MODULE(I, EL);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVING RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('RESULTS')


% PLOTTING 
figure()

for i = 1:(size(O.EL.fo,1)/3)
    u2_plot(i) =  O.EL.UP(2,i)*100; % Out-of-plane displacement u2
    u3_plot(i) =  O.EL.UP(3,i)*100; % Out-of-plane displacement u3
end
plot(u2_plot,u3_plot,...
    'LineWidth',2)
xlabel('u2 [cm]')
ylabel('u3 [cm]')
set(gca,'FontName','Times New Roman','fontsize', 20,'linewidth',1.15)
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca,'ticklength',3*get(gca,'ticklength'))