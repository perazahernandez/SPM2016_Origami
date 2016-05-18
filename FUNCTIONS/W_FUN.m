% CODE FOR SIMULATION OF ORIGAMI SHEETS WITH SMOOTH ELASTIC FOLDS
% AUTHOR: Edwin A. Peraza Hernandez
% 05-16-2016
% PAPER DOI: http://dx.doi.org/10.1016/j.cad.2016.05.010

function We = W_FUN(I, EL, fo, tet)
% WORK EXERTED BY EXTERNALLY APPLIED LOADS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUTS
%
% fo: Current force vectors (3 x # forces)
% tet: Array with all fold angles at the current increment (N_F x 1)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% OUTPUTS
%
% We: Work exerted by externally applied loads
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


We = 0; % INITIALIZING

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADDING CONTRIBUTIONS FOR EACH LOAD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if norm(EL.gravity,2) == 0 % No gravity loads
    for i = 1:size(fo,2)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % OBTAINING DISPLACEMENT AT POINT OF APPLICATION
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        MF = MAP_FACE(I, tet, I.PFolds(EL.FACE0(1,i),1:end)'); % OBTAINING MAPPING 
        temp = MF * [EL.X0(1:3,i) ; 1];
        
        u = temp(1:3,1) - EL.X0(1:3,i); % DISPLACEMENT
        
        We = We  +   dot(u, fo(1:3,i)); % DETERMINING WORK
    end
else
    % SPEEDING UP IF GRAVITY IS CONSIDERED, CALCULATING MAPS ONLY ONCE
    MFS = zeros(4,4,I.N_P);
    
    for i = 1:I.N_P
        MFS(1:end,1:end,i) = MAP_FACE(I, tet, I.PFolds(i,1:end)');
    end
    
    
    for i = 1:size(fo,2)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % OBTAINING DISPLACEMENT AT POINT OF APPLICATION
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        temp = MFS(1:end,1:end,EL.FACE0(1,i)) * [EL.X0(1:3,i) ; 1];
        
        u = temp(1:3,1) - EL.X0(1:3,i); % DISPLACEMENT
        
        We = We  +   dot(u, fo(1:3,i)); % DETERMINING WORK
    end
    
end
