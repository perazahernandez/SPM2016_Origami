% CODE FOR SIMULATION OF ORIGAMI SHEETS WITH SMOOTH ELASTIC FOLDS
% AUTHOR: Edwin A. Peraza Hernandez
% 05-16-2016
% PAPER DOI: http://dx.doi.org/10.1016/j.cad.2016.05.010

function O = ELASTIC_MODULE_P(I, EL)
% SIMULATION OF AN ORIGAMI SHEET WITH SMOOTH ELASTIC FOLDS 
% (PARALLEL for LOOPS FOR JACOBIAN OPTION)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Quasi-Newton parameter (how many increments to update Jacobian)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
number_mode = 1;

% INITIALIZING
O.ninc = EL.incs;
O.tet = zeros(I.N_F, O.ninc+1);
O.EL.rnormal = zeros(O.ninc+1,1);
O.EL.dtetnormal = zeros(O.ninc+1,1);
O.EL.exitFlag = zeros(O.ninc+1,1);
O.EL.niter = zeros(O.ninc+1,1);
O.EL.dres = zeros(O.ninc+1,I.N_I);

nloads = size(EL.f,2); % NUMBER OF APPLIED LOADS (NOT COUNTING GRAVITY)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% APPLYING GRAVITY LOADS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if norm(EL.gravity,2) == 0 % NO GRAVITY LOADS
else 
    for i = 1:I.N_P
        if i == 1  && isempty(EL.f) == 1
            EL.X0 = [I.Pcentroid(i,1); I.Pcentroid(i,2); 0];
            EL.FACE0 = [i];
            EL.f = EL.Pdensity(i) * EL.Ph(i) * I.Parea(i) * ...
                EL.gravity;
        else
            EL.X0 = [EL.X0, [I.Pcentroid(i,1); I.Pcentroid(i,2); 0]];
            EL.FACE0 = [EL.FACE0, i];
            EL.f = [EL.f, EL.Pdensity(i) * EL.Ph(i) * I.Parea(i) * ...
                EL.gravity];
        end
    end
    for i = 1:I.N_F
        EL.X0 = [EL.X0, [I.Fcentroid(i,1); I.Fcentroid(i,2); 0]];
        EL.FACE0 = [EL.FACE0, I.PconF(i)];
        EL.f = [EL.f, EL.Fdensity(i) * EL.h(i) * I.Farea(i) * ...
            EL.gravity];
        
    end
end


O.EL.N_f = size(EL.f,2); % NUMBER OF APPLIED FORCES

% STORING DATA OF FORCE VECTORS AND LOCATIONS
O.EL.fo = zeros(3*(O.ninc+1), O.EL.N_f);
O.EL.xf = zeros(3*(O.ninc+1), O.EL.N_f);
O.EL.UP = zeros(3*size(EL.FACEUP,2),(O.ninc+1));

% OBTAINING WHICH FOLDS BELONG TO THE MAPPING OF ANY CONSIDERED FACE UNDER LOADING
thetaMap = zeros(I.N_F,1); % Fold angles that are active in mappings of applied load faces
thetaMapC = zeros(I.N_F,I.N_F); % If Fold angles are coupled in mappings
for i = 1:size(I.PFolds,1)
    for j = 1:size(I.PFolds,2)
        if I.PFolds(i,j) ~= 0
            thetaMap(abs(I.PFolds(i,j))) = 1;
            
            for k = 1:size(I.PFolds,2)
                if I.PFolds(i,k) ~= 0
                    thetaMapC(abs(I.PFolds(i,j)),abs(I.PFolds(i,k))) = 1;
                end
            end
        end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INCREMENT LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:(O.ninc+1)
    disp(['Increment: ',num2str(i-1), ' out of ',num2str(O.ninc)])
    
    if i == 1 % INITIAL STATE
        O.tet(1:end,1) = I.tet0;
    elseif i == 2
        O.tet(1:end,i) = O.tet(1:end,i-1);
    elseif i == 3 % Use speed up constant
        O.tet(1:end,i) = O.tet(1:end,i-1) + ...
            (O.tet(1:end,i-1)-O.tet(1:end,i-2));
    elseif i == 4 % Use speed up linear
        O.tet(1:end,i) = O.tet(1:end,i-1) + ...
            2*(O.tet(1:end,i-1)-O.tet(1:end,i-2)) - ...
            (O.tet(1:end,i-2)-O.tet(1:end,i-3));
    else % Use speed up quadratic
        O.tet(1:end,i) = O.tet(1:end,i-1) + ...
            3*(O.tet(1:end,i-1)-O.tet(1:end,i-2)) - ...
            3*(O.tet(1:end,i-2)-O.tet(1:end,i-3)) + ...
            (O.tet(1:end,i-3)-O.tet(1:end,i-4));
    end
    
    
    if norm(EL.gravity,2) == 0 % LOADING UNDER NO GRAVITY
        
        % UNCOMMENT IF LOADING ONLY
        fo = EL.f * ( (1/O.ninc)*i - (1/O.ninc) );
        
        
        % UNCOMMENT IF LOADING/UNLOADING
        %     if i~=1
        %         if i <= O.ninc/2+1
        %             fo = fo + 2*EL.f *(1/O.ninc);
        %         else
        %             fo = fo - 2*EL.f *(1/O.ninc);
        %         end
        %     else
        %         fo = 0*EL.f;
        %     end
    
    else % APPLYING GRAVITY BEFORE LOADING
        if i~=1
            if i <= O.ninc/2+1 % APPLY GRAVITY LOAD
                fo(1:3,(nloads+1):end) = fo(1:3,(nloads+1):end) + 2*EL.f(1:3,(nloads+1):end)*(1/O.ninc);
            else % APPLY OTHER LOADS
                fo(1:3,1:nloads) = fo(1:3,1:nloads) + 2*EL.f(1:3,1:nloads)*(1/O.ninc);
            end
        else
            fo = 0*EL.f;
        end
    end
    
    O.EL.fo( (3*i-2):3*i, 1:end) = fo; % SAVING CURRENT LOADS
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ITERATION LOOP
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for it = 1:EL.maxiter
        
        O.EL.niter(i,1) = it;
        

        r = zeros(I.N_F, 1); % RESIDUAL VECTOR
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SPEED UP, CALCULATE JACOBIAN ONLY EVERY SELECT # OF ITERS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if mod(it,number_mode) == 1 || number_mode == 1
            drdt = zeros(I.N_F, I.N_F); % JACOBIAN
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % CONTRIBUTIONS TO RESIDUAL VECTOR FROM CONSTRAINTS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % FOLD ANGLES ASSOCIATED WITH EACH FOLD INTERSECTION
        ttet = zeros(I.N_I,max(I.n_j));
        tw = zeros(I.N_I,max(I.n_j));
        
        % MAPPING FOLD ANGLES AND FOLD WIDTHS ASSOCIATED WITH EACH FOLD INTERSECTION
        for j = 1:I.N_I
            for k = 1:I.n_j(j,1)
                ttet(j,k) = O.tet(abs(I.CF(j,k)),i);
                tw(j,k) = I.w(abs(I.CF(j,k)));
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % OBTAINING R AND d, AND DERIVATIVES
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for j = 1:I.N_I
            % FACE ANGLES OF THE CONSIDERED INTERSECTION
            alpv = I.alp(j,1:I.n_j(j,1));
            % FOLD ANGLES OF THE CONSIDERED INTERSECTION
            tetv = ttet(j,1:I.n_j(j,1));
            twv = tw(j,1:I.n_j(j,1));
            R = RMAT(alpv, tetv);
            d = dVEC(alpv, tetv, twv, I.lvecs(1:end,(2*j-1):2*j), I.ord);
            
            O.EL.dres(i,j) = norm(d,2); % Saving norm of d-vector
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % FIRST DERIVATIVES OF RESIDUAL
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for k = 1:I.n_j(j)
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % CALCULATING FIRST DERIVATIVES theta_k  (DRk, Ddk)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                tetvM = tetv;
                tetvM(1,k) = tetvM(1,k) - I.Dtet;
                dM = dVEC(alpv, tetvM, twv, I.lvecs(1:end,(2*j-1):2*j), I.ord);
                
                tetvP = tetv;
                tetvP(1,k) = tetvP(1,k) + I.Dtet;
                dP = dVEC(alpv, tetvP, twv, I.lvecs(1:end,(2*j-1):2*j), I.ord);
                
                % COMPONENT-WISE DERIVATIVES
                DRk = DRMAT(alpv, tetv, k); 
                Ddk = (-0.5*dM + 0.5*dP) / I.Dtet;
                
                % ADDING TO RESIDUAL VECTOR
                r(abs(I.CF(j,k))) = r(abs(I.CF(j,k))) + ...
                    EL.lamR * DRk(2,3)* R(2,3);
                r(abs(I.CF(j,k))) = r(abs(I.CF(j,k))) + ...
                    EL.lamR * DRk(3,1)* R(3,1);
                r(abs(I.CF(j,k))) = r(abs(I.CF(j,k))) + ...
                    EL.lamR * DRk(1,2)* R(1,2);
                
                r(abs(I.CF(j,k))) = r(abs(I.CF(j,k))) + ...
                    EL.lamd * Ddk(1)* d(1);
                r(abs(I.CF(j,k))) = r(abs(I.CF(j,k))) + ...
                    EL.lamd * Ddk(2)* d(2);
                r(abs(I.CF(j,k))) = r(abs(I.CF(j,k))) + ...
                    EL.lamd * Ddk(3)* d(3);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % SECOND DERIVATIVES OF RESIDUAL
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % SPEED UP, CALCULATE JACOBIAN ONLY EVERY SELECT # OF ITERS
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if mod(it,number_mode) == 1 || number_mode == 1
                    
                    tempJ = zeros(I.n_j(j),1);
                    
                    parfor l = 1:I.n_j(j)
                        if l < k
                        else
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % CALCULATING FIRST DERIVATIVE theta_l (DRl, Ddl)
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            tetvM = tetv;
                            tetvM(1,l) = tetvM(1,l) - I.Dtet;
                            dM = dVEC(alpv, tetvM, twv, I.lvecs(1:end,(2*j-1):2*j), I.ord);
                            
                            tetvP = tetv;
                            tetvP(1,l) = tetvP(1,l) + I.Dtet;
                            dP = dVEC(alpv, tetvP, twv, I.lvecs(1:end,(2*j-1):2*j), I.ord);
                            
                            % COMPONENT-WISE DERIVATIVES
                            DRl = DRMAT(alpv, tetv, l);
                            Ddl = (-0.5*dM + 0.5*dP) / I.Dtet;
                            
                            tempJ(l) = tempJ(l) + EL.lamR * DRk(2,3)* DRl(2,3) + ...
                                EL.lamR * DRk(3,1)* DRl(3,1) + ...
                                EL.lamR * DRk(1,2)* DRl(1,2) + ...
                                EL.lamd * Ddk(1) * Ddl(1) + ...
                                EL.lamd * Ddk(2) * Ddl(2) + ...
                                EL.lamd * Ddk(3) * Ddl(3);
                            
                            
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % CALCULATING SECOND DERIVATIVES theta_k theta_l
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            tetvM = tetv;
                            tetvM(1,k) = tetvM(1,k) - I.Dtet;
                            tetvM(1,l) = tetvM(1,l) - I.Dtet;
                            dM = dVEC(alpv, tetvM, twv, I.lvecs(1:end,(2*j-1):2*j), I.ord);
                            
                            tetvP = tetv;
                            tetvP(1,k) = tetvP(1,k) + I.Dtet;
                            tetvP(1,l) = tetvP(1,l) - I.Dtet;
                            dP = dVEC(alpv, tetvP, twv, I.lvecs(1:end,(2*j-1):2*j), I.ord);
                            
                            % COMPONENT-WISE DERIVATIVES
                            Ddkm = (-0.5*dM + 0.5*dP) / I.Dtet; 
                            
                            tetvM = tetv;
                            tetvM(1,k) = tetvM(1,k) - I.Dtet;
                            tetvM(1,l) = tetvM(1,l) + I.Dtet;
                            dM = dVEC(alpv, tetvM, twv, I.lvecs(1:end,(2*j-1):2*j), I.ord);
                            
                            tetvP = tetv;
                            tetvP(1,k) = tetvP(1,k) + I.Dtet;
                            tetvP(1,l) = tetvP(1,l) + I.Dtet;
                            dP = dVEC(alpv, tetvP, twv, I.lvecs(1:end,(2*j-1):2*j), I.ord);
                            
                            % COMPONENT-WISE DERIVATIVES
                            Ddkp = (-0.5*dM + 0.5*dP) / I.Dtet; 
                            
                            DDRkl = DDRMAT(alpv, tetv, k, l); 
                            DDdkl = (-0.5*Ddkm + 0.5*Ddkp) / I.Dtet;
                            
                            tempJ(l) = tempJ(l) + EL.lamR * DDRkl(2,3)* R(2,3) + ...
                                EL.lamR * DDRkl(3,1)* R(3,1) + ...
                                EL.lamR * DDRkl(1,2)* R(1,2) + ...
                                EL.lamd * DDdkl(1) * d(1) + ...
                                EL.lamd * DDdkl(2) * d(2) + ...
                                EL.lamd * DDdkl(3) * d(3);
                           
                        end
                    end
                    
                    for l = 1:I.n_j(j)
                        if k == l
                            drdt(abs(I.CF(j,k)),abs(I.CF(j,l))) = ...
                                drdt(abs(I.CF(j,k)),abs(I.CF(j,l))) ...
                                + tempJ(l);
                        else
                            drdt(abs(I.CF(j,k)),abs(I.CF(j,l))) = ...
                                drdt(abs(I.CF(j,k)),abs(I.CF(j,l))) ...
                                + tempJ(l);
                            drdt(abs(I.CF(j,l)),abs(I.CF(j,k))) = ...
                                drdt(abs(I.CF(j,k)),abs(I.CF(j,l))) ;
                        end
                    end
                    
                    
                    
                end
            end
            
        end % END OF FOLD INTERSECTIONS LOOP
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% CONTRIBUTIONS TO RESIDUAL VECTOR FROM TOTAL POTENTIAL
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
        tetv = O.tet(1:end, i);
        

        for j = 1:I.N_F    
            % FIRST AND SECOND DERIVATIVES OF U and We
            tetvM = tetv;
            tetvM(j) = tetvM(j) - I.Dtet;
            tetvP = tetv;
            tetvP(j) = tetvP(j) + I.Dtet;
            
            UiM = Ui_FUN(I, EL, tetvM, j);
            
            UiP = Ui_FUN(I, EL, tetvP, j);
            
            r(j) = r(j) + ...
                (-0.5*UiM + 0.5*UiP)/I.Dtet;
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % SPEED UP, CALCULATE JACOBIAN ONLY EVERY SELECT # OF ITERS
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if mod(it,number_mode) == 1 || number_mode == 1
                Ui =  Ui_FUN(I, EL, tetv, j);
                drdt(j,j) = drdt(j,j) + ...
                    (UiM - 2*Ui + UiP)/(I.Dtet^2);

            end

            
            % OBTAINING DERIVATIVES ONLY IF theta_i BELONGS TO THE MAPPING OF
            % ANY CONSIDERED FACE UNDER LOADING
            if thetaMap(j,1) ~= 0
                tetvM = tetv;
                tetvM(j) = tetvM(j) - I.Dtet;
                tetvP = tetv;
                tetvP(j) = tetvP(j) + I.Dtet;
                WeM = W_FUN(I, EL, fo, tetvM);
                
                WeP = W_FUN(I, EL, fo, tetvP);
                
                r(j) = r(j) - ...
                    (-0.5*WeM + 0.5*WeP)/I.Dtet;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % SPEED UP, CALCULATE JACOBIAN ONLY EVERY SELECT # OF ITERS
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if mod(it,number_mode) == 1 || number_mode == 1
                    We = W_FUN(I, EL, fo, tetv);
                    drdt(j,j) = drdt(j,j) - ...
                        (WeM - 2*We + WeP)/(I.Dtet^2);
                    
                    tempU = zeros(I.N_F,I.N_F);
                    
                    parfor k = 1:I.N_F
                        
                        if k <= j || thetaMapC(j,k) == 0 
                        else
                            tetvM = tetv;
                            tetvM(k) = tetvM(k) - I.Dtet;
                            tetvM(j) = tetvM(j) - I.Dtet;
                            tetvP = tetv;
                            tetvP(k) = tetvP(k) + I.Dtet;
                            tetvP(j) = tetvP(j) - I.Dtet;
                            
                            WeM = W_FUN(I, EL, fo, tetvM);
                            WeP = W_FUN(I, EL, fo, tetvP);
                            
                            DWeM = (-0.5*WeM + 0.5*WeP)/I.Dtet;
                            
                            tetvM = tetv;
                            tetvM(k) = tetvM(k) - I.Dtet;
                            tetvM(j) = tetvM(j) + I.Dtet;
                            tetvP = tetv;
                            tetvP(k) = tetvP(k) + I.Dtet;
                            tetvP(j) = tetvP(j) + I.Dtet;
                            
                            WeM = W_FUN(I, EL, fo, tetvM);
                            WeP = W_FUN(I, EL, fo, tetvP);
                            
                            DWeP = (-0.5*WeM + 0.5*WeP)/I.Dtet;
                            
                            DDWe = (-0.5*DWeM + 0.5*DWeP)/I.Dtet;
                            tempU(j,k) = tempU(j,k) - DDWe;
                            
                        end
                    end 
                    
                    for k = 1:I.N_F
                        drdt(j,k) = drdt(j,k) - tempU(j,k);
                        drdt(k,j) = drdt(j,k);
                    end
                    
                    
                end
            end
        end 

        O.EL.rnormal(i,1) = norm(r,2)/size(r,1);
        disp(['    Iteration: ',num2str(it), ', max: ',num2str(EL.maxiter), '  | Residual: ', num2str(O.EL.rnormal(i,1))])
        %disp(['    Iteration: ',num2str(it), ', max: ',num2str(EL.maxiter)])
		
        if O.EL.rnormal(i,1) < EL.tolr % NORM BELOW TOLERANCE
            O.EL.exitFlag(i,1) = 1;
            break
        else % NORM NOT BELOW TOLERANCE, REQUIRES FOLD ANGLE CORRECTIONS

            dtet = - inv(drdt)*r; % FOLD ANGLE CORRECTIONS
            O.EL.dtetnormal(i,1) = norm(dtet,2)/size(dtet,1);
            
            % SCALING IF FOLD ANGLE INCRMENT IS TOO LARGE
            if max(dtet) > EL.maxdtet
                disp('Correction too large')
                disp(['Max current change: ',num2str(max(dtet)), ', max change allowed: ',num2str(EL.maxdtet)])
                dtet = dtet*EL.maxdtet/max(dtet);
            end

            O.tet(1:end,i) = O.tet(1:end,i) + dtet; % UPDATING FOLD ANGLES
            
            if O.EL.dtetnormal(i,1) < EL.totdtet;
                O.EL.exitFlag(i,1) = 2;
                break
            end
            
        end
    end % END OF ITERATIONS LOOP
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OBTAINING LOCATIONS OF APPLIED LOADS AND
    % DISPLACEMENT OF POINTS OF INTEREST
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j = 1:size(fo,2)
        MF = MAP_FACE(I, O.tet(1:end,i), I.PFolds(EL.FACE0(1,j),1:end)'); % OBTAINING MAPPING 
        temp = MF * [EL.X0(1:3,j) ; 1]; 
        O.EL.xf( (3*i-2):3*i, j) = temp(1:3,1); % LOCATION
    end
    
    for j = 1:size(EL.FACEUP,2)
        MF = MAP_FACE(I, O.tet(1:end,i), I.PFolds(EL.FACEUP(1,j),1:end)'); % OBTAINING MAPPING
        temp = MF * [EL.XUP(1:3,j) ; 1];
        temp = temp(1:3,1);
        O.EL.UP( (3*j-2):3*j, i) = temp - EL.XUP(1:3,j); % DISPLACEMENT
        
    end
    
end % END OF INCREMENTS LOOP