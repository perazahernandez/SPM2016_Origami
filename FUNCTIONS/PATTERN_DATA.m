% CODE FOR SIMULATION OF ORIGAMI SHEETS WITH SMOOTH ELASTIC FOLDS
% AUTHOR: Edwin A. Peraza Hernandez
% 05-16-2016
% PAPER DOI: http://dx.doi.org/10.1016/j.cad.2016.05.010

function I = PATTERN_DATA(I,EL)
% FOLD PATTERN DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R3 = @(t) [cos(t),  -sin(t),    0;
    sin(t),     cos(t),     0;
    0,      0,      1]; % MATRIX FOR ROTATION ABOUT 3-AXIS


I.N_F = size(I.CM,1); % NUMBER OF FOLDS
I.N_V = size(I.v,1); % NUMBER OF VERTICES
I.N_I = size(I.CF,1); % NUMBER OF FOLD INTERSECTIONS
I.N_P = size(I.P_FoldVert,1); % NUMBER OF FACES


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NUMBER OF FOLDS ASSOCIATED WITH EACH INTERSECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I.n_j = zeros(I.N_I,1); % NUMBER OF FOLDS ASSOCIATED WITH EACH INTERSECTION
for i = 1:I.N_I
    for j = 1:size(I.CF,2)
        if I.CF(i,j) == 0
            break
        end
        I.n_j(i) = j;
    end
end

% INITIALIZING ARRAYS
I.m = zeros(I.N_F,2); % FOLD VECTORS
I.halp = zeros(I.N_F,1); % COUNTERCLOCKWISE ANGLE WITH RESPECT TO 1-AXIS

% FOLD CORNER POINTS
I.p1 = zeros(I.N_F,2);
I.p2 = zeros(I.N_F,2);
I.p3 = zeros(I.N_F,2);
I.p4 = zeros(I.N_F,2);

for i = 1:I.N_F
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FOLD VECTORS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    I.m(i,1) = I.v(I.CM(i,2),1) - I.v(I.CM(i,1),1);
    I.m(i,2) = I.v(I.CM(i,2),2) - I.v(I.CM(i,1),2);

	
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ANGLE BETWEEN FOLD VECTORS AND e1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    if I.m(i,2) >= 0
        I.halp(i,1) = acos(I.m(i,1)/norm(I.m(i,1:2),2));
    else
        I.halp(i,1) = 2*pi - acos(I.m(i,1)/norm(I.m(i,1:2),2));
    end
 
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FOLD CORNER POITNS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    I.p1(i,1:2) = I.v(I.CM(i,1),1:2) +(1/norm(I.m(i,1:2),2))*(...
        - I.w(i,1)/2 * [-I.m(i,2), I.m(i,1)] ...
        + I.d(i,1) * I.m(i,1:2)  );
    I.p2(i,1:2) = I.v(I.CM(i,1),1:2) + I.m(i,1:2) +(1/norm(I.m(i,1:2),2))*(...
        - I.w(i,1)/2 * [-I.m(i,2), I.m(i,1)] ...
        - I.d(i,2) * I.m(i,1:2)  );    
    I.p3(i,1:2) = I.v(I.CM(i,1),1:2) + I.m(i,1:2) +(1/norm(I.m(i,1:2),2))*(...
        + I.w(i,1)/2 * [-I.m(i,2), I.m(i,1)] ...
        - I.d(i,2) * I.m(i,1:2)  ); 
    I.p4(i,1:2) = I.v(I.CM(i,1),1:2) +(1/norm(I.m(i,1:2),2))*(...
        + I.w(i,1)/2 * [-I.m(i,2), I.m(i,1)] ...
        + I.d(i,1) * I.m(i,1:2)  );    
end


I.mI = zeros(max(I.n_j), 2*I.N_I);

I.halpI = zeros(I.N_I,max(I.n_j));

for i = 1:I.N_I
    for j = 1:I.n_j(i,1)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % FOLD VECTORS PER INTERSECTION
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        I.mI(j,(2*i-1):(2*i)) = sign(I.CF(i,j)) * I.m(abs(I.CF(i,j)),1:2);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ANGLE BETWEEN FOLD INTERSECTION VECTORS AND e1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if I.mI(j,(2*i)) >= 0
            I.halpI(i,j) = acos(I.mI(j,(2*i-1))/norm(I.mI(j,(2*i-1):(2*i)),2));
        else
            I.halpI(i,j) = 2*pi - acos(I.mI(j,(2*i-1))/norm(I.mI(j,(2*i-1):(2*i)),2));
        end
        
    end
end


I.alp = zeros(I.N_I, max(I.n_j));

for i = 1:I.N_I
    for j = 1:I.n_j(i,1)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % FACE CORNER ANGLES
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if j == I.n_j(i,1)
            I.alp(i,j) = acos( dot(I.mI(j,(2*i-1):(2*i)), I.mI(1,(2*i-1):(2*i)))/...
                (norm(I.mI(j,(2*i-1):(2*i)),2)*norm(I.mI(1,(2*i-1):(2*i)),2)) );            
        else
            I.alp(i,j) = acos( dot(I.mI(j,(2*i-1):(2*i)), I.mI(j+1,(2*i-1):(2*i)))/...
                (norm(I.mI(j,(2*i-1):(2*i)),2)*norm(I.mI(j+1,(2*i-1):(2*i)),2)) );
        end
    end
end



% MAPPING FOLD WIDTHS ASSOCIATED WITH EACH FOLD INTERSECTION
I.tw = zeros(I.N_I,max(I.n_j));

% MAPPING CORNERS POINTS ASSOCIATED WITH EACH FOLD INTERSECTION
I.BL = zeros(max(I.n_j),  2*I.N_I);
I.BR = zeros(max(I.n_j),  2*I.N_I);

for j = 1:I.N_I
    for k = 1:I.n_j(j,1)
        I.tw(j,k) = I.w(abs(I.CF(j,k)));
        if I.CF(j,k) > 0
            I.BL(k, (2*j-1):2*j) = I.p1( abs(I.CF(j,k)), 1:2);
            I.BR(k, (2*j-1):2*j) = I.p4( abs(I.CF(j,k)), 1:2);
        else
            I.BL(k, (2*j-1):2*j) = I.p4( abs(I.CF(j,k)), 1:2);
            I.BR(k, (2*j-1):2*j) = I.p1( abs(I.CF(j,k)), 1:2);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ljk VECTORS (SEE PAPERS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I.lvecs = zeros(max(I.n_j),  2*I.N_I);
for j = 1:I.N_I
    for k = 1:I.n_j(j,1)
        
        if k == I.n_j(j,1)
            temp1 = I.BL(1, (2*j-1):2*j) - I.BR(k, (2*j-1):2*j);
        else
            temp1 = I.BL(k+1, (2*j-1):2*j) - I.BR(k, (2*j-1):2*j);
        end
        
        temp2 = [temp1';0];
        temp3 = ( R3(-I.halpI(j,k)) * temp2 )';
        I.lvecs(k, (2*j-1):2*j) = temp3(1:2);
        
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FACE CORNER POINTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I.Pcorners = zeros(I.N_P, 2*2*size(I.P_FoldVert,2));
I.Pnpoints = zeros(I.N_P,1);

for i = 1:I.N_P 
    cnt = 1; % Counter 
    for j = 1:size(I.P_FoldVert,2)
        if I.P_FoldVert(i,j) == 0
            break
        end
        
        if I.P_FoldVert(i,j) > I.N_F % Connected to isolated vertex
            I.Pcorners(i, (2*cnt-1):2*cnt) = I.v(I.P_FoldVert(i,j) - I.N_F,1:2);
            cnt = cnt + 1;
        else
            if I.P_FoldVert(i,j) > 0 % Connected to fold in positive side
                I.Pcorners(i, (2*cnt-1):2*cnt) = I.p4( I.P_FoldVert(i,j), 1:2);
                cnt = cnt + 1;
                I.Pcorners(i, (2*cnt-1):2*cnt) = I.p3( I.P_FoldVert(i,j), 1:2);
                cnt = cnt + 1;                
            elseif I.P_FoldVert(i,j) < 0 % Connected to fold in negative side
                I.Pcorners(i, (2*cnt-1):2*cnt) = I.p2( abs(I.P_FoldVert(i,j)), 1:2);
                cnt = cnt + 1;
                I.Pcorners(i, (2*cnt-1):2*cnt) = I.p1( abs(I.P_FoldVert(i,j)), 1:2);
                cnt = cnt + 1;                  
            end
        end
    end
    I.Pnpoints(i,1) = cnt - 1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATING AREA OF FACES AND CENTROID POINTS 
% (USED WHEN GRAVITY IS APPLIED)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I.Parea = zeros(I.N_P,1);
I.Pcentroid = zeros(I.N_P,2);

for i = 1:I.N_P % FOR each face
    temp = 0;
    tempCx = 0;
    tempCy = 0;
    for j = 1:I.Pnpoints(i,1) % FOR each face vertex
        if j < I.Pnpoints(i,1)
            temp = temp + 0.5*...
                ( I.Pcorners(i, 2*j-1)*I.Pcorners(i, 2*(j+1)) - ...
                I.Pcorners(i, 2*j)*I.Pcorners(i, 2*(j+1)-1) );
            tempCx = tempCx + (I.Pcorners(i, 2*j-1) + I.Pcorners(i, 2*(j+1)-1))*...
                ( I.Pcorners(i, 2*j-1)*I.Pcorners(i, 2*(j+1)) - ...
                I.Pcorners(i, 2*j)*I.Pcorners(i, 2*(j+1)-1) );  
            tempCy = tempCy + (I.Pcorners(i, 2*j) + I.Pcorners(i, 2*(j+1)))*...
                ( I.Pcorners(i, 2*j-1)*I.Pcorners(i, 2*(j+1)) - ...
                I.Pcorners(i, 2*j)*I.Pcorners(i, 2*(j+1)-1) );            
        else
            temp = temp + 0.5*...
                ( I.Pcorners(i, 2*j-1)*I.Pcorners(i, 2*(1)) - ...
                I.Pcorners(i, 2*j)*I.Pcorners(i, 2*(1)-1) );  
            tempCx = tempCx + (I.Pcorners(i, 2*j-1) + I.Pcorners(i, 2*(1)-1))*...
                ( I.Pcorners(i, 2*j-1)*I.Pcorners(i, 2*(1)) - ...
                I.Pcorners(i, 2*j)*I.Pcorners(i, 2*(1)-1) );             
            tempCy = tempCy + (I.Pcorners(i, 2*j) + I.Pcorners(i, 2*(1)))*...
                ( I.Pcorners(i, 2*j-1)*I.Pcorners(i, 2*(1)) - ...
                I.Pcorners(i, 2*j)*I.Pcorners(i, 2*(1)-1) );               
        end
    end
    
    I.Parea(i) = temp;
    I.Pcentroid(i,1) = tempCx/(6*I.Parea(i));
    I.Pcentroid(i,2) = tempCy/(6*I.Parea(i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATING AREA OF FOLDS AND CENTROID POINTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I.Farea = zeros(I.N_F,1);
I.Fcentroid = zeros(I.N_F,2);

for i = 1:I.N_F % FOR each fold
    
    I.Farea(i) = ( ( (I.p1(i,1)-I.p2(i,1))^2 + (I.p1(i,2)-I.p2(i,2))^2)^(1/2)...
        *( (I.p1(i,1)-I.p4(i,1))^2 + (I.p1(i,2)-I.p4(i,2))^2)^(1/2));
    I.Fcentroid(i,1) = 0.25*(I.p1(i,1)+I.p2(i,1)+I.p3(i,1)+I.p4(i,1));
    I.Fcentroid(i,2) = 0.25*(I.p1(i,2)+I.p2(i,2)+I.p3(i,2)+I.p4(i,2));
end




%------------------------------------------------------------------
% CHECKING IF FOLDS OVERLAP
%------------------------------------------------------------------
% I.overlap = 0;
% I.overlap_Pairs = [0,0];
% I.overlap_stimated_l = [0];
% 
% for i = 1:I.N_F-1
%     
%     projv1 = [ - (I.p2(i,2)-I.p1(i,2)); I.p2(i,1)-I.p1(i,1) ];
%     projv2 = [ - (I.p3(i,2)-I.p2(i,2)); I.p3(i,1)-I.p2(i,1) ];
%     projv3 = [ - (I.p4(i,2)-I.p3(i,2)); I.p4(i,1)-I.p3(i,1) ];
%     projv4 = [ - (I.p1(i,2)-I.p4(i,2)); I.p1(i,1)-I.p4(i,1) ];
%     % normalizing projection vectors
%     projv1 = projv1/norm(projv1,2);     projv2 = projv2/norm(projv2,2);
%     projv3 = projv3/norm(projv3,2);     projv4 = projv4/norm(projv4,2);
%     for j = i+1:I.N_F
%         
%         projv1o = [ - (I.p2(j,2)-I.p1(j,2)); I.p2(j,1)-I.p1(j,1) ];
%         projv2o = [ - (I.p3(j,2)-I.p2(j,2)); I.p3(j,1)-I.p2(j,1) ];
%         projv3o = [ - (I.p4(j,2)-I.p3(j,2)); I.p4(j,1)-I.p3(j,1) ];
%         projv4o = [ - (I.p1(j,2)-I.p4(j,2)); I.p1(j,1)-I.p4(j,1) ];
%         % normalizing projection vectors
%         projv1o = projv1o/norm(projv1o,2);     projv2o = projv2o/norm(projv2o,2);
%         projv3o = projv3o/norm(projv3o,2);     projv4o = projv4o/norm(projv4o,2);
%         
%         
%         temp = 1;
%         templen = 0;
%         
%         %------------------------------------------------------------------
%         % SEARCHING FOR SEPARATION AXIS IN FOLD i
%         %------------------------------------------------------------------
%         if temp ~= 0
%             tempv1 = [I.p1(j,1) - I.p1(i,1); I.p1(j,2) - I.p1(i,2)];
%             chck1 = dot(projv1,tempv1);
%             tempv2 = [I.p2(j,1) - I.p1(i,1); I.p2(j,2) - I.p1(i,2)];
%             chck2 = dot(projv1,tempv2);
%             tempv3 = [I.p3(j,1) - I.p1(i,1); I.p3(j,2) - I.p1(i,2)];
%             chck3 = dot(projv1,tempv3);
%             tempv4 = [I.p4(j,1) - I.p1(i,1); I.p4(j,2) - I.p1(i,2)];
%             chck4 = dot(projv1,tempv4);
%             if chck1 <= 0 && chck2 <= 0 && chck3 <= 0 && chck4 <= 0
%                 temp = 0;
% %             else % OBTAINING STIMATED OVERLAP LENGTH
% %                 if chck1 > 0
% %                     templen = min([chck1, templen]);
% %                 end
% %                 if chck2 > 0
% %                     templen = min([chck2, templen]);
% %                 end      
% %                 if chck3 > 0
% %                     templen = min([chck3, templen]);
% %                 end
% %                 if chck4 > 0
% %                     templen = min([chck4, templen]);
% %                 end                    
%             end
%         end
%         
%         if temp ~= 0
%             tempv1 = [I.p1(j,1) - I.p2(i,1); I.p1(j,2) - I.p2(i,2)];
%             chck1 = dot(projv2,tempv1);
%             tempv2 = [I.p2(j,1) - I.p2(i,1); I.p2(j,2) - I.p2(i,2)];
%             chck2 = dot(projv2,tempv2);
%             tempv3 = [I.p3(j,1) - I.p2(i,1); I.p3(j,2) - I.p2(i,2)];
%             chck3 = dot(projv2,tempv3);
%             tempv4 = [I.p4(j,1) - I.p2(i,1); I.p4(j,2) - I.p2(i,2)];
%             chck4 = dot(projv2,tempv4);
%             if chck1 <= 0 && chck2 <= 0 && chck3 <= 0 && chck4 <= 0
%                 temp = 0;
% %             else % OBTAINING STIMATED OVERLAP LENGTH
% %                 if chck1 > 0
% %                     templen = min([chck1, templen]);
% %                 end
% %                 if chck2 > 0
% %                     templen = min([chck2, templen]);
% %                 end      
% %                 if chck3 > 0
% %                     templen = min([chck3, templen]);
% %                 end
% %                 if chck4 > 0
% %                     templen = min([chck4, templen]);
% %                 end                    
%             end
%         end
%         
%         
%         if temp ~= 0
%             tempv1 = [I.p1(j,1) - I.p3(i,1); I.p1(j,2) - I.p3(i,2)];
%             chck1 = dot(projv3,tempv1);
%             tempv2 = [I.p2(j,1) - I.p3(i,1); I.p2(j,2) - I.p3(i,2)];
%             chck2 = dot(projv3,tempv2);
%             tempv3 = [I.p3(j,1) - I.p3(i,1); I.p3(j,2) - I.p3(i,2)];
%             chck3 = dot(projv3,tempv3);
%             tempv4 = [I.p4(j,1) - I.p3(i,1); I.p4(j,2) - I.p3(i,2)];
%             chck4 = dot(projv3,tempv4);
%             if chck1 <= 0 && chck2 <= 0 && chck3 <= 0 && chck4 <= 0
%                 temp = 0;
% %             else % OBTAINING STIMATED OVERLAP LENGTH
% %                 if chck1 > 0
% %                     templen = min([chck1, templen]);
% %                 end
% %                 if chck2 > 0
% %                     templen = min([chck2, templen]);
% %                 end      
% %                 if chck3 > 0
% %                     templen = min([chck3, templen]);
% %                 end
% %                 if chck4 > 0
% %                     templen = min([chck4, templen]);
% %                 end                    
%             end
%         end
%         
%         
%         if temp ~= 0
%             tempv1 = [I.p1(j,1) - I.p4(i,1); I.p1(j,2) - I.p4(i,2)];
%             chck1 = dot(projv4,tempv1);
%             tempv2 = [I.p2(j,1) - I.p4(i,1); I.p2(j,2) - I.p4(i,2)];
%             chck2 = dot(projv4,tempv2);
%             tempv3 = [I.p3(j,1) - I.p4(i,1); I.p3(j,2) - I.p4(i,2)];
%             chck3 = dot(projv4,tempv3);
%             tempv4 = [I.p4(j,1) - I.p4(i,1); I.p4(j,2) - I.p4(i,2)];
%             chck4 = dot(projv4,tempv4);
%             if chck1 <= 0 && chck2 <= 0 && chck3 <= 0 && chck4 <= 0
%                 temp = 0;
% %             else % OBTAINING STIMATED OVERLAP LENGTH
% %                 if chck1 > 0
% %                     templen = min([chck1, templen]);
% %                 end
% %                 if chck2 > 0
% %                     templen = min([chck2, templen]);
% %                 end      
% %                 if chck3 > 0
% %                     templen = min([chck3, templen]);
% %                 end
% %                 if chck4 > 0
% %                     templen = min([chck4, templen]);
% %                 end                    
%             end
%         end
%         
%         
% 
%         %------------------------------------------------------------------
%         % SEARCHING FOR SEPARATION AXIS IN FOLD j
%         %------------------------------------------------------------------
%         if temp ~= 0
%             tempv1 = [I.p1(i,1) - I.p1(j,1); I.p1(i,2) - I.p1(j,2)];
%             chck1 = dot(projv1o,tempv1);
%             tempv2 = [I.p2(i,1) - I.p1(j,1); I.p2(i,2) - I.p1(j,2)];
%             chck2 = dot(projv1o,tempv2);
%             tempv3 = [I.p3(i,1) - I.p1(j,1); I.p3(i,2) - I.p1(j,2)];
%             chck3 = dot(projv1o,tempv3);
%             tempv4 = [I.p4(i,1) - I.p1(j,1); I.p4(i,2) - I.p1(j,2)];
%             chck4 = dot(projv1o,tempv4);
%             if chck1 <= 0 && chck2 <= 0 && chck3 <= 0 && chck4 <= 0
%                 temp = 0;
% %             else % OBTAINING STIMATED OVERLAP LENGTH
% %                 if chck1 > 0
% %                     templen = min([chck1, templen]);
% %                 end
% %                 if chck2 > 0
% %                     templen = min([chck2, templen]);
% %                 end      
% %                 if chck3 > 0
% %                     templen = min([chck3, templen]);
% %                 end
% %                 if chck4 > 0
% %                     templen = min([chck4, templen]);
% %                 end                    
%             end
%         end
%         
%         if temp ~= 0
%             tempv1 = [I.p1(i,1) - I.p2(j,1); I.p1(i,2) - I.p2(j,2)];
%             chck1 = dot(projv2o,tempv1);
%             tempv2 = [I.p2(i,1) - I.p2(j,1); I.p2(i,2) - I.p2(j,2)];
%             chck2 = dot(projv2o,tempv2);
%             tempv3 = [I.p3(i,1) - I.p2(j,1); I.p3(i,2) - I.p2(j,2)];
%             chck3 = dot(projv2o,tempv3);
%             tempv4 = [I.p4(i,1) - I.p2(j,1); I.p4(i,2) - I.p2(j,2)];
%             chck4 = dot(projv2o,tempv4);
%             if chck1 <= 0 && chck2 <= 0 && chck3 <= 0 && chck4 <= 0
%                 temp = 0;
% %             else % OBTAINING STIMATED OVERLAP LENGTH
% %                 if chck1 > 0
% %                     templen = min([chck1, templen]);
% %                 end
% %                 if chck2 > 0
% %                     templen = min([chck2, templen]);
% %                 end      
% %                 if chck3 > 0
% %                     templen = min([chck3, templen]);
% %                 end
% %                 if chck4 > 0
% %                     templen = min([chck4, templen]);
% %                 end                    
%             end
%         end
%         
%         
%         if temp ~= 0
%             tempv1 = [I.p1(i,1) - I.p3(j,1); I.p1(i,2) - I.p3(j,2)];
%             chck1 = dot(projv3o,tempv1);
%             tempv2 = [I.p2(i,1) - I.p3(j,1); I.p2(i,2) - I.p3(j,2)];
%             chck2 = dot(projv3o,tempv2);
%             tempv3 = [I.p3(i,1) - I.p3(j,1); I.p3(i,2) - I.p3(j,2)];
%             chck3 = dot(projv3o,tempv3);
%             tempv4 = [I.p4(i,1) - I.p3(j,1); I.p4(i,2) - I.p3(j,2)];
%             chck4 = dot(projv3o,tempv4);
%             if chck1 <= 0 && chck2 <= 0 && chck3 <= 0 && chck4 <= 0
%                 temp = 0;
% %             else % OBTAINING STIMATED OVERLAP LENGTH
% %                 if chck1 > 0
% %                     templen = min([chck1, templen]);
% %                 end
% %                 if chck2 > 0
% %                     templen = min([chck2, templen]);
% %                 end      
% %                 if chck3 > 0
% %                     templen = min([chck3, templen]);
% %                 end
% %                 if chck4 > 0
% %                     templen = min([chck4, templen]);
% %                 end                    
%             end
%         end
%         
%         
%         if temp ~= 0
%             tempv1 = [I.p1(i,1) - I.p4(j,1); I.p1(i,2) - I.p4(j,2)];
%             chck1 = dot(projv4o,tempv1);
%             tempv2 = [I.p2(i,1) - I.p4(j,1); I.p2(i,2) - I.p4(j,2)];
%             chck2 = dot(projv4o,tempv2);
%             tempv3 = [I.p3(i,1) - I.p4(j,1); I.p3(i,2) - I.p4(j,2)];
%             chck3 = dot(projv4o,tempv3);
%             tempv4 = [I.p4(i,1) - I.p4(j,1); I.p4(i,2) - I.p4(j,2)];
%             chck4 = dot(projv4o,tempv4);
%             if chck1 <= 0 && chck2 <= 0 && chck3 <= 0 && chck4 <= 0
%                 temp = 0;
% %             else % OBTAINING STIMATED OVERLAP LENGTH
% %                 if chck1 > 0
% %                     templen = min([chck1, templen]);
% %                 end
% %                 if chck2 > 0
% %                     templen = min([chck2, templen]);
% %                 end      
% %                 if chck3 > 0
% %                     templen = min([chck3, templen]);
% %                 end
% %                 if chck4 > 0
% %                     templen = min([chck4, templen]);
% %                 end                    
%             end
%         end        
%         
%         
%         
%         if temp == 1
%             I.overlap = 1;
%             I.overlap_Pairs = [I.overlap_Pairs;
%                 i,j];
%         end
%         
%         
%     end
% end
% 
% I.overlap_nPairs = size(I.overlap_Pairs,1)-1;
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LENGTH AND AREA MOMENT OF INERTIA OF FOLDS
% NORMALIZED VECTOR m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I.Li = zeros(I.N_F,1); 
I.Ii = zeros(I.N_F,1); 
I.unitmi = zeros(I.N_F,2); 
for i = 1:I.N_F
    I.Li(i,1) = norm(I.m(i,1:2)) - I.d(i,1) - I.d(i,2); % LENGTH OF FOLDS
    I.Ii(i,1) = (1/12) * I.Li(i,1) * (EL.h(i))^3; % AREA MOMENT OF INTERTIA OF FOLDS
    I.unitmi(i,1:2) = [I.m(i,1),I.m(i,2)]/norm(I.m(i,1:2)); % NORMALIZED VECTOR m
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INDICES OF FOLDS CROSSED BY THE PATHS CONNECTING THE FIXED FACE TO
% EVERY OTHER FACE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I.PFolds = zeros(I.N_P, I.N_F);

fdone = zeros(I.N_P,1);
fdone(I.FixedFace) = 1;
all_done = 0;

Pcur = [I.FixedFace];

countr = 1;

while all_done ~= 1
    Pnew = [];
    
    % CHECKING NEIGHBORS OF EACH FACE IN Pcur
    for i = 1:size(Pcur,1)
        for j = 1:size(I.P_FoldVert,2)
            if abs(I.P_FoldVert(Pcur(i),j)) > 0 && abs(I.P_FoldVert(Pcur(i),j)) <= I.N_F % IF CURRENT FACE HAS A NEIGHBORING FOLD
                % CHECK WHICH OTHER FACE HAS THIS FOLD AS NEIGHBOR
                fcur = I.P_FoldVert(Pcur(i),j);
                for k = 1:I.N_P
                    % CHECK IF FACE PATH IS ALREADY CONTRUCTED, DO NOT CONSIDER IF TRUE
                    if fdone(k) == 1
                    else
                        for l = 1:size(I.P_FoldVert,2)
                            if abs(I.P_FoldVert(k,l)) == abs(fcur) % FACES ARE NEIGHBORS AND SHARE FOLD fcur
                                
                                I.PFolds(k,1:end) = I.PFolds(Pcur(i),1:end);
                                I.PFolds(k,countr) = -fcur;
                                
                                fdone(k) = 1;
                                Pnew = [Pnew; k];
                            end
                        end
                    end
                end
            end
        end
    end

    Pcur = Pnew;
    countr = countr + 1;

    % CHECKING IF ALL FACES ARE ALREADY CONSIDERED 
    all_done = 1;
    if any(fdone == 0) == 1
        all_done = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OBTAINING FACE INDICES THAT ARE NEIGHBORS TO EACH FOLD IN THE p1-p2
% BOUNDARY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I.PconF = zeros(I.N_F,1);
for i = 1:size(I.P_FoldVert,1)
    for j = 1:size(I.P_FoldVert,2)
        if I.P_FoldVert(i,j) < 0
            I.PconF(abs(I.P_FoldVert(i,j))) = i;
        end
    end
end
