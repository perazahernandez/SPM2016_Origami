% CODE FOR SIMULATION OF ORIGAMI SHEETS WITH SMOOTH ELASTIC FOLDS
% AUTHOR: Edwin A. Peraza Hernandez
% 05-16-2016
% PAPER DOI: http://dx.doi.org/10.1016/j.cad.2016.05.010

function PLOTCONF(I,O)
% PLOTTING CURRENT CONFIGURATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Fcolors = [0.85 0.85 0.85]; 
fcolors = [0.21,1,0.21];
Load_Scale = 0.05; % SCALE FORCE ARROW TO VISUALIZE LOADS

if I.saveGIF == 1 % IF A .gif IMAGE OF THE OUTPUT FRAMES SHOULD BE WRITTEN
    filename = 'output.gif';
end


for i = 1:size(I.plotInc,1) % REPEAT FOR ALL REQUESTED INCREMENTS
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OBTAINING FOLD ANGLES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tet = O.tet(1:end, I.plotInc(i)+1);
    
    if I.saveGIF == 1 % IF A .gif IMAGE OF THE OUTPUT FRAMES SHOULD BE WRITTEN
        figure('units','normalized','position',[.1 .1 1 1],...
            'Color',[1 1 1])
    else
        figure('units','normalized','position',[.15 .15 0.5 0.5])
    end
    
    hold on  
    title(['Increment: ', num2str(I.plotInc(i))])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOTTING FACES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j = 1:I.N_P
        
        MF = MAP_FACE(I, tet, I.PFolds(j,1:end)'); % OBTAINING MAPPING FOR FACE
        
        px = zeros(I.Pnpoints(j),1);
        py = zeros(I.Pnpoints(j),1);
        pz = zeros(I.Pnpoints(j),1);
        for k = 1:I.Pnpoints(j)
            temp = MF * [I.Pcorners(j,2*k-1); I.Pcorners(j,2*k); 0; 1];
            px(k) = temp(1);
            py(k) = temp(2);
            pz(k) = temp(3);
        end
        
        fill3(px,py,pz,(Fcolors))
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOTTING FOLDS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j = 1:I.N_F
        
        MF = MAP_FACE(I, tet, I.PFolds( I.PconF(j,1),1:end)'); % OBTAINING MAPPING FOR ADJACENT FACE
        
        MFE = MAP_FOLDEX(I, tet, j); % FURTHER TRANSFORMATIONS
        
        MF = MF*MFE;
        
        FSP = FOLD_SURF_POINTS(I, tet, j, MF); % OBTAINING SURFACE POINTS 
        px = FSP(1:end,1:end,1);
        py = FSP(1:end,1:end,2);
        pz = FSP(1:end,1:end,3);
        
        surf(px(1:end,1:end),...
            py(1:end,1:end),...
            pz(1:end,1:end),...
            'MeshStyle','column',...
            'FaceColor',(fcolors),...
            'FaceLighting','gouraud')
    end
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOTTING APPLIED LOADS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j = 1:O.EL.N_f
        quiver3(O.EL.xf(3*(I.plotInc(i)+1)-2,j) , ...- Load_Scale*O.EL.fo(3*(I.plotInc(i)+1)-2,j), ...
            O.EL.xf(3*(I.plotInc(i)+1)-1,j) , ...- Load_Scale*O.EL.fo(3*(I.plotInc(i)+1)-1,j), ...
            O.EL.xf(3*(I.plotInc(i)+1),j) , ...- Load_Scale*O.EL.fo(3*(I.plotInc(i)+1),j), ...
            Load_Scale*O.EL.fo(3*(I.plotInc(i)+1)-2,j),...
            Load_Scale*O.EL.fo(3*(I.plotInc(i)+1)-1,j),...
            Load_Scale*O.EL.fo(3*(I.plotInc(i)+1),j),...
            0,...
            'r','LineWidth',2,'MaxHeadSize',0.3);
    end
    
    %view(40, 16) % VIEW ANGLE

	% LIGHTNING
    %light('Position',[-1 -1 1])
    %light('Position',[0 1 1])
	
    xlabel('X_1')
    ylabel('X_2')
    zlabel('X_3')
    set(gca,'FontName','Times New Roman','fontsize', 20,'linewidth',1.15)
    set(gca,'XMinorTick','on','YMinorTick','on')
    set(gca,'ticklength',3*get(gca,'ticklength'))
    axis equal
    %camproj('perspective')
    
    
    
    if I.saveGIF == 1 % IF A .gif IMAGE OF THE OUTPUT FRAMES SHOULD BE WRITTEN
        title ''
        axis([-0.05 0.85 0 0.4  -0.2 0.05])
        axis off
		
        zoom(1.6)  % ZOOM FACTOR
		
        drawnow
        frame = getframe;
        
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if i == 1
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf,...
                'DelayTime',0.125);
        else
            if i == size(I.plotInc,1)
                imwrite(imind,cm,filename,'gif','WriteMode','append',...
                    'DelayTime',1.5);
            else
                imwrite(imind,cm,filename,'gif','WriteMode','append',...
                    'DelayTime',0.125);
            end
        end
    end

end