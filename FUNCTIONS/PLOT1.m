% CODE FOR SIMULATION OF ORIGAMI SHEETS WITH SMOOTH ELASTIC FOLDS
% AUTHOR: Edwin A. Peraza Hernandez
% 05-16-2016
% PAPER DOI: http://dx.doi.org/10.1016/j.cad.2016.05.010

function PLOT1(I)
% INITIAL PLOTS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fcolors = [0.85 0.85 0.85]; 
fcolors = [0.21,1,0.21];

wnt_textP = 1; % 1: text in plots (for faces)
wnt_textF = 1; % 1: text in plots (for folds)
wnt_textV = 1; % 1: text in plots (for vertices)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING VERTICES AND FOLD LINES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
title('Vertices and fold centerlines')
box on
hold on

for i = 1:I.N_F
    plot([I.v(I.CM(i,1),1), I.v(I.CM(i,2),1)],...
        [I.v(I.CM(i,1),2), I.v(I.CM(i,2),2)],...
        'k','LineWidth',2)
    if wnt_textF == 1
        text(0.5*(I.v(I.CM(i,1),1)+I.v(I.CM(i,2),1)) + 0.02*max(I.v(1:end,1)),...
            0.5*(I.v(I.CM(i,1),2)+I.v(I.CM(i,2),2)) + 0.02*max(I.v(1:end,1)), ...
            num2str(i), 'FontSize',15, 'Color', 'k')
    end  
end

for i = 1:I.N_V
    plot(I.v(i,1),I.v(i,2),'o','MarkerEdgeColor','none',...
        'MarkerFaceColor','r','MarkerSize',10)
    if wnt_textV == 1
        text(I.v(i,1) + 0.02*max(I.v(1:end,1)),I.v(i,2) + 0.02*max(I.v(1:end,1)), ...
            num2str(i), 'FontSize',15, 'Color', 'r')
    end
end

xlabel('X_1')
ylabel('X_2')
set(gca,'FontName','Times New Roman','fontsize', 20,'linewidth',1.15)
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca,'ticklength',3*get(gca,'ticklength'))
axis equal


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING VERTICES AND FOLD VECTORS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure()
% box on
% hold on
% 
% for i = 1:I.N_V
%     plot(I.v(i,1),I.v(i,2),'o','MarkerEdgeColor','none',...
%         'MarkerFaceColor','r','MarkerSize',10)
%     if wnt_textV == 1    
%         text(I.v(i,1) + 0.04*max(I.v(1:end,1)),I.v(i,2) + 0.04*max(I.v(1:end,1)), ...
%             num2str(i), 'FontSize',15, 'Color', 'r')
%     end
% end
% 
% for i = 1:I.N_F
%     quiver(I.v(I.CM(i,1),1), I.v(I.CM(i,1),2),...
%         I.v(I.CM(i,2),1)-I.v(I.CM(i,1),1), I.v(I.CM(i,2),2)-I.v(I.CM(i,1),2),0,...
%         'k','LineWidth',2,'MaxHeadSize',0.3)
%     
%     if wnt_textF == 1
%         text(0.5*(I.v(I.CM(i,1),1)+I.v(I.CM(i,2),1)) + 0.04*max(I.v(1:end,1)),...
%             0.5*(I.v(I.CM(i,1),2)+I.v(I.CM(i,2),2)) + 0.04*max(I.v(1:end,1)), ...
%             num2str(i), 'FontSize',15, 'Color', 'k')
%     end   
% end
% 
% xlabel('X_1')
% ylabel('X_2')
% set(gca,'FontName','Times New Roman','fontsize', 20,'linewidth',1.15)
% set(gca,'XMinorTick','on','YMinorTick','on')
% set(gca,'ticklength',3*get(gca,'ticklength'))
% axis equal


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING SURFACES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
title('Reference configuration')
box on
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING FOLDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:I.N_F
    px = [I.p1(i,1)  I.p2(i,1)  I.p3(i,1)  I.p4(i,1)]';
    py = [I.p1(i,2)  I.p2(i,2)  I.p3(i,2)  I.p4(i,2)]';
    fill(px, py, fcolors)
    if wnt_textF == 1    
        text(mean(px),mean(py), ...
            num2str(i),'HorizontalAlignment','center','FontSize',15,'Color', 'b')
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING FACES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:I.N_P
    px = zeros(I.Pnpoints(i,1),1);
    py = zeros(I.Pnpoints(i,1),1);
    for j = 1:I.Pnpoints(i,1)
        px(j,1) = I.Pcorners(i,2*j-1);
        py(j,1) = I.Pcorners(i,2*j);
    end
    fill(px, py, Fcolors)
    if wnt_textP == 1
        text(mean(px),mean(py), ...
            num2str(i),'HorizontalAlignment','center','FontSize',15,'Color', 'k')
    end   
end
xlabel('X_1')
ylabel('X_2')
set(gca,'FontName','Times New Roman','fontsize', 20,'linewidth',1.15)
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca,'ticklength',3*get(gca,'ticklength'))
axis equal



drawnow 