% CODE FOR SIMULATION OF ORIGAMI SHEETS WITH SMOOTH ELASTIC FOLDS
% AUTHOR: Edwin A. Peraza Hernandez
% 05-16-2016
% PAPER DOI: http://dx.doi.org/10.1016/j.cad.2016.05.010

function PLOTANG(I,O)
% PLOTTING REQUESTED FOLD ANGLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(I.plotAng,1) > 0
	CM = hsv(size(I.plotAng,1));
    figure()
    box on
    hold on
    
    incvec = linspace(0,size(O.tet,2)-1,size(O.tet,2)); % ARRAY WITH INCREMENT NUMBERS

    for i = 1:size(I.plotAng,1)
        plot(incvec, 180*O.tet(I.plotAng(i,1),1:end)/pi,...
            'LineWidth',2,'color',CM(i,:))
        legendInfo{i} = ['$\hat{\theta}_' num2str(I.plotAng(i,1)) '$'];
    end
    
    legend(legendInfo,'Interpreter','latex')
    xlabel('Increment')
    ylabel('Angle [deg]')
    set(gca,'FontName','Times New Roman','fontsize', 20,'linewidth',1.15)
    set(gca,'XMinorTick','on','YMinorTick','on')
    set(gca,'ticklength',3*get(gca,'ticklength'))
end