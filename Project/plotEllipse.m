function plotEllipse( parameters , AP , ellipsePoints , UE , TYPE )
x_ell = ellipsePoints( 1 , : );
y_ell = ellipsePoints( 2 , : );

%fig = figure(); 
hold on
fig.WindowState = 'maximized';
% patch( [parameters.xmin parameters.xmax parameters.xmax parameters.xmin parameters.xmin],[parameters.ymin parameters.ymin parameters.ymax parameters.ymax parameters.ymin], [54 114 91]./255 , 'FaceAlpha', .1)
plot( AP(:,1) , AP(:,2) , '^','MarkerSize',10,'MarkerEdgeColor',[0.64,0.08,0.18],'MarkerFaceColor' , [0.64,0.08,0.18])
fill( x_ell , y_ell , [.9 .95 1],'edgecolor',[0, 0.4470, 0.7410],'linewidth',2);alpha(.5)
plot( UE(1) , UE(2) , 'o','MarkerSize',10,'MarkerEdgeColor',[0.30,0.75,0.93],'MarkerFaceColor' , [0.30,0.75,0.93])
plot( [AP(:,1) , repmat(UE(1),size(AP,1),1)]' , [AP(:,2),repmat(UE(2),size(AP,1),1)]','--','LineWidth',1,'color',[220,220,220]./255)
% legend('Localization Area','AP','Ellipse','UE','location','best')
% legend('AP','Ellipse','UE','link','location','best')
xticks( [parameters.xmin:20:parameters.xmax] )  , yticks( [parameters.xmin:20:parameters.xmax] )
xlim( [parameters.xmin-15 parameters.xmax+15] ) , ylim( [parameters.ymin-15 parameters.ymax+15] )
xlabel('meter','FontSize',26) , ylabel('meter','FontSize',26)
grid on
box on
axis equal

% axis equal
switch TYPE
    case 'TOA'
        title([' ',num2str(TYPE),', $N_{AP}$ = ',num2str(parameters.numberOfAP),' , $\sigma $ = ',num2str(parameters.sigmaTOA),' m'],'Interpreter','Latex')
    case 'RSS'
        title([' ',num2str(TYPE),', $N_{AP}$ = ',num2str(parameters.numberOfAP),' , $\sigma $ = ',num2str(parameters.sigmaRSS),' dB'],'Interpreter','Latex')
    case 'AOA'
        title([' ',num2str(TYPE),', $N_{AP}$ = ',num2str(parameters.numberOfAP),' , $\sigma $ = ',num2str(rad2deg(parameters.sigmaAOA)),' deg '],'Interpreter','Latex')
    case 'TDOA'
        title([' ',num2str(TYPE),', $N_{AP}$ = ',num2str(parameters.numberOfAP),' , $\sigma $ = ',num2str(parameters.sigmaTDOA),' m'],'Interpreter','Latex')
end


end