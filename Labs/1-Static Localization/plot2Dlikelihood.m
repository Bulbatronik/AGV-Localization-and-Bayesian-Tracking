function plot2Dlikelihood(parameters, AP , UE , x , y , likelihood , TYPE)

fig = figure();
fig.WindowState = 'maximized';
for a = 1:parameters.numberOfAP
    
    if parameters.numberOfAP>2
        subplot(round(parameters.numberOfAP/2),round(parameters.numberOfAP/2),a)
    else
        subplot(2,1,a)
    end
    imagesc( x , y ,   (   ( squeeze(  likelihood(a,:,:) )' )  ) )   ;  hold on

    set(gca,'YDir','normal')

    plot( AP(:,1) , AP(:,2) , '^','MarkerSize',10,'MarkerEdgeColor',[0.64,0.08,0.18],'MarkerFaceColor',[0.64,0.08,0.18] );    hold on
    plot( AP(a,1) , AP(a,2) , '^','MarkerSize',10,'MarkerEdgeColor',[102,254,0]./255,'MarkerFaceColor',[102,254,0]./255 )
    plot( UE(:,1) , UE(:,2) , 'o','MarkerSize',10,'MarkerEdgeColor',[0.30,0.75,0.93],'MarkerFaceColor',[0.30,0.75,0.93] )


    colorbar;
    axis equal
    xlim([parameters.xmin parameters.xmax]) , ylim([parameters.ymin parameters.ymax])
    xlabel('[m]','FontSize',26), ylabel('[m]','FontSize',26);
   

    switch TYPE
        case 'TOA'
            title(['Likelihood ',num2str(TYPE),', ${AP}$ = ',num2str(a),' , $\sigma $ = ',num2str(parameters.sigmaTOA),' m '],'Interpreter','Latex')
        case 'RSS'
            title(['Likelihood ',num2str(TYPE),', $N_{AP}$ = ',num2str(a),' , $\sigma $ = ',num2str(parameters.sigmaRSS),' m '],'Interpreter','Latex')
       case 'AOA'
            title(['Likelihood ',num2str(TYPE),', ${AP}$ = ',num2str(a),' , $\sigma $ = ',num2str(rad2deg(parameters.sigmaAOA)),' deg '],'Interpreter','Latex')
        case 'TDOA'
            title(['Likelihood ',num2str(TYPE),', ${AP}$ = 1-',num2str(a),' ,   $\sigma $ = ',num2str(parameters.sigmaTDOA),' dB '],'Interpreter','Latex')
    end
     pause
end