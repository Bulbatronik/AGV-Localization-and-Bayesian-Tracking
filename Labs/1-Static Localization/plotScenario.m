function plotScenario( parameters , AP , UE )

nAP = parameters.numberOfAP;
xmin = parameters.xmin;
ymin = parameters.ymin;
xmax = parameters.xmax;
ymax = parameters.ymax;

fig = figure(); hold on
fig.WindowState = 'maximized';
% Localization area
patch( [xmin xmax xmax xmin] , [ymin ymin ymax ymax], [54 114 91]./255 , 'FaceAlpha', .1)

% AP
plot( AP(:,1) , AP(:,2) , '^' , 'MarkerSize', 10 , 'MarkerEdgeColor' , [ 0.64 , 0.08 , 0.18 ] , 'MarkerFaceColor' , [ 0.64 , 0.08 , 0.18 ] )
legend('AP','location','northwest')
axis equal

for a=1:nAP
   text( AP(a,1)+2 , AP(a,2) , sprintf('AP %d ', a) )
end

% UE
plot( UE(1) , UE(2) , 'o' , 'MarkerSize' , 10 , 'MarkerEdgeColor' , [0.30 , 0.75 , 0.93 ] , 'MarkerFaceColor' , [0.30 , 0.75 , 0.93] )

grid on
% grid minor
box on
axis equal
legend( 'Localization Area' , 'AP' , 'UE' , 'location' , 'best' )
xticks( [xmin:20:xmax] )  , yticks( [xmin:20:xmax] )
xlim( [xmin-15 xmax+15] ) , ylim( [ymin-15 ymax+15] )
xlabel( 'meter' , 'FontSize' , 26 ) , ylabel( 'meter' , 'FontSize' , 26 )



end