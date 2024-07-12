function plotNLSiterations( parameters , AP , UE , uHat , TYPE )

xmin = parameters.xmin;
ymin = parameters.ymin;
xmax = parameters.xmax;
ymax = parameters.ymax;

fig = figure(); hold on
fig.WindowState = 'maximized';
plot( AP(:,1) , AP(:,2) , '^','MarkerSize',10,'MarkerEdgeColor',[0.64,0.08,0.18],'MarkerFaceColor',[0.64,0.08,0.18] )
plot( UE(:,1) , UE(:,2) , 'o','MarkerSize',10,'MarkerEdgeColor',[0.30,0.75,0.93],'MarkerFaceColor',[0.30,0.75,0.93] )
plot( uHat(:,1),uHat(:,2),'-s','MarkerEdgeColor',[0, 0, 160]./255,'MarkerFaceColor',[0, 0, 160]./255)
legend('AP','UE','NLS')
axis equal
xticks([xmin:20:xmax])  , yticks([xmin:20:xmax])
xlim([xmin-15 xmax+15]) , ylim([ymin-15 ymax+15])
xlabel('meter','FontSize',26) , ylabel('meter','FontSize',26)
grid on
box on
for i = 1:size(uHat,1)
   text(uHat(i,1),uHat(i,2),sprintf('Iter %d ', i-1), 'fontsize',12)
end

end