function plotParticles(parameters,AP,PR,UE)


% fig = figure(); 
% fig.WindowState = 'maximized';
clf
hold on
plot( AP(:,1) , AP(:,2) , '^','MarkerSize',10,'MarkerEdgeColor',[0.64,0.08,0.18],'MarkerFaceColor',[0.64,0.08,0.18] )
plot( UE(:,1) , UE(:,2) , 'o','MarkerSize',10,'MarkerEdgeColor',[0.30,0.75,0.93],'MarkerFaceColor',[0.30,0.75,0.93] )
scatter( PR.samples(1,:),PR.samples(2,:),10,'filled','markerfacecolor',[0.75 0.75 0.75])
legend('AP','UE','PR')
xlabel('[m]'), ylabel('[m]');
xlim([parameters.xmin parameters.xmax])
ylim([parameters.ymin parameters.ymax])
axis equal
grid on
box on
end