function plotParticles3DSurface( PR  )

figure
for i = 1 :size(PR,2)
    
    x= PR(i).samples(1,:);
    y= PR(i).samples(2,:);
    z= PR(i).weights;
    
    % interpolation on regular grid
    xlin=linspace(min(x),max(x),100);
    ylin=linspace(min(y),max(y),100);
    [X,Y]=meshgrid(xlin,ylin);
    Z=griddata(x,y,z,X,Y,'cubic');
    
    % visualization
    mesh(X,Y,Z);
    axis tight; hold on
    % plot3(x,y,z,'.', 'MarkerSize',15)
    % xlabel('x')
    % ylabel('y')
    % surf(X,Y,Z),shading flat
end
zlim([-0.1 1.1])
xlabel('m')
ylabel('m')
zlabel('weight')
end