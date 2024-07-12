function plotParticles3D( PR  )

figure
for i = 1 : size(PR,2)
    plot3( PR(i).samples(1,:) , PR(i).samples(2,:) , PR(i).weights , '.')
    hold on
end

grid on
zlim([-0.1 1.1])
xlabel('m')
ylabel('m')
zlabel('weight')
end