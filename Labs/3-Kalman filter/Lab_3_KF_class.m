clear all, clc, close all
set(0,'DefaultTextFontSize',14);
set(0,'DefaultAxesFontSize',14);
set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultTextInterpreter','latex')

%% scenario settings
parameters.xmin = 0; parameters.ymin = -10;
parameters.xmax =  50; parameters.ymax =  85;

%% parameters
parameters.simulationTime = 20; %s
parameters.samplingTime = 2; %s


%% load UE trajectory
%load('UE_trajectory_1.mat')
load('UE_trajectory_2.mat')

fig = figure(); hold on
fig.WindowState = 'maximized';
plot( UE(:,1) , UE(:,2) , 'o','MarkerSize',10,'MarkerEdgeColor',[0.30,0.75,0.93],'MarkerFaceColor',[0.30,0.75,0.93] )
legend('UE')
for time = 1:parameters.simulationTime
   text(UE(time,1)+1,UE(time,2),sprintf('Time %d ', time), 'fontsize',12)
end
xlabel('[m]'), ylabel('[m]');
xlim([parameters.xmin parameters.xmax])
ylim([parameters.ymin parameters.ymax])
axis equal
grid on
box on


%% generate ABSOLUTE POSITION measurements
% rho = u + n (true pos + noise)
parameters.sigmaPos = 3; %m
R = parameters.sigmaPos^2*eye(2);

rho = zeros( parameters.simulationTime , 2 );

for time = 1 : parameters.simulationTime
    rho( time , : ) =  UE(time,:) + parameters.sigmaPos*randn(1,2);
end

fig = figure(); hold on
fig.WindowState = 'maximized';
plot( UE(:,1) , UE(:,2) , 'o','MarkerSize',10,'MarkerEdgeColor',[0.30,0.75,0.93],'MarkerFaceColor',[0.30,0.75,0.93] )
plot( rho(:,1) , rho(:,2) , '^k')
for time = 1:parameters.simulationTime
   plot( [UE(time,1),rho(time,1)] , [UE(time,2),rho(time,2)] ,'-','LineWidth', 0.5, 'color','k' )
end
legend('True','UE meas')
xlabel('[m]'), ylabel('[m]');
xlim([parameters.xmin parameters.xmax])
ylim([parameters.ymin parameters.ymax])
axis equal
grid on
box on


%% Tracking by KF
% Single measurement (e.g. GPS on position), the state is only position
% assumption: stationary UE

%initialization
UE_init = [0, 0];% don't know anything
UE_init_COV = 1000^2*eye(2);% high uncertainty
x_hat = NaN( parameters.simulationTime , 2);
P_hat = NaN( 2 , 2 , parameters.simulationTime );

fig = figure(11); hold on
fig.WindowState = 'maximized';
title('Time: ' , num2str(0) )
plot( UE(:,1) , UE(:,2) , 'o','MarkerSize',10,'MarkerEdgeColor',[0.30,0.75,0.93],'MarkerFaceColor',[0.30,0.75,0.93] )
plot( UE_init(1) , UE_init(2) , 'x','MarkerSize',10 )
plotCovariance( UE_init_COV , UE_init(1,1) , UE_init(1,2)  , 3 ,'Initialization');
axis equal
legend('True UE','UE est = init.','Cov. prior')


H = [1, 0;...
     0, 1];
Q = zeros(2,2);%5*eye(2);% enlarges the uncertainty you have at t-1; TUNE (FAST/SLOW)
%make the target move. Take it to zero -> you are averaging
F = [1, 0;...
     0, 1];

pause
close all

%update over time
for time = 1 : parameters.simulationTime

    %prediction
    if time == 1
        x_pred = UE_init';
        P_pred = UE_init_COV;
    else
        x_pred = F * x_hat(time -1, :)';
        P_pred = F * P_hat(:, :, time-1) * F' + Q;
    end
    %update
    G = P_pred * H' * inv(H * P_pred * H' + R);
    x_hat(time,:) =  x_pred + G * (rho(time, :)' - H * x_pred);
    P_hat(:,:,time) = P_pred - G * H * P_pred;

    %plot evolution
    fig = figure(11);
    fig.WindowState = 'maximized';
    plot( UE(:,1) , UE(:,2) , 'o','MarkerSize',10,'MarkerEdgeColor',[0.30,0.75,0.93],'MarkerFaceColor',[0.30,0.75,0.93] ), hold on
    plot( UE(time,1) , UE(time,2) , 'o','MarkerSize',10,'MarkerEdgeColor',[0.30,0.75,0.93],'MarkerFaceColor',[0.50,0,0] ),
    plotCovariance( P_pred  , x_pred(1) , x_pred(2)  , 3 , 'Prior');
    axis equal
    xlim([parameters.xmin parameters.xmax]) , ylim([parameters.ymin parameters.ymax])
    xlabel('[m]'), ylabel('[m]');
    legend('True UE (all)','True UE (current)','Cov. pred.')
    title('Time: ' , num2str(time) )
    pause(0.1)
    plot( rho(time,1) , rho(time,2) , '^k')
    plot( x_hat(:,1) , x_hat(:,2) , '-g','Marker','s')
    plotCovariance( P_hat(:,:,time)  , x_hat(time,1) , x_hat(time,2)  , 3 , 'Update');
    legend('True UE (all)','True UE (current)','Cov. pred.','meas.','KF est.','Cov. upd.')

    pause(0.1)
    hold off

end


% plot estimated trajectory
figure,hold on
plot( UE(:,1) , UE(:,2) , 'o','MarkerSize',10,'MarkerEdgeColor',[0.30,0.75,0.93],'MarkerFaceColor',[0.30,0.75,0.93] )
plot( x_hat(:,1) , x_hat(:,2) , '-g','Marker','s')
legend('UE true','KF est.')
xlabel('[m]'), ylabel('[m]');
xlim([parameters.xmin parameters.xmax])
ylim([parameters.ymin parameters.ymax])
axis equal
grid on


% distance error
DeltaPosition_KF = UE(:,1:2) - x_hat(:,1:2);
err_KF = sqrt( sum ( DeltaPosition_KF.^2,2 ) );

figure
plot( parameters.samplingTime:parameters.samplingTime:parameters.samplingTime*parameters.simulationTime , err_KF ), hold on
xlabel('time [s]'), ylabel('m'), title('Position error')


%% Tracking by KF
% Single measurement (e.g. GPS on position), the state is position and velocity
% assumption: constant velocity UE


%initialization
UE_init = [0, 0, 0, 0];
UE_init_COV = 100^2*eye(4);
x_hat = NaN( parameters.simulationTime , 4);
P_hat = NaN( 4 , 4 , parameters.simulationTime );

fig = figure(11); hold on
fig.WindowState = 'maximized';
title('Time: ' , num2str(0) )
plot( UE(:,1) , UE(:,2) , 'o','MarkerSize',10,'MarkerEdgeColor',[0.30,0.75,0.93],'MarkerFaceColor',[0.30,0.75,0.93] )
plot( UE_init(1) , UE_init(2) , 'x','MarkerSize',10 )
plotCovariance( UE_init_COV(1:2,1:2) , UE_init(1,1) , UE_init(1,2)  , 3 ,'Initialization');
axis equal
legend('True UE','UE est = init.','Cov. prior')


H = [1, 0, 0, 0;...
     0, 1, 0, 0]; 
driving_noise_sigma = 20;% LARGE -> BETTER FOR SHARP ANGLES(MORE REACTIVE)

L = [parameters.samplingTime^2/2, 0;...
               0,            parameters.samplingTime^2/2;...
     parameters.samplingTime,     0;...
               0,            parameters.samplingTime];
%Q = diag([2, 2, 0, 0]);
Q = driving_noise_sigma^2 * L*L';

F = [1, 0, parameters.samplingTime, 0;...
     0, 1,      0,             parameters.samplingTime;...
     0, 0,      1,                  0;...
     0, 0,      0,                  1];

pause
close all

%update over time
for time = 1 : parameters.simulationTime

    %prediction
    if time == 1
        x_pred = UE_init';
        P_pred = UE_init_COV;
    else
       x_pred = F*x_hat(time-1,:)';
       P_pred = F * P_hat(:,:, time-1)*F' + Q;
    end
    %update
    G = P_pred * H' * inv(H * P_pred * H' + R);
    x_hat(time,:) =  x_pred + G * (rho(time, :)' - H * x_pred);
    P_hat(:,:,time) = P_pred - G * H * P_pred;
    

    %plot evolution
    fig = figure(11);
    fig.WindowState = 'maximized';
    plot( UE(:,1) , UE(:,2) , 'o','MarkerSize',10,'MarkerEdgeColor',[0.30,0.75,0.93],'MarkerFaceColor',[0.30,0.75,0.93] ), hold on
    plot( UE(time,1) , UE(time,2) , 'o','MarkerSize',10,'MarkerEdgeColor',[0.30,0.75,0.93],'MarkerFaceColor',[0.50,0,0] ),
    plotCovariance( P_pred(1:2,1:2)  , x_pred(1) , x_pred(2)  , 3 , 'Prior');
    axis equal
    xlim([parameters.xmin parameters.xmax]) , ylim([parameters.ymin parameters.ymax])
    xlabel('[m]'), ylabel('[m]');
    legend('True UE (all)','True UE (current)','Cov. pred.')
    title('Time: ' , num2str(time) )
    pause(0.1)
    plot( rho(time,1) , rho(time,2) , '^k')
    plot( x_hat(:,1) , x_hat(:,2) , '-g','Marker','s')
    plotCovariance( P_hat(1:2,1:2,time)  , x_hat(time,1) , x_hat(time,2)  , 3 , 'Update');
    legend('True UE (all)','True UE (current)','Cov. pred.','meas.','KF est.','Cov. upd.')

    pause(0.1)
    hold off

end