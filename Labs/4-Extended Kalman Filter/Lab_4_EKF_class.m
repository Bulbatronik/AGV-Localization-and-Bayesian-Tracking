clear all, clc, close all
set(0,'DefaultTextFontSize',20);
set(0,'DefaultAxesFontSize',24);
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultTextInterpreter','latex')

%% scenario settings
parameters.xmin = -50; parameters.ymin = -50;
parameters.xmax =  100; parameters.ymax =  100;

%% position of the Access Points
parameters.numberOfAP = 4;
AP(1,1) = -50 ; AP(1,2) = -50;
AP(2,1) = 50 ; AP(2,2) = -50;
AP(3,1) = 50 ; AP(3,2) = 100;
AP(4,1) = -50 ; AP(4,2) = 100;

%% parameters
parameters.simulationTime = 20; %s
parameters.samplingTime = 2; %s
%parameters.UEvelocity = 2; %m/s

%% generate UE trajectory
load('UE_trajectory_1.mat')

 
fig = figure(); hold on
fig.WindowState = 'maximized';
plot( UE(:,1) , UE(:,2) , 'o','MarkerSize',10,'MarkerEdgeColor',[0.30,0.75,0.93],'MarkerFaceColor',[0.30,0.75,0.93] ), hold on
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
plot( AP(:,1) , AP(:,2) , '^' , 'MarkerSize', 10 , 'MarkerEdgeColor' , [ 0.64 , 0.08 , 0.18 ] , 'MarkerFaceColor' , [ 0.64 , 0.08 , 0.18 ] )
legend('UE','AP')

%% generate TOA measurements
TYPE = 'TOA';
parameters.sigmaTOA = 3; %m
rho = generateTOAMeasurements(parameters,AP,UE);

% TYPE = 'AOA';
% parameters.sigmaAOA = deg2rad(3); %deg
% rho = generateAOAMeasurements(parameters,AP,UE);

R = buildCovarianceMatrix( parameters , TYPE );


%% Tracking by EKF
% TOA/AOA measurement 

%initialization
UE_init = [0,0];
UE_init_COV = diag([100^2,100^2]);
x_hat = NaN( parameters.simulationTime , 2);
P_hat = NaN( 2 , 2 , parameters.simulationTime );
% 
% fig = figure(11); hold on
% fig.WindowState = 'maximized';
% title('Time: ' , num2str(0) )
% plot( UE(:,1) , UE(:,2) , 'o','MarkerSize',10,'MarkerEdgeColor',[0.30,0.75,0.93],'MarkerFaceColor',[0.30,0.75,0.93] )
% plot( UE_init(1) , UE_init(2) , 'x','MarkerSize',10 )
% plotCovariance( UE_init_COV , UE_init(1,1) , UE_init(1,2)  , 3 ,'Initialization');
% axis equal
% legend('True UE','UE est = init.','Cov. prior')


Q = diag([ 1 , 1 ] ); % 100, 10, 2, 0.5
F = [1 0 ; 0 1];

% pause
% close all

%update over time
for time = 1 : parameters.simulationTime

    %prediction
    if time == 1
        x_pred =  UE_init';
        P_pred =  UE_init_COV;
    else
        x_pred = F * x_hat(time-1,:)';
        P_pred = F * P_hat(:,:,time-1)*F' + Q;
    end
    H = buildJacobianMatrixH(parameters, x_pred' , AP , TYPE);

    %update
    G = P_pred * H' * inv( H*P_pred*H' + R);
    x_hat(time,:) = x_pred + G * ( rho(time,:)' - measurementModel ( parameters , x_pred' , AP , TYPE)' ) ;
    P_hat(:,:,time) = P_pred - G * H * P_pred;


    %plot evolution
%     fig = figure(11);
%     fig.WindowState = 'maximized';
%     plot( UE(:,1) , UE(:,2) , 'o','MarkerSize',10,'MarkerEdgeColor',[0.30,0.75,0.93],'MarkerFaceColor',[0.30,0.75,0.93] ), hold on
%     plot( UE(time,1) , UE(time,2) , 'o','MarkerSize',10,'MarkerEdgeColor',[0.30,0.75,0.93],'MarkerFaceColor',[0.50,0,0] ),
%     plotCovariance( P_pred  , x_pred(1) , x_pred(2)  , 3 , 'Prior');
%     axis equal
%     xlim([parameters.xmin parameters.xmax]) , ylim([parameters.ymin parameters.ymax])
%     xlabel('[m]'), ylabel('[m]');
%     legend('True UE (all)','True UE (current)','Cov. pred.')
%     title('Time: ' , num2str(time) )
%     pause
%     plot( x_hat(:,1) , x_hat(:,2) , '-g','Marker','s')
%     plotCovariance( P_hat(:,:,time)  , x_hat(time,1) , x_hat(time,2)  , 3 , 'Update');
%     legend('True UE (all)','True UE (current)','Cov. pred.','KF est.','Cov. upd.')
% 
%     pause
%     hold off

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

%% Part 2 - EKF for pos and vel estimate

clear all, clc, close all
set(0,'DefaultTextFontSize',14);
set(0,'DefaultAxesFontSize',14);
set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultTextInterpreter','latex')

parameters.numberOfAP = 4;
AP(1,1) = -5000 ; AP(1,2) = -5000;
AP(2,1) = 5000 ; AP(2,2) = -5000;
AP(3,1) = 5000 ; AP(3,2) = 1000;
AP(4,1) = -5000 ; AP(4,2) = 1000;

parameters.simulationTime = 20; %s
parameters.samplingTime = 2; %s

load('UE_trajectory_3.mat')

fig = figure(11); hold on
fig.WindowState = 'maximized';
plot( UE(:,1) , UE(:,2) , 'o','MarkerSize',10,'MarkerEdgeColor',[0.30,0.75,0.93],'MarkerFaceColor',[0.30,0.75,0.93] ), hold on
txt = "Start" ;
text( UE(1,1),UE(1,2),txt)
txt2 = "End" ;
text( UE(end,1),UE(end,2),txt2)
plot( AP(:,1) , AP(:,2) , '^','MarkerSize',10,'MarkerEdgeColor',[147,0,0]./255,'MarkerFaceColor',[147,0,0]./255)
xlabel('[m]'), ylabel('[m]');
legend('UE true','AP','location','best')
axis equal


figure,hold on
plot( UE(:,3:4).*3.6)
xlabel('time')
ylabel('km/h')
legend('$v_x$','$v_y$')
title('UE velocity')

%% generate TOA measurements
parameters.simulationTime = size(UE,1);
TYPE = 'TOA';
parameters.sigmaTOA = 0.001; %m
rho = generateTOAMeasurements( parameters, AP, UE(:,1:2) );
R = buildCovarianceMatrix( parameters , 'TOA' );


%introduce BIAS
% figure
% plot(rho,'--k'),hold on
% rho(250:350,[1,2,3]) = rho(250:350,[1,2,3])+100;
% plot(rho,'-r')


%initialization
UE_init = [0, 0, 0, 0];
UE_init_COV = diag([10000^2, 10000^2, 100^2, 100^2]);
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
pause
close all

% motion model parameters
sigma_driving = 2;
L = [0.5* parameters.samplingTime^2, 0; 
                  0, parameters.samplingTime^2;
                  parameters.samplingTime,0; 
                  0, parameters.samplingTime];
Q = sigma_driving^2 *L * L';
F = [1, 0, parameters.samplingTime, 0;
     0, 1, 0, parameters.samplingTime;
     0, 0, 1, 0;
     0, 0, 0, 1];


%update over time
for time = 1 : parameters.simulationTime

    %prediction
    if time == 1
        x_pred =  UE_init';
        P_pred = UE_init_COV;
    else
        x_pred = F * x_hat(time-1,:)';
        P_pred = F * P_hat(:,:,time-1) *F' + Q;
    end
    H = buildJacobianMatrixH(parameters, x_pred(1:2)' , AP , TYPE);
    H = [H zeros(parameters.numberOfAP, 2)];%no direct measurements of the velocity

    %update
    G = P_pred * H' * inv( H*P_pred*H' + R);
    x_hat(time,:) = x_pred + G * (rho( time , : )'  - measurementModel(parameters,x_pred(1:2)',AP,TYPE)' );
    P_hat(:,:,time) = P_pred - G * H * P_pred;


%     %plot evolution
%     fig = figure(11);
%     fig.WindowState = 'maximized';
%     plot( AP(:,1) , AP(:,2) , '^','MarkerSize',10,'MarkerEdgeColor',[147,0,0]./255,'MarkerFaceColor',[147,0,0]./255), hold on
%     plot( UE(:,1) , UE(:,2) , 'o','MarkerSize',10,'MarkerEdgeColor',[0.30,0.75,0.93],'MarkerFaceColor',[0.30,0.75,0.93] ), 
%     plot( UE(time,1) , UE(time,2) , 'o','MarkerSize',10,'MarkerEdgeColor',[0.30,0.75,0.93],'MarkerFaceColor',[0.50,0,0] ),
%     plotCovariance( P_pred(1:2,1:2)  , x_pred(1) , x_pred(2)  , 3 , 'Prior');
% %     axis equal
%     xlim([UE(time,1)-100 UE(time,1)+100])
%     ylim([UE(time,2)-50 UE(time,2)+50])
%     xlabel('[m]'), ylabel('[m]');
%     legend('AP','True UE (all)','True UE (current)','Cov. pred.')
%     title('Time: ' , num2str(time) )
%     pause(0.3)
%     plot( x_hat(:,1) , x_hat(:,2) , '-g','Marker','s')
%     plotCovariance( P_hat(1:2,1:2,time)  , x_hat(time,1) , x_hat(time,2)  , 3 , 'Update');
%     legend('AP','True UE (all)','True UE (current)','Cov. pred.','KF est.','Cov. upd.')
% %     axis equal
%     xlim([UE(time,1)-100 UE(time,1)+100])
%     ylim([UE(time,2)-50 UE(time,2)+50])
%     
%     pause(0.3)
%     hold off

end

% plot estimated trajectory
figure,hold on
plot( AP(:,1) , AP(:,2) , '^','MarkerSize',10,'MarkerEdgeColor',[147,0,0]./255,'MarkerFaceColor',[147,0,0]./255)
plot( UE(:,1) , UE(:,2) , 'o','MarkerSize',10,'MarkerEdgeColor',[0.30,0.75,0.93],'MarkerFaceColor',[0.30,0.75,0.93] )
plot( x_hat(:,1) , x_hat(:,2) , '-g','Marker','s')
legend('AP','UE true','KF est.')
xlabel('[m]'), ylabel('[m]');
axis equal
grid on


% distance error
DeltaPosition_KF = UE(:,1:2) - x_hat(:,1:2);
err_KF = sqrt( sum ( DeltaPosition_KF.^2,2 ) );

figure
plot( parameters.samplingTime:parameters.samplingTime:parameters.samplingTime*parameters.simulationTime , err_KF ), hold on
xlabel('time [s]'), ylabel('m'), title('Position error')


figure
subplot(221)
plot( parameters.samplingTime:parameters.samplingTime:parameters.samplingTime*parameters.simulationTime , UE(:,1)), hold on
plot( parameters.samplingTime:parameters.samplingTime:parameters.samplingTime*parameters.simulationTime , x_hat(:,1)), hold on
xlabel('time [s]'), ylabel('m'), title('Position error $x$'), legend('true','KF est.')
subplot(222)
plot( parameters.samplingTime:parameters.samplingTime:parameters.samplingTime*parameters.simulationTime , UE(:,2)), hold on
plot( parameters.samplingTime:parameters.samplingTime:parameters.samplingTime*parameters.simulationTime , x_hat(:,2)), hold on
xlabel('time [s]'), ylabel('m'), title('Position error $y$'), legend('true','KF est.')
subplot(223)
plot( parameters.samplingTime:parameters.samplingTime:parameters.samplingTime*parameters.simulationTime , UE(:,3)), hold on
plot( parameters.samplingTime:parameters.samplingTime:parameters.samplingTime*parameters.simulationTime , x_hat(:,3)), hold on
xlabel('time [s]'), ylabel('m/s'), title('Velocity error $v_x$'), legend('true','KF est.')
subplot(224)
plot( parameters.samplingTime:parameters.samplingTime:parameters.samplingTime*parameters.simulationTime , UE(:,4)), hold on
plot( parameters.samplingTime:parameters.samplingTime:parameters.samplingTime*parameters.simulationTime , x_hat(:,4)), hold on
xlabel('time [s]'), ylabel('m/s'), title('Velocity error $v_y$'), legend('true','KF est.')
%CDF


