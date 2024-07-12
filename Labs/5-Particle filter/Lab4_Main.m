clear all, clc, close all
set(0,'DefaultTextFontSize',20);
set(0,'DefaultAxesFontSize',24);
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultTextInterpreter','latex')

%% scenario settings
parameters.xmin = -100; parameters.ymin = -100;
parameters.xmax =  100; parameters.ymax =  100;

%% position of the Access Points
parameters.numberOfAP = 4;
[ AP ] = generatePositionOfAP(parameters);


%% parameters
parameters.simulationTime = 20; %s
parameters.samplingTime = 2; %s
%parameters.UEvelocity = 2; %m/s


%% generate UE trajectory
% [ UE ]= generateUEtrajectory(parameters);
load('UE_trajectory.mat')
 
fig = figure(); hold on
fig.WindowState = 'maximized';
plot( AP(:,1) , AP(:,2) , '^','MarkerSize',10,'MarkerEdgeColor',[0.64,0.08,0.18],'MarkerFaceColor',[0.64,0.08,0.18] )
plot( UE(:,1) , UE(:,2) , 'o','MarkerSize',10,'MarkerEdgeColor',[0.30,0.75,0.93],'MarkerFaceColor',[0.30,0.75,0.93] )
legend('AP','UE')
xlabel('[m]'), ylabel('[m]');
xlim([parameters.xmin parameters.xmax])
ylim([parameters.ymin parameters.ymax])
axis equal
grid on


%% generate TOA measurements
parameters.sigmaTOA = 3; %m
% [ rhoUEAP ] = generateTOAMeasurements(parameters,AP,UE);
load('rhoUEAP_4AP.mat')


%% Tracking by PF
parameters.numberOfParticles = 1000;

%i.  Initialization of particles (prior)
PR = generatePriorOfParticles(parameters);
% plotParticles( parameters , AP , PR , UE )
% plotParticles3D( PR )
% plotParticles3DSurface( PR  )
 
%Tracking
uHat = zeros(parameters.simulationTime,2);

for time=1:parameters.simulationTime
    likelihood = ones(parameters.numberOfParticles,1);
    
    %ii. evaluate the likelihood for all particles
    for a = 1:parameters.numberOfAP 
        likelihood = likelihood.*evaluateLikelihoodTOA( parameters, rhoUEAP(time,a) , AP(a,:) , PR.samples' );
    end
        %normalization
        likelihood = likelihood./sum(likelihood);
        PR.weights = likelihood';
    %      plotParticles( parameters , AP , PR , UE )
    %     plotParticles3D( PR  )
    %     plotParticles3DSurface( PR  )
    
    %iii. resampling
    indexes = resamplingAlgorithm( PR.weights , parameters.numberOfParticles );
    PR.samples = PR.samples(1:2,indexes);
%     plotParticles( parameters , AP , PR , UE )

    %iv. UE estimate
    uHat(time,:) = mean(PR.samples,2);
    
    %v. Propagation of particles
    PR = propagateParticles(parameters,PR);
%     plotParticles( parameters , AP , PR , UE )

% pause(1)

end


%% plot estimated trajectory
figure,hold on
plot( AP(:,1) , AP(:,2) , '^','MarkerSize',10,'MarkerEdgeColor',[0.64,0.08,0.18],'MarkerFaceColor',[0.64,0.08,0.18] )
plot( UE(:,1) , UE(:,2) , 'o','MarkerSize',10,'MarkerEdgeColor',[0.30,0.75,0.93],'MarkerFaceColor',[0.30,0.75,0.93] )
plot( uHat(:,1) , uHat(:,2) , '-r','Marker','s')
legend('AP','UE GT','UE est')
xlabel('[m]'), ylabel('[m]');
xlim([parameters.xmin parameters.xmax])
ylim([parameters.ymin parameters.ymax])
axis equal
grid on


%% distance error
DeltaPosition = UE(:,1:2)-uHat(:,1:2);
err = sqrt( sum ( DeltaPosition.^2,2 ));

figure
plot(parameters.samplingTime:parameters.samplingTime:parameters.samplingTime*parameters.simulationTime,err)
xlabel('time [s]'), ylabel('m'), title('Position error')
