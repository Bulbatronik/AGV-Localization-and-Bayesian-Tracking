clear all, clc, close all
set(0,'DefaultTextFontSize',22)
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultAxesFontSize',16)

%% 1. Define the localization scenario
% scenario settings
parameters.xmin = -100; parameters.ymin = -100;
parameters.xmax =  100; parameters.ymax =  100;

% position of the Access Points
parameters.numberOfAP = 4 ; % set 2, 3, 4 or 6 for specific positions, any other number will return random positions
[ AP ] = generatePositionOfAP(parameters);

% position of the UE
UE = [ -20 , 20 ];

% plot
plotScenario( parameters , AP , UE )


%% 2. Generate noisy measurements 
TYPE = 'TDOA'; 

parameters.np = 2;
parameters.Pref = 10;

% build observation vector/matrix
h = measurementModel( parameters , UE , AP , TYPE );

% generate noisy measurements
parameters.sigmaTOA = 5; %m
parameters.sigmaAOA = deg2rad( 5 ); %deg
parameters.sigmaTDOA = 5; %m
parameters.sigmaRSS = 5; % dB


rho = generateNoisyMeasurements( parameters , h , TYPE );
R = buildCovarianceMatrix(parameters , TYPE);

%% 3. Implement NLS 
parameters.NiterMax = 100;
[ uHat , numberOfPerformedIterations ] = iterativeNLS( parameters , AP , TYPE , R , rho ); %this function should return a NiterMax x 2 vector of the estimates provided by the NLS
uHat = uHat( 1:numberOfPerformedIterations , : ); % this is the final estimate of NLS at the last iteration

%% 4. plot the NLS iterations 
plotNLSiterations( parameters , AP , UE , uHat , TYPE )
numberOfPerformedIterations

%% 5. perform MonteCarlo simulations
parameters.numberOfMC = 1000;
uHatMC = zeros( parameters.numberOfMC , 2 );
for imc = 1 : parameters.numberOfMC

    rho = generateNoisyMeasurements( parameters , h , TYPE );
    [ uHat , numberOfPerformedIterations ] = iterativeNLS( parameters , AP , TYPE , R , rho );
    uHatMC( imc , : ) = uHat( numberOfPerformedIterations , : );

end

%plot all the estimated values
close all
fig = figure; hold on
fig.WindowState = 'maximized';
plot( AP(:,1) , AP(:,2) , '^','MarkerSize',10,'MarkerEdgeColor',[0.64,0.08,0.18],'MarkerFaceColor',[0.64,0.08,0.18] )
plot( UE(:,1) , UE(:,2) , 'o','MarkerSize',10,'MarkerEdgeColor',[0.30,0.75,0.93],'MarkerFaceColor',[0.30,0.75,0.93] )
plot( uHatMC(:,1) , uHatMC(:,2) ,'.','color',[255,174,0]./255)
legend( 'AP' , 'UE' , 'NLS' )
axis equal
xticks( [parameters.xmin:20:parameters.xmax] )  , yticks( [parameters.xmin:20:parameters.xmax] )
xlim( [parameters.xmin-15 parameters.xmax+15] ) , ylim( [parameters.ymin-15 parameters.ymax+15] )
xlabel('meter','FontSize',26) , ylabel('meter','FontSize',26)
grid on
box on
title(['Estimated UE positions by ',num2str(parameters.numberOfMC),' Monte Carlo for ',num2str(TYPE),' measurements '],'Interpreter','Latex')

%% 6. CHECK CRB
%remove fig=figure() from function "plotEllipse"
k = 3; %3sigma
H = buildJacobianMatrixH( parameters , UE , AP , TYPE ); 
calculateEllipse( parameters , H , R , UE , AP , TYPE , k );


%% Extra for fun
NPos = 5;
parameters.numberOfMC = 1000;
uHatMC = zeros( parameters.numberOfMC , 2 );

close all
fig = figure; hold on
fig.WindowState = 'maximized';
% plot( AP(:,1) , AP(:,2) , '^','MarkerSize',10,'MarkerEdgeColor',[147,0,0]./255,'MarkerFaceColor',[147,0,0]./255)
for idx = 1:NPos
    UE(1) = randi( [parameters.xmin,parameters.xmax] );
    UE(2) = randi( [parameters.ymin,parameters.ymax] );
    for imc = 1:parameters.numberOfMC
        [ h ] = measurementModel( parameters , UE , AP , TYPE );
        [ rho ] = generateNoisyMeasurements( parameters , h , TYPE );
        [ uHat , numberOfPerformedIterations ] = iterativeNLS( parameters , AP , TYPE , R , rho );
        uHatMC( imc , : ) = uHat( numberOfPerformedIterations , : );
    end
    plot( UE(1) , UE(2) , 'o','MarkerSize',10,'MarkerEdgeColor',[0.30,0.75,0.93],'MarkerFaceColor',[0.30,0.75,0.93])
    plot( uHatMC(:,1) , uHatMC(:,2) ,'.','color',[255,174,0]./255)
    xticks([parameters.xmin:20:parameters.xmax])  , yticks([parameters.xmin:20:parameters.xmax])
    xlim([parameters.xmin-15 parameters.xmax+15]) , ylim([parameters.ymin-15 parameters.ymax+15])
    xlabel('meter','FontSize',26) , ylabel('meter','FontSize',26)
    grid on
    box on
    axis equal
    k = 3; %3sigma
    [ H ] = buildJacobianMatrixH( parameters , UE , AP , TYPE ); 
    calculateEllipse( parameters , H , R , UE , AP , TYPE , k );
pause
end