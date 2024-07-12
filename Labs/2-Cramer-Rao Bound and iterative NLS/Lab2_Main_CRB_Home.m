clear all, clc, close all
set(0,'DefaultTextFontSize',22)
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultTextInterpreter','latex')
%% 
set(0,'DefaultAxesFontSize',16)

%% 1. Define the localization scenario
% scenario settings
parameters.xmin = -100; parameters.ymin = -100;
parameters.xmax =  100; parameters.ymax =  100;

% position of the Access Points
parameters.numberOfAP = 4 ; % set 2, 3, 4 or 6 for specific positions, any other number will return random positions
[ AP ] = generatePositionOfAP(parameters);

% position of the UE
UE = [ 0 , 0 ];

% plot
plotScenario( parameters , AP , UE )


%% 2. Generate noisy measurements 
TYPE = 'RSS'; 

parameters.np = 2;
parameters.Pref = 10;

% build observation vector/matrix
h = measurementModel( parameters , UE , AP , TYPE ); %this function should return a row vector of numberOfAP columns of values corresponding to the true measurement rho

% measurement accuracy
parameters.sigmaTOA = 5; %m
parameters.sigmaAOA = deg2rad( 5 ); %deg
parameters.sigmaTDOA = 5; %m
parameters.sigmaRSS = 1; % dB


rho = generateNoisyMeasurements( parameters , h , TYPE ); %this function should return a row vector of numberOfAP columns of values corresponding to the true measurement rho + noise

%% 3. build covariance matrix
R = buildCovarianceMatrix( parameters , TYPE ); % this function should return a numberOfAP x numberOfAP matrix with elements on the main diagonal only. The values of such elements are the variance used to generate the noisy measurement in the previous function

%% 4. build matrix H
H = buildJacobianMatrixH( parameters , UE , AP , TYPE ); %this should function return the Jacobian matrix H of numberOfAP rows and 2 columns (x and y) 

%% 5. calculate and plot the CRB ellipse
k=3; % #sigmas
calculateEllipse( parameters , H , R , UE , AP , TYPE , k );
legend('AP','Ellipse','UE','location','best')

%% extra funny things
%% 1
close all
%remove fig=figure() from function "plotEllipse"
for i=1:7
    UE = [(rand-0.5)*2*parameters.xmax,(rand-0.5)*2*parameters.ymax];
    TYPE = 'RSS';
    H = buildJacobianMatrixH( parameters , UE , AP , TYPE );
    parameters.sigmaTOA = randi( 10 );
    R = buildCovarianceMatrix( parameters , TYPE );
    HDOP(i) = sqrt( trace( inv(H'*inv(R)*H) ) );
    calculateEllipse( parameters , H , R , UE , AP , TYPE , 3 );
    legend('AP','Ellipse','UE','location','best')
    pause(2)
end

%% 2 
%remove patch from function "plotEllipse"
close all
parameters.sigmaTOA = 10;
vecY = -100:8:100;
for i=1:length(vecY)
    UE = [ 0 , vecY(i) ];
    TYPE = 'TDOA';
    H = buildJacobianMatrixH( parameters , UE , AP , TYPE );
    R = buildCovarianceMatrix( parameters , TYPE );
    HDOP(i) = sqrt( trace( inv( H'*inv(R)*H ) ) );
    calculateEllipse( parameters , H , R , UE , AP , TYPE , 3 );
    legend('AP','Ellipse','UE','location','best')
    pause(2)
end