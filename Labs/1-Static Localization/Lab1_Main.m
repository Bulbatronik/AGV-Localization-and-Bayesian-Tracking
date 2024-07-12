clear all, clc, close all
set(0,'DefaultTextFontSize',22)
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultAxesFontSize',16)

%% 1. Define the localization scenario
% scenario settings
parameters.xmin = -200; parameters.ymin = -200;
parameters.xmax =  200; parameters.ymax =  200;

% position of the Access Points
parameters.numberOfAP = 4 ; % set 2, 3, 4 or 6 for specific positions, any other number will return random positions
[ AP ] = generatePositionOfAP(parameters);

% position of the UE
UE = [ 0 , -30 ];

% plot scenario
plotScenario( parameters , AP , UE )


%% 2. Create a grid of evaluation points to sample the scenario
x = linspace( parameters.xmin , parameters.xmax , 100 ) ;
y = linspace( parameters.ymin , parameters.ymax , 100 ) ;

%% 3. Compute the likelihood of the measurement
% definition of measurement accuracies
parameters.sigmaTOA = 10; % m
parameters.sigmaAOA = deg2rad(5); % deg
parameters.sigmaTDOA = 5; % m
parameters.sigmaRSS = 5; % dB

% compute likelihood for each AP in each evaluation point
likelihood = zeros(parameters.numberOfAP,length(x),length(y));
TYPE = 'TDOA';

for a = 1:parameters.numberOfAP
    % compute rho
    switch TYPE
        case 'TOA'
            rho_True = sqrt( sum([ UE - AP(a,:) ].^2 , 2 ) ) ;
        case 'AOA'
            rho_True = atan2( ( UE(2)-AP(a,2) ) , ( UE(1)-AP(a,1) ) );
        case 'RSS'
            np = 2;
            Pref = 10;
            rho_True = Pref - 10*np*log10( sqrt( sum([ UE - AP(a,:) ].^2 , 2 ) ) );
        case 'TDOA'
            if a>1
                rho_True = sqrt( sum([ UE - AP(a,:) ].^2 , 2 ) ) - sqrt(sum([UE-AP(1,:)].^2 , 2 ) );
            end
    end
    % evaluate the likelihood in each evaluation point
    for i=1:1:length(x)
        for j=length(y):-1:1
            switch TYPE
                case 'TOA'
                    likelihood(a,i,j) = evaluateLikelihoodTOA( parameters, rho_True , AP(a,:) , [x(i),y(j)] );
                case 'AOA'
                    likelihood(a,i,j) = evaluateLikelihoodAOA( parameters, rho_True , AP(a,:) , [x(i),y(j)] );
                case 'RSS'
                    likelihood(a,i,j) = evaluateLikelihoodRSS( parameters, rho_True , AP(a,:) , [x(i),y(j)] , np , Pref);
                case 'TDOA'
                    if a>1
                        likelihood(a,i,j) = evaluateLikelihoodTDOA( parameters, rho_True , AP(1,:) , AP(a,:) , [x(i),y(j)] ); % AP1 is the master AP
                    end
            end
        end %j
    end %i
end %a

%% 4. Plot the likelihood for each AP
plot2Dlikelihood( parameters, AP , UE , x , y , likelihood , TYPE)
plot3Dlikelihood( parameters, AP , UE , x , y , likelihood , TYPE)

%% 5. Compute the maximum likelihood and plot it
% compute ML
maximumLikelihood = ones( length(x) , length(y) );
if strcmp(TYPE,'TOA') || strcmp(TYPE,'AOA') || strcmp(TYPE,'RSS')
    for a = 1:parameters.numberOfAP
        maximumLikelihood = maximumLikelihood.*squeeze(likelihood(a,:,:));
    end
elseif  strcmp(TYPE,'TDOA')
    for a = 2:parameters.numberOfAP
        maximumLikelihood = maximumLikelihood.*squeeze(likelihood(a,:,:));
    end
end
maximumLikelihood = maximumLikelihood./sum(sum(maximumLikelihood)); %normalization

%plot
plotMaximumlikelihood( parameters, AP , UE , x , y , maximumLikelihood , TYPE)

%% 6. Evaluate the UE estimated position
maxValue = max( maximumLikelihood(:) );
[x_ind, y_ind] = find( maximumLikelihood == maxValue );
u_est = [x(x_ind), y(y_ind)];
% it should be the closest evaluation point to UE as we have no noise
fprintf('\n The estimated UE position is [ %.3f %.3f ] \n', u_est)