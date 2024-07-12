function rho = generateNoisyMeasurements( parameters , h , TYPE )

rho = zeros(1, parameters.numberOfAP);

if strcmp(TYPE,'TDOA')
    rho = zeros(1, parameters.numberOfAP-1);
end

switch TYPE
    case 'TOA'
        rho = h + parameters.sigmaTOA * randn(1, parameters.numberOfAP);
    case 'AOA'
        rho = h + parameters.sigmaAOA * randn(1, parameters.numberOfAP);
    case 'RSS'
        rho = h + parameters.sigmaRSS * randn(1, parameters.numberOfAP);
    case 'TDOA'
        rho = h + parameters.sigmaTDOA * randn(1, parameters.numberOfAP-1);
end

