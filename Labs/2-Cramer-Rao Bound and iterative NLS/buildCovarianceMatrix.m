function R = buildCovarianceMatrix(parameters , TYPE)

% construct a diagonal matrix
switch TYPE
    case 'TOA'
        R = parameters.sigmaTOA^2 * eye(parameters.numberOfAP);
    case 'AOA'
        R = parameters.sigmaAOA^2 * eye(parameters.numberOfAP);
    case 'RSS'
        R = parameters.sigmaRSS^2 * eye(parameters.numberOfAP);
    case 'TDOA'
        R = parameters.sigmaTDOA^2 * eye(parameters.numberOfAP-1);
end

