function [rho] = generateTOAMeasurements(parameters, AP, UE)

T = parameters.simulationTime;%s
NA = parameters.numberOfAP;
sigma = parameters.sigmaTOA;

rho = zeros(T, NA);
for time=1:T
    rho(time,:) = sqrt(sum([UE(time,:)-AP].^2, 2))';
end

