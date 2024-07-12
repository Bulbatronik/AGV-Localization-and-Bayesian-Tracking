function [PR] = propagateParticles(parameters,PR)

sigma_v = 1; %m

%motion model M1
PR.samples = PR.samples + parameters.samplingTime.*sigma_v.*randn(2,parameters.numberOfParticles);
PR.weights = ones(1,parameters.numberOfParticles)./parameters.numberOfParticles;
end