function likelihood = evaluateLikelihoodRSS( parameters , rho , AP , evaluationPoint , np , Pref)

evaluationDistance = sqrt(sum([evaluationPoint-AP].^2,2)); 
evaluationRho =  Pref - 10*np*log10(evaluationDistance) ;
            
argument =  rho - evaluationRho  ;

likelihood = 1/sqrt(2*pi*parameters.sigmaRSS.^2)*exp(-0.5* argument^2 / parameters.sigmaRSS.^2);

end