function [uHat, numberOfPerformedIterations] = iterativeNLS( parameters , AP , TYPE , R , rho )
uHat = zeros(parameters.NiterMax, 2);

x = parameters.xmin/2 + (parameters.xmax-parameters.xmin)/2*rand; 
y = parameters.ymin/2 + (parameters.ymax-parameters.ymin)/2*rand; 

%% Step 1 - initial guess
uHat(1,:) = [x, y];  

eps = 10^-6; 

for iter = 2:parameters.NiterMax
    %% Step 2 - compute Jacobian matrix
    H = buildJacobianMatrixH( parameters , uHat(iter-1,:) , AP , TYPE );
    %% Step 3 - compute the observation matrix and evaluate the difference delta rho
    h = measurementModel( parameters , uHat(iter-1,:) , AP , TYPE );

    delta_rho = rho - h;
    %% Step 4 - compute the correction
    dUH = pinv(H)*delta_rho';
    %dUH = (H'*R^-1*H)^-1*H'*R^-1*delta_rho';
    %% Step 5 - update the estimate
    uHat(iter,:) = uHat(iter-1,:) + dUH';

    %% stopping criterion
    if abs(sum(uHat(iter ,:) - uHat(iter-1 ,:), 2)) < eps
        break
    end
end

numberOfPerformedIterations = iter;
