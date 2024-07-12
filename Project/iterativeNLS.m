function [uHat, numberOfPerformedIterations] = iterativeNLS( NiterMax , AP , rho, sample_number,R )
uHat = zeros(NiterMax, 2);%XYZ

x = 2;
y = 2;
z = 4;


uHat(1,:) = [x,y]; %init  , z


eps = 10^-6; 

for iter = 2:NiterMax
    H = buildJacobianMatrixH(uHat(iter-1,:) , AP);%CHANGED
    h = measurementModel(uHat(iter-1,:) , AP);
    
    delta_rho = rho - h';
    %dUH =(inv(H'*H))*H'*delta_rho;
    dUH = (H'*R^-1*H)^-1*H'*R^-1*delta_rho;
    uHat(iter,:) = uHat(iter-1,:) + dUH';
    
    flag = abs(sum(uHat(iter ,:) - uHat(iter-1 ,:), 2));
    
    if flag < eps
        
       
        
        break;
    end
end

numberOfPerformedIterations = iter;
