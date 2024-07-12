function H = buildJacobianMatrixH( parameters , UE , AP , TYPE ) %this should function return the Jacobian matrix H of numberOfAP rows and 2 columns (x and y) 

H = zeros(parameters.numberOfAP, 2);
if strcmp(TYPE,'TDOA')
    H = zeros(parameters.numberOfAP-1, 2);
end

for a = 1:parameters.numberOfAP
    switch TYPE
        case 'TOA'
            d = sqrt( sum([UE-AP(a,:)].^2 , 2 ) );
            H(a, :) = [(UE(1)- AP(a,1))/d, (UE(2)- AP(a,2))/d];
        case 'AOA'
            d2 = sum([UE-AP(a,:)].^2 , 2 );
            H(a, :) = [-(UE(2)- AP(a,2))/d2, (UE(1)- AP(a,1))/d2];
        case 'RSS'
            d2 = sum([UE-AP(a,:)].^2 , 2 );
            H(a, :) = -10*parameters.np/log(10)*[(UE(1)- AP(a,1))/d2, (UE(2)- AP(a,2))/d2];
        case 'TDOA'
            d1 = sqrt( sum([UE-AP(1,:)].^2 , 2 ) );
            if a>1
                d = sqrt( sum([UE-AP(a,:)].^2 , 2 ) );
                H(a-1, :) = [(UE(1) - AP(a,1))/d - (UE(1) - AP(1,1))/d1, (UE(2) - AP(a, 2))/d - (UE(2) - AP(1, 2))/d1]; 
            end
    end
end