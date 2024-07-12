function H = buildJacobianMatrixH(UE , AP) %this should function return the Jacobian matrix H of numberOfAP rows and 2 columns (x and y) 

numberOfAP = 6;
H = zeros(numberOfAP-1, 2); %, 3


%id = [1, 3, 4, 5];

t = 2;
for a = 1:numberOfAP
    %i = id(a);
    %if a ~= 2 % AP2 - master
    if a == 1
        d = sqrt( sum((UE-AP(a,:)).^2 ) );
        d2 = sqrt( sum((UE-AP(2,:)).^2  ) );
        
        H(a, :) = -[(UE(1) - AP(a,1))/d - (UE(1) - AP(2,1))/d2,...
            (UE(2) - AP(a, 2))/d - (UE(2) - AP(2, 2))/d2];
%             (UE(3) - AP(a, 3))/d - (UE(3) - AP(2, 3))/d2]; 
    elseif a>2
        d = sqrt( sum((UE-AP(a,:)).^2 ) );
        d2 = sqrt( sum((UE-AP(2,:)).^2) );
        H(t, :) = -[(UE(1) - AP(a,1))/d - (UE(1) - AP(2,1))/d2,...
            (UE(2) - AP(a, 2))/d - (UE(2) - AP(2, 2))/d2];
%             (UE(3) - AP(a, 3))/d - (UE(3) - AP(2, 3))/d2]; 
        t = t+1;
    end
    %disp(H)
end
%disp("DONE")