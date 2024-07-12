function [x_hat] = KFTrack(AP,rho,sampleNumber,init,R,sigma_driving)

numberOfAP = 6;
samplingTime = 0.1;
UE_init = [init(1), init(2), 0, 0, 0, 0];
UE_init_COV = diag([0.04^2, 0.04^2, 0.04^2, 0.04^2, 0.04^2, 0.04^2]);
%figure;
%hold on
%plotCovariance( UE_init_COV(1:2,1:2) , UE_init(1,1) , UE_init(1,2)  , 3 ,'Initialization');

x_hat = NaN( sampleNumber , 6);
P_hat = NaN( 6 , 6 , sampleNumber );



L = [(samplingTime^3)/6, 0; 
      0,(samplingTime^3)/6;
      (samplingTime^2)/2,0; 
      0,(samplingTime^2)/2;
      samplingTime, 0;
      0, samplingTime     ];
              
% Q = sigma_driving^2 *(L * L');
% 
% F = [1, 0, samplingTime, 0;
%      0, 1, 0, samplingTime;
%      0, 0, 1, 0;
%      0, 0, 0, 1];

 %update over time
 stopFlag = 0;

for time = 1 : sampleNumber

    if time > 12

        
        speedmean = norm(mean(x_hat(time-11:time-1,3:4)));
        
        
    else
       
        speedmean = 1000;
    end

 if speedmean < 1e-2 

     F = [1, 0, 0, 0, 0, 0;
          0, 1, 0, 0, 0, 0;
          0, 0, 0, 0, 0, 0;
          0, 0, 0, 0, 0, 0;
          0, 0, 0, 0, 0, 0;
          0, 0, 0, 0, 0, 0];

     Q = 0;
     stopFlag = 1;

 else
%      F = [1, 0, samplingTime, 0;
%      0, 1, 0, samplingTime;
%      0, 0, 1, 0;
%      0, 0, 0, 1];

     F = [1, 0, samplingTime, 0, (samplingTime^2)/2, 0;
          0, 1, 0, samplingTime, 0, (samplingTime^2)/2;
          0, 0, 1, 0, samplingTime, 0;
          0, 0, 0, 1, 0, samplingTime;
          0, 0, 0, 0, 1, 0;
          0, 0, 0, 0, 0, 1];
     
     Q = sigma_driving^2 *(L * L');

 end

    %prediction
    if time == 1
        x_pred =  UE_init';
        P_pred = UE_init_COV;
    else
        x_pred = F * x_hat(time-1,:)';
        P_pred = F * P_hat(:,:,time-1) *F' + Q;
    end

    H = buildJacobianMatrixH(x_pred(1:2)' , AP);
    H = [H zeros(numberOfAP-1, 4)];%no direct measurements of the velocity 

    

    %update
    G = P_pred * H' * inv( H*P_pred*H' + R);
    if stopFlag == 1
        G = G .* 0;
    end

    x_hat(time,:) = x_pred + G * (rho( : , time )'  - measurementModel(x_pred(1:2)',AP) )';
    P_hat(:,:,time) = P_pred - G * H * P_pred;


 %plotCovariance( P_pred(1:2,1:2)  , x_pred(1) , x_pred(2)  , 3 , 'Prior');
 %plotCovariance( P_hat(1:2,1:2,time)  , x_hat(time,1) , x_hat(time,2)  , 3 , 'Update');
 %hold on;
   
end
 


end

