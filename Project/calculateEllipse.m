function calculateEllipse( parameters , H , R , UE , AP , TYPE , k )
% input k is the k-sigma ellipse value (e.g., k=3 sigma)

% CRB
C = inv( H'*inv(R)*H ); 

if rcond(C) <= 1e-4
    damp_factor = 1e-4;
    C = inv( H'*inv(R)*H + damp_factor*eye(size(H'*inv(R)*H)) );
end

% singular value decomposition
[U , S , ~] = svd( C ); %C = U S V'
S = sqrt(S);              
[~ , pos_max] = max( diag(S) );

a = k * max(max(S));  % Major semiaxis
b = k * min(max(S));  % Minor semiaxis
theta = atan2( U(2,pos_max) , U(1,pos_max) );  % orientation of the ellipse

% ellipse equation
axisTheta = 0:.01:2*pi;
ellipsePoints(1,:) = a*cos(axisTheta)*cos(theta) - b*sin(axisTheta)*sin(theta) + UE(1);
ellipsePoints(2,:) = a*cos(axisTheta)*sin(theta) + b*sin(axisTheta)*cos(theta) + UE(2);

% Plot ellipse
plotEllipse( parameters , AP , ellipsePoints , UE , TYPE )
