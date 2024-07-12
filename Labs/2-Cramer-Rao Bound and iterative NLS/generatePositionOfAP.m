function [ AP ] = generatePositionOfAP(parameters)

nAP = parameters.numberOfAP;
xmin = parameters.xmin/2;
ymin = parameters.ymin/2;
xmax = parameters.xmax/2;
ymax = parameters.ymax/2;

if nAP == 2
    AP(1,1) = xmin ; AP(1,2) = ymin;
    AP(2,1) = xmax ; AP(2,2) = ymax;
elseif nAP == 3
    AP(1,1) = xmin ; AP(1,2) = ymin;
    AP(2,1) = xmax ; AP(2,2) = ymin;
    AP(3,1) = (xmin+xmax)/2 ; AP(3,2) = ymax;
elseif nAP == 4
    AP(1,1) = xmin ; AP(1,2) = ymin;
    AP(2,1) = xmax ; AP(2,2) = ymin;
    AP(3,1) = xmax ; AP(3,2) = ymax;
    AP(4,1) = xmin ; AP(4,2) = ymax;
elseif nAP == 6
    AP(1,1) = xmin ; AP(1,2) = ymin;
    AP(2,1) = xmax ; AP(2,2) = ymin;
    AP(3,1) = xmax ; AP(3,2) = ymax;
    AP(4,1) = xmin ; AP(4,2) = ymax;
    AP(5,1) = (xmin+xmax)/2 ; AP(5,2) = ymin;
    AP(6,1) = (xmin+xmax)/2 ; AP(6,2) = ymax;
else
    warning('AP randomly placed')
    AP(1:nAP,1) = 4*( randi(xmax,nAP,1) - (xmax)/2 );
    AP(1:nAP,2) = 4*( randi(ymax,nAP,1) - (ymax)/2 );
end


end