clear all, clc, close all
set(0,'DefaultTextFontSize',18)
set(0,'DefaultLineLineWidth',1.2);
set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultAxesFontSize',16)


%% Create a Satellite Scenario
% startTime = datetime( 2022 , 5 , 13 , 15 , 10 , 0 );
startTime = datetime( 'now' );
numHours = 24;
stopTime = startTime + hours( numHours );
sampleTime = 60*15; %s
satScenario = satelliteScenario( startTime , stopTime , sampleTime );

%% play
%satScenario.play

%% Add Satellites to the Scenario
%https://www.celestrak.com/NORAD/elements/

%sat = satellite( satScenario , "GPSsatellites.tle" );
% sat = satellite( satScenario , "GNSSsatellites.tle" );
sat = satellite( satScenario , "Galileo.tle" );
%sat = satellite( satScenario , "Beidou.tle" );

%show(sat)
%groundTrack( sat )
%satScenario.play


%% Add Ground Stations
%Add the Madrid and Canberra Deep Space Communications Complexes as the ground stations of interest, and specify their latitudes and longitudes.

name = ["Madrid Deep Space Communications Complex", "Canberra Deep Space Communications Complex"];
GroundStation_lat = [40.43139, -35.40139];
GroundStation_lon = [-4.24806, 148.98167];
groundStationPos = groundStation( satScenario , "Name", name , "Latitude", GroundStation_lat , "Longitude", GroundStation_lon );


%% check satellites in visibility now
time =  datetime( 'now' )+hours(0.5);
[ satPos_now , satVel_now ] = states( sat , time , 'CoordinateFrame', 'ECEF' );

% Calculate satellite visibilities from a given receiver position
receiverCoordinates = [45.47868, 9.23253, 15]; %DEIB coordinates (Lat, Lon, Alt)
receiverVelocity = [0 0 0];
% figure, geoplot( receiverCoordinates(1) , receiverCoordinates(2) , 'x')

maskAngle = 10; %deg on elevation
[ satAz , satEl , VisibleSatFLAG ] = lookangles( receiverCoordinates , squeeze(satPos_now)' , maskAngle);

% number and names of all satellites in visibility
clc
fprintf('%d satellites visible at %s \n',nnz( VisibleSatFLAG ) , time );
fprintf('Names: \n')
fprintf('%s  \n', [ sat( VisibleSatFLAG ).Name ] )

satIndexes = 1:numel(VisibleSatFLAG);
figure
sp = skyplot( satAz(VisibleSatFLAG) , satEl(VisibleSatFLAG) , LabelData=satIndexes(VisibleSatFLAG) );

%% plot visibility condition in next hours
secondsPerHour = 3600;
timeElapsed = 0:sampleTime:(secondsPerHour*numHours);
timeAxis = startTime + seconds(timeElapsed);

numSats = size(sat,2);
numSamples = numel(timeAxis);
satPos_time = zeros(numSamples,numSats);
satVel_time = zeros(numSamples,numSats);
satAz = zeros(numSamples,numSats);
satEl = zeros(numSamples,numSats);
VisibleSatFLAG = false(numSamples,numSats);
satIDs = 1:numSats;
satellitePosOverTime = zeros(numSats,3,numSamples);
satelliteVelOverTime = zeros(numSats,3,numSamples);


for ii = 1:numSamples
    [ satPos_time , satVel_time ] = states( sat , timeAxis(ii) , "CoordinateFrame","ECEF");
    [ satAz(ii,:) , satEl(ii,:) , VisibleSatFLAG(ii,:) ] = lookangles( receiverCoordinates , squeeze(satPos_time)' , maskAngle);
    satellitePosOverTime(:,:,ii) = squeeze(satPos_time)';
    satelliteVelOverTime(:,:,ii) = squeeze(satVel_time)';
end

visPlotData = double(VisibleSatFLAG);
visPlotData( visPlotData == false ) = NaN; % Hide invisible satellites.
visPlotData = visPlotData + (0:numSats-1); 

figure
subplot(121)
plot( timeAxis , visPlotData ,".k" )
yticks(1:numSats)
grid on
ylabel( "Satellite ID" ), ylim( [0 , numSats] )
xlabel( "Time" )
title( "Satellite Visibility Chart" )

subplot(122)
plot( timeAxis , sum(VisibleSatFLAG,2) , '-k')
grid on
xlabel( "Time" )
ylabel( "Number of visible satellites" )
ylim( [0 max( sum(VisibleSatFLAG,2) )+1] )


%% calculate pseudoranges
allPseudoranges = NaN(numSats,numSamples);
allPseudovelocities = NaN(numSats,numSamples);

for ii = 1:numSamples
    satPos_time = satellitePosOverTime(:,:,ii);
    satVel_time = satelliteVelOverTime(:,:,ii);
    rangeStd = 5; %[m]
    rangeRateStd = 1; %[m/s]
    [ allPseudoranges(:,ii) , allPseudovelocities(:,ii) ] = pseudoranges(receiverCoordinates , satPos_time, receiverVelocity, satVel_time , 'RangeAccuracy', rangeStd,  'RangeRateAccuracy', rangeRateStd );

end

figure()
plot( timeAxis , allPseudoranges' )
legend( sat.Name ) 
ylabel('distance [m]')
box on

%% Estimate receiver position from pseudoranges
estimatedPosition = NaN(numSamples,3);
estimatedVelocity = NaN(numSamples,3);
hdop = NaN(numSamples,1);
vdop = NaN(numSamples,1);

for ii = 1:numSamples
    currentPseudoranges = allPseudoranges(:,ii);
    currentPseudovelocities = allPseudovelocities(:,ii);
    VisibleSatFLAG_time = VisibleSatFLAG(ii,:);
    satPos_time = satellitePosOverTime(:,:,ii);
    satVel_time = satelliteVelOverTime(:,:,ii);
    
    [ estimatedPosition(ii,:) , estimatedVelocity(ii,:) , hdop(ii,:) , vdop(ii,:) ] = ...
        receiverposition( currentPseudoranges(VisibleSatFLAG_time) , satPos_time(VisibleSatFLAG_time,:) , currentPseudovelocities(VisibleSatFLAG_time) , satVel_time(VisibleSatFLAG_time,:) );

end

figure
geoplot( receiverCoordinates(1) , receiverCoordinates(2) ,'o'), hold on
geoplot( estimatedPosition(:,1) , estimatedPosition(:,2) , '.')
legend("Ground Truth","Estimate")


estimatedPosition_NED  = lla2ned( estimatedPosition , receiverCoordinates , "ellipsoid" );
truePosition_NED = lla2ned( receiverCoordinates , receiverCoordinates , "ellipsoid");

figure
plot( timeAxis , abs(estimatedPosition_NED-truePosition_NED) )
legend("x","y","z")
ylabel("Absolute error (m)")
title("Position (NED) Error")


figure
plot3( truePosition_NED(1) , truePosition_NED(2) , truePosition_NED(3) ,'ok'), hold on
plot3( estimatedPosition_NED(:,1) , estimatedPosition_NED(:,2) , estimatedPosition_NED(:,3) ,'.r')
box on
legend('True','Est.')
xlabel('[m]'), ylabel('[m]'), zlabel('[m]')

figure
plot(  timeAxis , hdop ), hold on
plot(  timeAxis , vdop )
legend("HDOP","VDOP")
title("GDOP")

%% error metrics
truePos_NED_vec = truePosition_NED.*ones(size(estimatedPosition_NED));
%MAE
MAE = mae (estimatedPosition_NED-truePos_NED_vec);
%RMSE - mode 1
RMSE_1 = sqrt( immse(estimatedPosition_NED,truePos_NED_vec) );
%RMSE - mode 2
RMSE_2 = sqrt ( mean( (estimatedPosition_NED(:) - truePos_NED_vec(:)).^2 ) );

clc
fprintf('MAE %d \n', MAE)
fprintf('RMSE %d \n', RMSE_1)
fprintf('RMSE %d \n', RMSE_2)


err = sqrt(sum( (estimatedPosition_NED - truePos_NED_vec).^2 , 2 ) );
figure
cdfplot( err )
xlabel('error [m]')
ylabel('Probability')
title('CDF of location error')

figure
boxplot( estimatedPosition_NED-truePos_NED_vec , 'Labels',{'x','y','z'} )
ylabel('Location error [m]')


%% comparison of accuracy on the Earth 

LatVec =  -80:10:80 ;
receiverCoordinates = [ LatVec ; zeros(size(LatVec)) ;  zeros(size(LatVec)) ]';
figure, geoplot( receiverCoordinates(:,1) , receiverCoordinates(:,2) , 'x')

avg_hdop = NaN(1,size(receiverCoordinates,1));
avg_vdop = NaN(1,size(receiverCoordinates,1));
for i = 1:size(receiverCoordinates,1)
    hdop =  NaN(1,numSamples);
    vdop = NaN(1,numSamples);
for ii = 1:numSamples
    
    satPos_time = satellitePosOverTime(:,:,ii);
    satVel_time = satelliteVelOverTime(:,:,ii);
    rangeStd = 5; %[m]
    rangeRateStd = 1; %[m/s]
    [ currentPseudoranges , currentPseudovelocities ] = pseudoranges(receiverCoordinates(i,:) , satPos_time, receiverVelocity, satVel_time , 'RangeAccuracy', rangeStd,  'RangeRateAccuracy', rangeRateStd );

    VisibleSatFLAG_time = VisibleSatFLAG(ii,:);
    satPos_time = satellitePosOverTime(:,:,ii);
    satVel_time = satelliteVelOverTime(:,:,ii);
    
    [ ~ , ~ , hdop(ii) , vdop(ii) ] = ...
        receiverposition( currentPseudoranges(VisibleSatFLAG_time) , satPos_time(VisibleSatFLAG_time,:) , currentPseudovelocities(VisibleSatFLAG_time) , satVel_time(VisibleSatFLAG_time,:) );

end
avg_hdop(i) = mean(hdop);
avg_vdop(i) = mean(vdop);

end

figure
plot( receiverCoordinates(:,1) , avg_hdop ), hold on
plot( receiverCoordinates(:,1) , avg_vdop )
xlabel('Latitude')
ylabel('GDOP')
legend('Horizontal','Vertical')


%% effect of satellite obstruction / malfunctioning
estimatedPosition = NaN(numSamples,3);
estimatedVelocity = NaN(numSamples,3);
hdop = NaN(numSamples,1);
vdop = NaN(numSamples,1);
estimatedPosition2 = NaN(numSamples,3);
estimatedVelocity2 = NaN(numSamples,3);
hdop2 = NaN(numSamples,1);
vdop2 = NaN(numSamples,1);
allPseudoranges = NaN(numSats,numSamples);
allPseudovelocities = NaN(numSats,numSamples);
numVisibleSat = NaN(numSamples);
numVisibleSat2 = NaN(numSamples);


receiverCoordinates = [45.47868, 9.23253, 15]; %DEIB coordinates (Lat, Lon, Alt)
receiverVelocity = [0 0 0];


for ii = 1:numSamples
    rng(ii)
    satPos_time = satellitePosOverTime(:,:,ii);
    satVel_time = satelliteVelOverTime(:,:,ii);
    rangeStd = 5; %[m]
    rangeRateStd = 1; %[m/s]
    [ currentPseudoranges , currentPseudovelocities ] = pseudoranges(receiverCoordinates , satPos_time, receiverVelocity, satVel_time , 'RangeAccuracy', rangeStd,  'RangeRateAccuracy', rangeRateStd );

    VisibleSatFLAG_time = VisibleSatFLAG(ii,:);
    VisibleSatFLAG_time2 = VisibleSatFLAG_time;
    for i = 1 : size( VisibleSatFLAG_time2 , 2) 
        if VisibleSatFLAG_time2(i) == 1 && nnz(VisibleSatFLAG_time2)>3
            if rand > 0.8
                VisibleSatFLAG_time2(i) = 0;
            end
        end
    end
    numVisibleSat(ii) = sum(VisibleSatFLAG_time);
    numVisibleSat2(ii) = sum(VisibleSatFLAG_time2);
    satPos_time = satellitePosOverTime(:,:,ii);
    satVel_time = satelliteVelOverTime(:,:,ii);
    
     [ estimatedPosition(ii,:) , estimatedVelocity(ii,:) , hdop(ii,:) , vdop(ii,:) ] = ...
        receiverposition( currentPseudoranges(VisibleSatFLAG_time) , satPos_time(VisibleSatFLAG_time,:) , currentPseudovelocities(VisibleSatFLAG_time) , satVel_time(VisibleSatFLAG_time,:) );

    satPos_time2 = satellitePosOverTime(:,:,ii);
    satVel_time2 = satelliteVelOverTime(:,:,ii);
    
     [ estimatedPosition2(ii,:) , estimatedVelocity2(ii,:) , hdop2(ii,:) , vdop2(ii,:) ] = ...
        receiverposition( currentPseudoranges(VisibleSatFLAG_time2) , satPos_time(VisibleSatFLAG_time2,:) , currentPseudovelocities(VisibleSatFLAG_time2) , satVel_time(VisibleSatFLAG_time2,:) );

end


estimatedPosition_NED  = lla2ned( estimatedPosition , receiverCoordinates , "ellipsoid" );
truePosition_NED = lla2ned( receiverCoordinates , receiverCoordinates , "ellipsoid");
estimatedPosition_NED2 = lla2ned( estimatedPosition2 , receiverCoordinates , "ellipsoid" );

err = sqrt(sum( (estimatedPosition_NED - truePos_NED_vec).^2 , 2 ) );
err2 = sqrt(sum( (estimatedPosition_NED2 - truePos_NED_vec).^2 , 2 ) );

figure
cdfplot( err ), hold on
cdfplot( err2)
xlabel('error [m]')
ylabel('Probability')
title('CDF of location error')
legend('Full Vis.','Reduced Vis.')

figure 
plot( timeAxis , numVisibleSat ), hold on
plot( timeAxis , numVisibleSat2 )
legend('Full Vis.','Reduced Vis.')
title('Number of visible satellites')