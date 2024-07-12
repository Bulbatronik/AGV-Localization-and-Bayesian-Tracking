clear;
clc;
load('rho_TDOA_experiment1.mat');
%% 
%samples per tag
samples{1} = 1:800;
samples{2}= 1:2008;
samples{3} = 1:2006;
samples{4} = 1:2014;

%% 
close all

for j = 1:4
    for i = 1:5 %measurements
        signal = rho{j}(i,:);
        if j == 4 && i == 2
        figure;
        subplot(2,1,1);
        plot(signal);
        title('Original Signal');
        xlim([1 2010]);
        end

        signal =inpaint_nans(signal);
        rhoIm{j}(i,:) = signal;
    
        [B tf]= rmoutliers(signal,'movmedian',50,'SamplePoints',samples{j});
        signal(tf) = NaN;
        signal =inpaint_nans(signal);
        
            

        [B tf]= rmoutliers(signal,'movmedian',10,'SamplePoints',samples{j});
        signal(tf) = NaN;
        signal =inpaint_nans(signal);
      
    
        order = 3;
        framelen = 5; 
    
        sgf = sgolayfilt(signal,order,framelen);
    
        [B tf]= rmoutliers(sgf,'movmedian',3,'SamplePoints',samples{j});
        sgf(tf) = NaN;
        sgf = inpaint_nans(sgf);
        
    
        filtered{j}(i,:) = sgf;%smoothdata(sgf,'sgolay');

        if j == 4 && i == 1
            
            subplot(2,1,2);
            plot(filtered{j}(i,:));
            title('Impainted and Spike Removed Signal');
            xlim([1 2010]);
        end
  
    end
end


%% NLS and EKF
load('AP.mat')
load('smoothData.mat');
NiterMax = 100;
sigma_driving = 0.01:0.01:0.1; %0.08 optimum

%for k = 1:10
for j = 1:4
    
    position_ekf{j} = zeros(length(samples{j}), 6);
    position_nls{j} = zeros(length(samples{j}), 2);

    
    signal_dummy = filtered{j};
   
    R = buildCovarianceMatrix(signal_dummy,numel(signal_dummy(1,:)));
   

    for i = 1:length(samples{j})
        
        [ uHat{j} , numberOfPerformedIterations{j}] = iterativeNLS(NiterMax, AP(:,1:2),filtered{j}(:,i),i ,R); %this function should return a NiterMax x 2 vector of the estimates provided by the NLS
        iterations{j}(i) = numberOfPerformedIterations{j};
        position_nls{j}(i,:) = uHat{j}(numberOfPerformedIterations{j}, :);
        
        [ uHatg{j} , numberOfPerformedIterationsg{j}] = iterativeNLS(NiterMax, AP(:,1:2),smoothData{j}(:,i),i ,R); %this function should return a NiterMax x 2 vector of the estimates provided by the NLS
        iterationsg{j}(i) = numberOfPerformedIterationsg{j};
        position_groundroot{j}(i,:) = uHatg{j}(numberOfPerformedIterationsg{j}, :);
       
    end



    position_ekf{j} = KFTrack(AP(:,1:2),filtered{j},length(samples{j}),position_groundroot{j}(1,:),R,0.08);%sigma_driving(4)
    
%     if j == 1
%         error_dum = sqrt(sum((position_groundroot{j}(1:200,1:2) - position_ekf{j}(1:200,1:2)).^2,2));
%         xi1 = linspace(1,200,500);
%         error{k}(j,:) = interp1([1:200],error_dum,xi1);
% 
%     else
%         error{k}(j,:) = sqrt(sum((position_groundroot{j}(1:500,1:2) - position_ekf{j}(1:500,1:2)).^2,2));
%     end


    
end

% error1(1,:) = error{k}(1,:);
% error1(2,:) = error{k}(2,:);
% error1(3,:) = error{k}(3,:);
% error1(4,:) = error{k}(4,:);
% 
% 
% 
% 
% 
% averageError = (error1(1,:)+ error1(2,:) +error1(3,:) + error1(4,:))/4;
% cdfplot(averageError);
% hold on;
% 
% 
% 
% end
% legend("Sigma Driving 0.01m","Sigma Driving 0.02m","Sigma Driving  0.03m","Sigma Driving  0.04m","Sigma Driving 0.05m",...
%     "Sigma Driving  0.06m","Sigma Driving  0.07m","Sigma Driving  0.08m","Sigma Driving 0.09m","Sigma Driving 0.1m");
% 



figure;
plot(position_nls{1}(:,1),position_nls{1}(:,2),'b');
hold on
plot(position_ekf{1}(:,1),position_ekf{1}(:,2),'b--');
title('Tag 1');
legend('NLS', 'EKF');
grid on;
axis square;
xlim([5 12]);
ylim([16 18.5]);
xlabel(" X (m)");
ylabel(" Y (m)");


figure;
plot(position_nls{2}(:,1),position_nls{2}(:,2),'r');
hold on
plot(position_ekf{2}(:,1),position_ekf{2}(:,2),'r--');
title('Tag 2');
legend('NLS', 'EKF');
grid on;
axis square;
xlim([5 12]);
ylim([16 18.5]);
xlabel(" X (m)");
ylabel(" Y (m)");

figure;
plot(position_nls{3}(:,1),position_nls{3}(:,2),'c');
hold on
plot(position_ekf{3}(:,1),position_ekf{3}(:,2),'c--');
title('Tag 3');
legend('NLS', 'EKF');
grid on;
axis square;
xlim([5 12]);
ylim([16 18.5]);
xlabel(" X (m)");
ylabel(" Y (m)");


figure;
plot(position_nls{4}(:,1),position_nls{4}(:,2),'g');
hold on
plot(position_ekf{4}(:,1),position_ekf{4}(:,2),'g--');
title('Tag 4');
legend('NLS', 'EKF');
grid on;
axis square;
xlim([5 12]);
ylim([16 18.5]);
xlabel(" X (m)");
ylabel(" Y (m)");


%% Noise 
noiseSamSum = 800 + 2008 + 2006 + 2014;
noisew1 = 800 / noiseSamSum;
noisew2 = 2008 / noiseSamSum;
noisew3 = 2006 / noiseSamSum;
noisew4 = 2014 / noiseSamSum;

for i = 1:4

noiseCal{i} = mean(abs(rhoIm{i} - smoothData{i}),1);

noiseMean(i) = mean(noiseCal{i});

noiseVar(i) = var(noiseCal{i});

end

meanNoiseVar = sqrt(noisew1 * noiseVar(1) + noisew2 * noiseVar(2)...
             + noisew3 * noiseVar(3) + noisew4 * noiseVar(4)) ;
%% Average Position 

pos_ekf1(:,1) = position_ekf{1}(:,1);
pos_ekf1(:,2) = position_ekf{1}(:,2);
pos_ekf2(:,1) = position_ekf{2}(:,1);
pos_ekf2(:,2) = position_ekf{2}(:,2);
pos_ekf3(:,1) = position_ekf{3}(:,1);
pos_ekf3(:,2) = position_ekf{3}(:,2);
pos_ekf4(:,1) = position_ekf{4}(:,1);
pos_ekf4(:,2) = position_ekf{4}(:,2);

xi1 = linspace(1,numel(samples{1}),numel(samples{4}));
xi2 = linspace(1,numel(samples{2}),numel(samples{4}));
xi3 = linspace(1,numel(samples{3}),numel(samples{4}));

pos_ekf1p(:,1) = interp1(samples{1},pos_ekf1(:,1),xi1);
pos_ekf1p(:,2) = interp1(samples{1},pos_ekf1(:,2),xi1);
pos_ekf2p(:,1) = interp1(samples{2},pos_ekf2(:,1),xi2);
pos_ekf2p(:,2) = interp1(samples{2},pos_ekf2(:,2),xi2);
pos_ekf3p(:,1) = interp1(samples{3},pos_ekf3(:,1),xi3);
pos_ekf3p(:,2) = interp1(samples{3},pos_ekf3(:,2),xi3);


averageEKF = (pos_ekf1p + pos_ekf2p + pos_ekf3p + pos_ekf4)/4;

pos_nls1(:,1) = position_nls{1}(:,1);
pos_nls1(:,2) = position_nls{1}(:,2);
pos_nls2(:,1) = position_nls{2}(:,1);
pos_nls2(:,2) = position_nls{2}(:,2);
pos_nls3(:,1) = position_nls{3}(:,1);
pos_nls3(:,2) = position_nls{3}(:,2);
pos_nls4(:,1) = position_nls{4}(:,1);
pos_nls4(:,2) = position_nls{4}(:,2);

pos_grn1(:,1) = position_groundroot{1}(:,1);
pos_grn1(:,2) = position_groundroot{1}(:,2);
pos_grn2(:,1) = position_groundroot{2}(:,1);
pos_grn2(:,2) = position_groundroot{2}(:,2);
pos_grn3(:,1) = position_groundroot{3}(:,1);
pos_grn3(:,2) = position_groundroot{3}(:,2);
pos_grn4(:,1) = position_groundroot{4}(:,1);
pos_grn4(:,2) = position_groundroot{4}(:,2);


pos_nls1p(:,1) = interp1(samples{1},pos_nls1(:,1),xi1);
pos_nls1p(:,2) = interp1(samples{1},pos_nls1(:,2),xi1);
pos_nls2p(:,1) = interp1(samples{2},pos_nls2(:,1),xi2);
pos_nls2p(:,2) = interp1(samples{2},pos_nls2(:,2),xi2);
pos_nls3p(:,1) = interp1(samples{3},pos_nls3(:,1),xi3);
pos_nls3p(:,2) = interp1(samples{3},pos_nls3(:,2),xi3);

pos_grn1p(:,1) = interp1(samples{1},pos_grn1(:,1),xi1);
pos_grn1p(:,2) = interp1(samples{1},pos_grn1(:,2),xi1);
pos_grn2p(:,1) = interp1(samples{2},pos_grn2(:,1),xi2);
pos_grn2p(:,2) = interp1(samples{2},pos_grn2(:,2),xi2);
pos_grn3p(:,1) = interp1(samples{3},pos_grn3(:,1),xi3);
pos_grn3p(:,2) = interp1(samples{3},pos_grn3(:,2),xi3);


averageNLS = (pos_nls1p + pos_nls2p + pos_nls3p + pos_nls4)/4;
averageGRN = (pos_grn1p + pos_grn2p + pos_grn3p + pos_grn4)/4;


plot(averageNLS(:,1),averageNLS(:,2),'LineWidth',1);
hold on
plot(averageEKF(:,1),averageEKF(:,2),'LineWidth',1);
plot(averageGRN(:,1),averageGRN(:,2),'LineWidth',1);
legend('NLS','EKF','Ground Root');
title('Average Positioning');
grid on;

%% Errors

errorEKFx = abs(averageGRN(1:500,1) - averageEKF(1:500,1)); %MSE EKF
errorEKFy = abs(averageGRN(1:500,2) - averageEKF(1:500,2));
errorNLSx = abs(averageGRN(1:500,1) - averageNLS(1:500,1)); %MSE NLS
errorNLSy = abs(averageGRN(1:500,2) - averageNLS(1:500,2));

errorEKFxs = abs(averageGRN(501:2014,1) - averageEKF(501:2014,1)); %MSE EKF
errorEKFys = abs(averageGRN(501:2014,2) - averageEKF(501:2014,2));
errorNLSxs = abs(averageGRN(501:2014,1) - averageNLS(501:2014,1)); %MSE NLS
errorNLSys = abs(averageGRN(501:2014,2) - averageNLS(501:2014,2));

deltaEKFx = mean(errorEKFx);
deltaEKFy = mean(errorEKFy);
deltaNLSx = mean(errorNLSx);
deltaNLSy = mean(errorNLSy);

deltaEKFxs = mean(errorEKFxs);
deltaEKFys = mean(errorEKFys);
deltaNLSxs = mean(errorNLSxs);
deltaNLSys = mean(errorNLSys);

drmsEKF = norm([deltaEKFx deltaEKFy]);
drmsNLS = norm([deltaNLSx deltaNLSy]);

errorEKF = [errorEKFx errorEKFy];
Cekf = errorEKF' * errorEKF ./2014;

errorNLS = [errorNLSx errorNLSy];
Cnls = errorNLS' * errorNLS ./2014;

Gekf = Cekf ./ meanNoiseVar;
Gnls = Cnls ./ meanNoiseVar;

DxEkf = sqrt(Gekf(1,1));
DyEkf = sqrt(Gekf(2,2));

DxNls = sqrt(Gnls(1,1));
DyNls = sqrt(Gnls(2,2));



%% Moving Plot
f = figure(1);
ax = axes;
hold on;
h = animatedline;
p = plot(ax,averageEKF(1,1),averageEKF(1,2),'s','LineStyle','none');
xlim([5 11]);
ylim([16.5 18.5]);
grid on;
%axis square;

samplenum = numel(samples{4});


filename = 'testnew51.gif';

for i = 1:samplenum
      
    p.XData = averageEKF(i,1);
    p.YData = averageEKF(i,2);
   
    
addpoints(h,averageEKF(i,1),averageEKF(i,2));
pp =text(0.01,0.01, strcat(num2str(i*0.1),' sec'), ...
     'Units', 'normalized', ...   % Not depending on the data
     'HorizontalAlignment', 'left', ...
     'VerticalAlignment', 'bottom');


frame = getframe(1);
      im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);
      if i == 1
          imwrite(imind,cm,filename,'gif','DelayTime',0.001, 'Loopcount',inf);
      else
          imwrite(imind,cm,filename,'gif','DelayTime',0.001,'WriteMode','append');
      end
    
   pause(0.02);

   if i~= 2014
        delete(pp);
   end
     
   
end
%% 

plot(rhoIm{2}(5,:),'LineWidth',1);
hold on
plot(smoothData{2}(5,:),'LineWidth',1);
xlim([1 2010]);
title('Original and Noise Removed Data');
legend('Original Signal (Impainted)' , 'Noise Removed Signal');


%% Position Errors Point by Point

figure;
subplot(4,1,1);
plot(errorEKFx);
title('Error EKF x');

subplot(4,1,3);
plot(errorEKFy);
title('Error EKF y');

subplot(4,1,2);
plot(errorNLSx);
title('Error NLS x');

subplot(4,1,4);
plot(errorNLSy);
title('Error NLS y');

%% 


errorEKFhist = sqrt(errorEKFx.^2 + errorEKFy.^2);
errorNLShist = sqrt(errorNLSx.^2 + errorNLSy.^2);

errorEKFhists = sqrt(errorEKFxs.^2 + errorEKFys.^2);
errorNLShists = sqrt(errorNLSxs.^2 + errorNLSys.^2);

figure;
cdfplot(errorEKFhist);
hold on
cdfplot(errorNLShist);
legend('EKF Horizontal Error','NLS Horizontal Error');
title('CDF of Horizontal Errors During Motion');

figure;
cdfplot(errorEKFhists);
hold on
cdfplot(errorNLShists);
legend('EKF Horizontal Error','NLS Horizontal Error');
title('CDF of Horizontal Errors During Static Phase');

figure;
histogram(errorEKFhist);
hold on;
histogram(errorNLShist);
legend('EKF Horizontal Error' , 'NLS Horizontal Error');
title('Histogram of Horizontal Errors');

figure;
histogram(errorEKFhists);
hold on;
histogram(errorNLShists);
legend('EKF Horizontal Error' , 'NLS Horizontal Error');
title('Histogram of Horizontal Errors During Static Phase');






