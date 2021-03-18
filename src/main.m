clc
clearvars -except gps_observations sbas_observations alpha beta a gps_navigation sbas_navigation
format longg
tic()
%% ---------- Choose a dataset ---------- %%
dataset = 1;
    
if dataset == 1 %RF2
    
    %Set the initial position estimate
    initialEstimate = llh2ecef([38.7166700,-9.1333300,0]); %Lisbon
    
    %Load pre-processed .mat file containing observations and navigation
    %messages
    load('data\7 day dataset\24h Splits\COM3_R2_RF2_25102020.mat');
    
    %Define the hatch filter smoothing constant
    hatchFilterParameter = 0.01;
    
    %Define the reference receiver position
    positionRef = [4918525.5233,-791212.0300,3969762.2262]; %RF2/COM3 - PPP 7 dias
    
    
elseif dataset == 2 %RF6
    
    %Set the initial position estimate
    initialEstimate = llh2ecef([38.7166700,-9.1333300,0]); %Lisbon
    
    %Load pre-processed .mat file containing observations and navigation
    %messages
    load('data\7 day dataset\24h Splits\COM4_R1_RF6_25102020.mat');
    
    %Define the hatch filter smoothing constant
    hatchFilterParameter = 0.01;
    
    %Define the reference receiver position
    positionRef = [4918532.1188,-791212.5264,3969754.7230]; %RF6/COM4 - PPP 7 dias 
    
elseif dataset == 5 %24h dataset at 1Hz from IGS([-34.6,-58.367,0] = [2756517.813191,-4474878.70094808,-3601428.67912596]
    
    %Set the initial position estimate
    initialEstimate = [2756517.813191,-4474878.70094808,-3601428.67912596]; %Buenos Aires, [-34.6,-58.367,0]
    
    %Load pre-processed .mat file containing observations and navigation
    %messages
    load('data\RINEX examples\MGUE00ARG_R_20202990000_01D_01S_MO.mat');
    %Define the hatch filter smoothing constant
    hatchFilterParameter = 0.01;
    
    %Define the reference receiver position
    positionRef = [1823327.9138,-4850352.4421,-3709085.5307];
end


%%
% ---------- Set the initial conditions ---------- %

%Set the mask angle for 10�
maskAngle = 10;

% Set the EKF parameters
h_0 = 2*10^-19; %Temperature-compensated crystal
h_m2 = 2*10^-20; %Temperature-compensated crystal
sigmaR = 10;
P_initial = [10^4,10^4,10^4,9*10^6,9*10^3];    

EKFStart = 1; %Number of iterations after which EKF is considered usable
UKFStart = 1; %Number of iterations after which UKF is considered usable

%Set the required Pfa for RAIM and set the RAIM test flag
Pfa = 8*10^-6;
enableRAIMTest = 1;

%---------- Initialize variables ---------- %

%Create an array to detect cycle slips and to store the Hatch filter
%weights for each satellite
hatchParameters = [ones(32,1)];

%Create an array to store the previous observations
observationPrev = zeros(32,2); %Column 1: last pseudorange; Column 2: last phase

navMessageCurrent = zeros(32,31);%Create an array to save the most current navigation message for all 32 satellites.
navMessageIndex = 1; %This variable is used as a counter for the current line of the parsed navigation messages

phaseHistory = zeros(32,1);                                 %DEBUG
pseudorangeHistory = zeros(32,gps_observations(end,3)+1);   %DEBUG

%Get all the RAIM thresholds
RAIMThreshold = RAIM('initialize',Pfa);


%%
% ---------- Run the simulator ---------- %


%Iterate for each observation epoch, t
j = 1; %Counter for the number of iterations
skipped = 0; %DEBUG

%Copy the initial estimates
initialEstimate_LS = initialEstimate;
initialEstimate_WLS = initialEstimate;
initialEstimate_EKF = initialEstimate;
initialEstimate_UKF = initialEstimate;


for i=1:16:size(gps_observations,1)
    t0 = gps_observations(1,3); %First observation defines the starting time
    
    if i == 1
        dt = gps_observations(17,3) - t0; %Get the time step in seconds
    elseif i <= size(gps_observations,1)-16
        dt = gps_observations(i+16,3) - gps_observations(i,3);
    end
    
    t = gps_observations(i,3); %All observations have the same epoch, just use the first one
    WN_LSF = gps_observations(i,2);
    
    if t == 18737       %DEBUG
        t = 18737;      %DEBUG
    end                 %DEBUG
    
    if dt > 1
        fprintf("Error! dt > 1!!\n");
        beep
        beep
        dt = dt;
    end
    
    %Get a slice of the GPS observations array corresponding to the current
    %receiver epoch
    gpsObservationCurrent = gps_observations(i:i+15,:);
    if t~=t0
        observationOld = observationFiltered;   %For the Hatch Filter
    end
    
    %Update the navigation message if there's a new message when
    %the current epoch equals the time of transmission of the navigation
    %message
    [navMessageCurrent,navMessageIndex] = updateNavMessage(gps_navigation,navMessageCurrent,navMessageIndex,t,WN_LSF);
    
    
    %Couple GPS observations with GPS navigation messages
    [navMessageFiltered,observationFiltered] = coupleNavObs(t,WN_LSF,gpsObservationCurrent,navMessageCurrent);
    
    %Apply the elevation mask
    if t ~= t0
        [observationFiltered,navMessageFiltered] = elevationMask(t,maskAngle,receiverPos_LS(j-1,1:3),observationFiltered,navMessageFiltered);
    end
    
    
    %Only continue if there's enough satellites
    if size(observationFiltered) >= 4
        
        if t == 10970       %DEBUG
            t = 10970;
        end
        
        if enableRAIMTest == 1
                %Add test biases to check RAIM performance
                if t >= 100 && t < 200      
                    ii = 3;
                    observationFiltered(ii,5) = observationFiltered(ii,5) + 15;
                end                 
        
                if t >= 3000 && t < 3500       
                    ii = 2;
                    observationFiltered(ii,5) = observationFiltered(ii,5) + 30;
                end                 
        
                
                if t > 10000 && t < 12000       
                    ii = 2;
                    observationFiltered(ii,5) = observationFiltered(ii,5) + (t-10000);
                end                 
        end
        
        %----Remove all the modeled errors----%
        [observationFiltered,satPos_tTX,tTX] = removeErrors(observationFiltered,navMessageFiltered,t,alpha,beta,initialEstimate_LS);
        
        for k=1:size(observationFiltered,1)
            currentSVN = observationFiltered(k,4);
            phaseHistory(currentSVN,j) = observationFiltered(k,6);
            pseudorangeHistory(currentSVN,j) = observationFiltered(k,5);
        end
        
        %Save the total number and PRN of visible satellites
        nSatsVisible(j) = size(navMessageFiltered,1);
        obsSatUsed(t-t0+1,:) = zeros(1,32);                     %DEBUG
        navSatUsed(t-t0+1,:) = zeros(1,32);                     %DEBUG
        for kk = 1:size(observationFiltered,1)
            currentSVN = observationFiltered(kk,4);
            obsSatUsed(t-t0+1,currentSVN) = currentSVN;
        end
        
        for kk = 1:size(navMessageFiltered,1)                   %DEBUG
            currentSVN = navMessageFiltered(kk,1);              %DEBUG
            navSatUsed(t-t0+1,currentSVN) = currentSVN;         %DEBUG
        end                                                     %DEBUG
        
        if (sum(obsSatUsed(t-t0+1,:)-navSatUsed(t-t0+1,:)) ~=0) %DEBUG
            disp ("nav and obs satellite mismatch!!!");         %DEBUG
        end                                                     %DEBUG
        
        
        %----Implement the Hatch Filter----%
        [hatchParameters] = detectCycleSlip(observationFiltered,hatchParameters);
        
        [observationFiltered(:,5),hatchParameters,observationPrev] = filterHatch(observationFiltered,hatchParameters,observationPrev,hatchFilterParameter);
        
        
        if t~= t0
            [RAIMLevel(j),observationFiltered,navMessageFiltered,satPos_tTX,testStatistic(j)] = RAIM('run',observationFiltered,navMessageFiltered,satPos_tTX,RAIMThreshold,initialEstimate_LS);
            %RAIMLevel(j) = 0; %DEBUG
        else
            RAIMLevel(j) = 0; %To allow first run
        end
        
        %Purge the epoch if the RAIM algorithm detected more than 1
        %problematic satellite
        if RAIMLevel(j) ~= -1
            %----Least Squares method----%
            [receiverPos_LS(j,:),residuals_LS,GDOP_LS(j),PDOP_LS(j),TDOP_LS(j),H_LS] = leastSquares(observationFiltered,satPos_tTX,initialEstimate_LS);
            
            initialEstimate_LS(1:3) = receiverPos_LS(j,1:3);
            
            %----Weighted Least Squares method----%
            [receiverPos_WLS(j,:),residuals_WLS,GDOP_WLS(j),PDOP_WLS(j),TDOP_WLS(j),H_WLS] = weightedLeastSquares(observationFiltered,navMessageFiltered,satPos_tTX,initialEstimate_LS);
            initialEstimate_WLS(1:3) = receiverPos_WLS(j,1:3);
            
            %----Extended Kalman Filter method----%
            if t == t0
                %Initialize the Extended Kalman Filter
                initialEstimate_EKF = receiverPos_LS(j,1:3); %Initialize the EKF by running the Least Squares method once to improve convergence speed
                [PHI,Qk_EKF,Hk_EKF,Pk_m_EKF,Xk_est_EKF,Xk_pred_EKF,Pk_pred_EKF,Kk,INOV] = extendedKF('initialize',h_0,h_m2,sigmaR,P_initial,observationFiltered,satPos_tTX,initialEstimate_EKF,dt);
            else
                [PHI,Qk_EKF,Hk_EKF,Pk_m_EKF,Xk_est_EKF,Xk_pred_EKF,Pk_pred_EKF,Kk,INOV] = extendedKF('recursive',h_0,h_m2,sigmaR,P_initial,observationFiltered,satPos_tTX,initialEstimate_EKF,dt,PHI,Qk_EKF,Xk_pred_EKF,Pk_pred_EKF);
            end
            receiverPos_EKF(j,:) = Xk_est_EKF(1:3);
            
            %----Unscented Kalman Filter method----%
            if t == t0
                %Initialize the Unscented Kalman Filter
                initialEstimate_UKF = receiverPos_LS(j,1:3); %Initialize the UKF by running the Least Squares method once to improve convergence speed
                [F_UKF,Qk_UKF,Xk_est_UKF,Pk_est_UKF] = unscentedKF('initialize',h_0,h_m2,sigmaR,P_initial,observationFiltered,satPos_tTX,initialEstimate_EKF,dt);
            else
                [F_UKF,Qk_UKF,Xk_est_UKF,Pk_est_UKF] = unscentedKF('recursive',h_0,h_m2,sigmaR,P_initial,observationFiltered,satPos_tTX,initialEstimate_EKF,dt,F_UKF,Qk_UKF,Xk_est_UKF,Pk_est_UKF);
            end
            receiverPos_UKF(j,:) = Xk_est_UKF(1:3);
            
            
            positionErrorLS(j) = norm(positionRef - receiverPos_LS(j,1:3));
            errorLS = norm(positionErrorLS(j));
            %fprintf("Position error LS: %f meters\n",errorLS);
            
            positionErrorWLS(j) = norm(positionRef - receiverPos_WLS(j,1:3));
            errorWLS = norm(positionErrorWLS(j));
            %fprintf("Position error LS: %f meters\n",errorLS);
            
            positionErrorEKF(j) = norm(positionRef - receiverPos_EKF(j,:));
            errorEKF = norm(positionErrorEKF(j));
            %fprintf("Position error EFK: %f meters\n",positionErrorEKF(j-1));
            
            positionErrorUKF(j) = norm(positionRef - receiverPos_UKF(j,:));
            errorEKF = norm(positionErrorEKF(j));
            %fprintf("Position error UFK: %f meters\n",positionErrorUKF(j-1));
            
            %Time series
            receiverPos_LS_time(t-t0+1,:) = receiverPos_LS(j,1:3);
            positionErrorLS_time(t-t0+1) = positionErrorLS(j);
            
            receiverPos_WLS_time(t-t0+1,:) = receiverPos_WLS(j,1:3);
            positionErrorWLS_time(t-t0+1) = positionErrorWLS(j);
            
            receiverPos_EKF_time(t-t0+1,:) = receiverPos_EKF(j,1:3);
            positionErrorEKF_time(t-t0+1) = positionErrorEKF(j);
            
            receiverPos_UKF_time(t-t0+1,:) = receiverPos_UKF(j,1:3);
            positionErrorUKF_time(t-t0+1) = positionErrorUKF(j);
            
            
            j=j+1;
        else
            skipped = skipped + 1;  %DEBUG
            fprintf("Skipped epoch %d due to faulty measurements (RAIM)\n",t);
            
            %Time series
            receiverPos_LS_time(t-t0+1,:) = NaN;
            positionErrorLS_time(t-t0+1) = NaN;
            
            receiverPos_WLS_time(t-t0+1,:) = NaN;
            positionErrorWLS_time(t-t0+1) = NaN;
            
            receiverPos_EKF_time(t-t0+1,:) = NaN;
            positionErrorEKF_time(t-t0+1) = NaN;
            
            receiverPos_UKF_time(t-t0+1,:) = NaN;
            positionErrorUKF_time(t-t0+1) = NaN;
            
            RAIMLevel(t-t0+1) = NaN;
            testStatistic(t-t0+1) = NaN;
            
        end
    else                        %DEBUG
        for kk = 1:size(observationFiltered,1)
            currentSVN = observationFiltered(kk,4);
            obsSatUsed(t-t0+1,currentSVN) = currentSVN;
        end
        skipped = skipped + 1;  %DEBUG
        fprintf("Skipped epoch %d due to insufficient number of satellites: n = %d\n",t,size(observationFiltered,1));
        
        %Time series
        receiverPos_LS_time(t-t0+1,:) = NaN;
        positionErrorLS_time(t-t0+1) = NaN;
        
        receiverPos_WLS_time(t-t0+1,:) = NaN;
        positionErrorWLS_time(t-t0+1) = NaN;
        
        receiverPos_EKF_time(t-t0+1,:) = NaN;
        positionErrorEKF_time(t-t0+1) = NaN;
        
        receiverPos_UKF_time(t-t0+1,:) = NaN;
        positionErrorUKF_time(t-t0+1) = NaN;
        
        RAIMLevel(t-t0+1) = NaN;
        testStatistic(t-t0+1) = NaN;
        
        
    end                         
    
    
    
    
    
    
    
end
%%
%-----------Error analysis-----------%
%clc
fprintf("========== Error Analysis ==========\n");

%-----Least Squares error analysis-----%
fprintf("----- Least Squares-----\n");
positionAvgLS(1) = sum(receiverPos_LS(:,1))/size(receiverPos_LS,1);
positionAvgLS(2) = sum(receiverPos_LS(:,2))/size(receiverPos_LS,1);
positionAvgLS(3) = sum(receiverPos_LS(:,3))/size(receiverPos_LS,1);
errorLS = norm(positionRef - positionAvgLS);
fprintf("The least-squares error is: %f meters\n",errorLS);

positionRMS_LS = sqrt(sum((receiverPos_LS(:,1)-positionRef(1)).^2 + (receiverPos_LS(:,2)-positionRef(2)).^2 + (receiverPos_LS(:,3)-positionRef(3)).^2)/size(receiverPos_LS,1));
fprintf("The least-squares RMS error is: %f meters\n",positionRMS_LS);

positionDRMS_LS  = accMetrics2d('drms',receiverPos_LS(:,1:3));
positionCEP_LS   = accMetrics2d('cep',receiverPos_LS(:,1:3));
positionR95_LS   = accMetrics2d('r95',receiverPos_LS(:,1:3));
positionMRSE_LS  = accMetrics3d('mrse',receiverPos_LS(:,1:3));
positionSEP_LS   = accMetrics3d('sep',receiverPos_LS(:,1:3));
positionSAS90_LS = accMetrics3d('sas90',receiverPos_LS(:,1:3));

fprintf("The least-squares DRMS error is: %f meters\n",positionDRMS_LS);
fprintf("The least-squares CEP error is: %f meters\n",positionCEP_LS);
fprintf("The least-squares R95 error is: %f meters\n",positionR95_LS);
fprintf("The least-squares MRSE error is: %f meters\n",positionMRSE_LS);
fprintf("The least-squares SEP error is: %f meters\n",positionSEP_LS);
fprintf("The least-squares SAS90 error is: %f meters\n\n",positionSAS90_LS);


%-----Weighted Least Squares error analysis-----%
fprintf("----- Weighted Least Squares-----\n");
positionAvgWLS(1) = sum(receiverPos_WLS(:,1))/size(receiverPos_WLS,1);
positionAvgWLS(2) = sum(receiverPos_WLS(:,2))/size(receiverPos_WLS,1);
positionAvgWLS(3) = sum(receiverPos_WLS(:,3))/size(receiverPos_WLS,1);
errorWLS = norm(positionRef - positionAvgWLS);
fprintf("The weighted least-squares error is: %f meters\n",errorWLS);

positionRMS_WLS = sqrt(sum((receiverPos_WLS(:,1)-positionRef(1)).^2 + (receiverPos_WLS(:,2)-positionRef(2)).^2 + (receiverPos_WLS(:,3)-positionRef(3)).^2)/size(receiverPos_WLS,1));
fprintf("The weighted least-squares RMS error is: %f meters\n",positionRMS_WLS);

positionDRMS_WLS  = accMetrics2d('drms',receiverPos_WLS(:,1:3));
positionCEP_WLS   = accMetrics2d('cep',receiverPos_WLS(:,1:3));
positionR95_WLS   = accMetrics2d('r95',receiverPos_WLS(:,1:3));
positionMRSE_WLS  = accMetrics3d('mrse',receiverPos_WLS(:,1:3));
positionSEP_WLS   = accMetrics3d('sep',receiverPos_WLS(:,1:3));
positionSAS90_WLS = accMetrics3d('sas90',receiverPos_WLS(:,1:3));

fprintf("The weighted least-squares DRMS error is: %f meters\n",positionDRMS_WLS);
fprintf("The weighted least-squares CEP error is: %f meters\n",positionCEP_WLS);
fprintf("The weighted least-squares R95 error is: %f meters\n",positionR95_WLS);
fprintf("The weighted least-squares MRSE error is: %f meters\n",positionMRSE_WLS);
fprintf("The weighted least-squares SEP error is: %f meters\n",positionSEP_WLS);
fprintf("The weighted least-squares SAS90 error is: %f meters\n\n",positionSAS90_WLS);


%-----Extended Kalman Filter error analysis-----%
fprintf("----- Extended Kalman Filter-----\n");
fprintf("Note: The first %d interations where ignored\n",EKFStart-1);
positionAvgEKF(1) = sum(receiverPos_EKF(EKFStart:end,1))/size(receiverPos_EKF(EKFStart:end,:),1);
positionAvgEKF(2) = sum(receiverPos_EKF(EKFStart:end,2))/size(receiverPos_EKF(EKFStart:end,:),1);
positionAvgEKF(3) = sum(receiverPos_EKF(EKFStart:end,3))/size(receiverPos_EKF(EKFStart:end,:),1);
errorEKF = norm(positionRef - positionAvgEKF);
fprintf("The Extended Kalman Filter error is: %f meters\n",errorEKF);

positionRMS_EKF = sqrt(sum((receiverPos_EKF(EKFStart:end,1)-positionRef(1)).^2 + (receiverPos_EKF(EKFStart:end,2)-positionRef(2)).^2 + (receiverPos_EKF(EKFStart:end,3)-positionRef(3)).^2)/size(receiverPos_EKF,1));
fprintf("The Extended Kalman Filter RMS error is: %f meters\n",positionRMS_EKF);

positionDRMS_EKF  = accMetrics2d('drms',receiverPos_EKF(:,1:3));
positionCEP_EKF   = accMetrics2d('cep',receiverPos_EKF(:,1:3));
positionR95_EKF   = accMetrics2d('r95',receiverPos_EKF(:,1:3));
positionMRSE_EKF  = accMetrics3d('mrse',receiverPos_EKF(:,1:3));
positionSEP_EKF   = accMetrics3d('sep',receiverPos_EKF(:,1:3));
positionSAS90_EKF = accMetrics3d('sas90',receiverPos_EKF(:,1:3));

fprintf("The Extended Kalman Filter DRMS error is: %f meters\n",positionDRMS_EKF);
fprintf("The Extended Kalman Filter CEP error is: %f meters\n",positionCEP_EKF);
fprintf("The Extended Kalman Filter R95 error is: %f meters\n",positionR95_EKF);
fprintf("The Extended Kalman Filter MRSE error is: %f meters\n",positionMRSE_EKF);
fprintf("The Extended Kalman Filter SEP error is: %f meters\n",positionSEP_EKF);
fprintf("The Extended Kalman Filter SAS90 error is: %f meters\n\n",positionSAS90_EKF);



%-----Unscented Kalman Filter error analysis-----%
fprintf("----- Unscented Kalman Filter-----\n");
fprintf("Note: The first %d interations where ignored\n",UKFStart-1);
positionAvgUKF(1) = sum(receiverPos_UKF(UKFStart:end,1))/size(receiverPos_UKF(UKFStart:end,:),1);
positionAvgUKF(2) = sum(receiverPos_UKF(UKFStart:end,2))/size(receiverPos_UKF(UKFStart:end,:),1);
positionAvgUKF(3) = sum(receiverPos_UKF(UKFStart:end,3))/size(receiverPos_UKF(UKFStart:end,:),1);
errorUKF = norm(positionRef - positionAvgUKF);
fprintf("The Unscented Kalman Filter error is: %f meters\n",errorUKF);

positionRMS_UKF = sqrt(sum((receiverPos_UKF(UKFStart:end,1)-positionRef(1)).^2 + (receiverPos_UKF(UKFStart:end,2)-positionRef(2)).^2 + (receiverPos_UKF(UKFStart:end,3)-positionRef(3)).^2)/size(receiverPos_UKF,1));
fprintf("The Unscented Kalman Filter RMS error is: %f meters\n",positionRMS_UKF);

positionDRMS_UKF  = accMetrics2d('drms',receiverPos_UKF(:,1:3));
positionCEP_UKF   = accMetrics2d('cep',receiverPos_UKF(:,1:3));
positionR95_UKF   = accMetrics2d('r95',receiverPos_UKF(:,1:3));
positionMRSE_UKF  = accMetrics3d('mrse',receiverPos_UKF(:,1:3));
positionSEP_UKF   = accMetrics3d('sep',receiverPos_UKF(:,1:3));
positionSAS90_UKF = accMetrics3d('sas90',receiverPos_UKF(:,1:3));

fprintf("The Unscented Kalman Filter DRMS error is: %f meters\n",positionDRMS_UKF);
fprintf("The Unscented Kalman Filter CEP error is: %f meters\n",positionCEP_UKF);
fprintf("The Unscented Kalman Filter R95 error is: %f meters\n",positionR95_UKF);
fprintf("The Unscented Kalman Filter MRSE error is: %f meters\n",positionMRSE_UKF);
fprintf("The Unscented Kalman Filter SEP error is: %f meters\n",positionSEP_UKF);
fprintf("The Unscented Kalman Filter SAS90 error is: %f meters\n\n",positionSAS90_UKF);

%
%%
%-----Error plots-----%
figure
tPlot = transpose(seconds(0:1:size(positionErrorLS_time,2)-1));
hold on
plot(tPlot,positionErrorLS_time);
hold on
plot(tPlot,positionErrorWLS_time);
hold on
plot(tPlot,positionErrorEKF_time);
hold on
plot(tPlot,positionErrorUKF_time);
%ylim([0 15])
%yticks(1:1:32)
xticks(seconds(0:7200:size(positionErrorLS_time,2)));
xlim([seconds(0) seconds(size(positionErrorLS_time,2))]);
xtickangle(45)
xtickformat('hh:mm:ss')
legend('Least-Squares','Weighted Least-Squares','Extended Kalman Filter','Unscented Kalman Filter');
ylabel('Absolute error (m)')

%%
%-----XYZ Plots-----%
figure
hold on
title('Error in X');
plot(receiverPos_LS(:,1)-positionRef(1));
plot(receiverPos_WLS(:,1)-positionRef(1));
plot(receiverPos_EKF(:,1)-positionRef(1));
plot(receiverPos_UKF(:,1)-positionRef(1));
legend('Least-Squares','Weighted Least-Squares','Extended Kalman Filter','Unscented Kalman Filter');

figure
hold on
title('Error in Y');
plot(receiverPos_LS(:,2)-positionRef(2));
plot(receiverPos_WLS(:,2)-positionRef(2));
plot(receiverPos_EKF(:,2)-positionRef(2));
plot(receiverPos_UKF(:,2)-positionRef(2));
legend('Least-Squares','Weighted Least-Squares','Extended Kalman Filter','Unscented Kalman Filter');

figure
hold on
title('Error in Z');
plot(receiverPos_LS(:,3)-positionRef(3));
plot(receiverPos_WLS(:,3)-positionRef(3));
plot(receiverPos_EKF(:,3)-positionRef(3));
plot(receiverPos_UKF(:,3)-positionRef(3));
legend('Least-Squares','Weighted Least-Squares','Extended Kalman Filter','Unscented Kalman Filter');

%% Pseudorange and Phase plots, for debugging only
% figure
% for i=1:32
%     plot(pseudorangeHistory(i,:),'+','LineStyle','none','MarkerSize',3)
%     title(['Pseudorange for SVN',num2str(i)]);
%     pause
% end
% figure
% title('Phase observations');
% plot(1:j-1,phaseHistory(:,:));
% for i=1:32
%     legendStr{i} = ['SVN' num2str(i)];
% end
% legend(legendStr);
% plotyy(1:j-1,positionErrorLS,1:j-1,phaseHistory(:,:));
% legend('Position Error Least-Squares','Position Error Extended Kalman Filter','Position Error Unscented Kalman Filter');

%% Compute the error in ENU components
positionRef_llh = ecef2llh(positionRef);
for i = 1:size(receiverPos_LS,1)
    errorENU_LS(i,:) = ecef2enu(positionRef,receiverPos_LS(i,1:3),positionRef_llh(1),positionRef_llh(2));
    errorENU_WLS(i,:) = ecef2enu(positionRef,receiverPos_WLS(i,1:3),positionRef_llh(1),positionRef_llh(2));
    errorENU_EKF(i,:) = ecef2enu(positionRef,receiverPos_EKF(i,1:3),positionRef_llh(1),positionRef_llh(2));
    errorENU_UKF(i,:) = ecef2enu(positionRef,receiverPos_UKF(i,1:3),positionRef_llh(1),positionRef_llh(2));
end

figure
hold on
title('N-S Error');
plot(errorENU_LS(:,2));
plot(errorENU_WLS(:,2));
plot(errorENU_EKF(:,2));
plot(errorENU_UKF(:,2));
legend('Least-Squares','Weighted Least-Squares','Extended Kalman Filter','Unscented Kalman Filter');

figure
hold on
title('E-W Error');
plot(errorENU_LS(:,1));
plot(errorENU_WLS(:,1));
plot(errorENU_EKF(:,1));
plot(errorENU_UKF(:,1));
legend('Least-Squares','Weighted Least-Squares','Extended Kalman Filter','Unscented Kalman Filter');

figure
hold on
title('U-D Error');
plot(errorENU_LS(:,3));
plot(errorENU_WLS(:,3));
plot(errorENU_EKF(:,3));
plot(errorENU_UKF(:,3));
legend('Least-Squares','Weighted Least-Squares','Extended Kalman Filter','Unscented Kalman Filter');

%%
figure;
tPlot = transpose(seconds(0:1:size(positionErrorLS_time,2)-1));
%title('DOP values');
%plot(tPlot,GDOP_LS,tPlot,PDOP_LS,tPlot,TDOP_LS);
plot(tPlot,GDOP_LS);
xticks(seconds(0:7200:size(positionErrorLS_time,2)));
xlim([seconds(0) seconds(size(positionErrorLS_time,2))]);
xtickangle(45)
xtickformat('hh:mm:ss')
%legend('GDOP','PDOP','TDOP');

%% Plot the RAIM statistics and threshold

for kk = 1:size(testStatistic,2)
    n = nnz(obsSatUsed(kk,:));
    
    if n ~= 0
        threshold(kk) = RAIMThreshold(n);
    else
        threshold(kk) = NaN;
    end
end
plot(testStatistic,'+','LineStyle','none','MarkerSize',5)
hold on
plot(threshold)


%% Plot the visible satellites chart
tPlot = transpose(seconds(0:1:size(obsSatUsed,1)-1));
for i=1:size(obsSatUsed,1)
    for ii=1:32
        if obsSatUsed(i,ii) == 0
            obsSatUsed(i,ii) = NaN;
        end
    end
end
figure;
%subplot(2,1,1);
plot(tPlot,obsSatUsed,'+','LineStyle','none','MarkerSize',3);
yticks(1:1:32)
ylim([0 33])
xticks(seconds(0:7200:size(obsSatUsed,1)));
xlim([seconds(0) seconds(size(obsSatUsed,1))]);
xtickangle(45)
xtickformat('hh:mm:ss')
grid on
%ax=gca;
%yticks(1:1:32)
%set(gca,'xtick',0:10000:gca.XAxis.Limits(2))
%set(gca,'XRuler.Exponent',0)
figure
%subplot(2,1,2);
plot(tPlot,nSatsVisible);
xticks(seconds(0:7200:size(obsSatUsed,1)));
xlim([seconds(0) seconds(size(obsSatUsed,1))]);
xtickangle(45)
xtickformat('hh:mm:ss')
%ylabel('Number of satellites in use')


toc()
%% Finished sound
audio = load('train');
sound(audio.y,audio.Fs)
