%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation of glucose dynamics parameters per subject
% Using the Least Squares Method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Cleaning
clear all; close all; clc;

%% Read data
dataTable = readtable('Dataset_glucosio_finale.csv');
    
%% Model parameters and estimation
v = 5;         % Delay for insulin (minutes)
w = 5;         % Delay for carbohydrates (minutes)
hist_len = 6;  % Length of history for beta and gamma

try
    uniquePatients = unique(dataTable.PatientID); 
catch ME
     error('Error: %s', ME.message);
end

numPatients = length(uniquePatients);
% Structure for storing results
estimatedParams = struct('PatientID', [], 'Age', [], 'AgeGroup', [], ...
                         'phi1', [], 'phi2', [], 'beta', [], 'gamma', []);
estimatedParams(numPatients).PatientID = []; % Preallocate the structure
estimationSuccess = false(numPatients, 1); % Flags to track successful estimates

%% Patient Estimation Loop
for i = 1:numPatients
    patientID = uniquePatients(i);
    
    % Data extraction for the current patient by managing different types of PatientID
    patientID_str = string(patientID); % Converting Current ID to String for Comparisons
        
    if iscell(dataTable.PatientID)
        % Table Column is Cell Array 
       patientData = dataTable(strcmp(dataTable.PatientID, patientID_str), :); 
    end

    % Time order
    patientData = sortrows(patientData, 'Minute');

    % Extract time series
    time = patientData.Minute;
    glucose = patientData.Glucose;
    insulin = patientData.Insulin;
    carbs = patientData.CarbsIntake;
    
    % Age management: take the earliest age found for the patient
    if ~isempty(patientData.Age) && isnumeric(patientData.Age) && ~isnan(patientData.Age(1))
       age = patientData.Age(1); 
    end
    
    n_obs = height(patientData); % Number of observations for this patient
    
    % Determine the starting index for the estimate
    start_idx = max([2, v + hist_len, w + hist_len]);
   
    % Construction of the matrix of regressors (X) and of the output vector (y)
    num_estimation_points = n_obs - start_idx + 1;
    y = glucose(start_idx : n_obs); 
    X = zeros(num_estimation_points, 2 + 2 * hist_len); 
    row_idx = 0; 
    valid_row = true(num_estimation_points, 1); % Trace valid lines
   
    for t_loop = start_idx : n_obs % Index for y(t_loop)
        t_reg = t_loop - 1; % Index for regressors G(t_reg), G(t_reg-1), etc.
        row_idx = row_idx + 1; 
        
        % Autoregressive terms
        if t_reg < 1 || (t_reg-1) < 1
            valid_row(row_idx) = false; continue; % Invalid indexes
        end
        X(row_idx, 1) = glucose(t_reg);      % G(t)
        X(row_idx, 2) = glucose(t_reg-1);    % G(t-1)
        
        % Insulin Terms (beta)
        for j = 1:hist_len
            idx_ins = t_reg - v - (j-1);
             if idx_ins > 0 && idx_ins <= n_obs % Check upper limit n_obs
                 X(row_idx, 2 + j) = insulin(idx_ins);
             else
                 valid_row(row_idx) = false; 
                 break; % Exit the inner loop j for this row
             end
        end
                
        % Carbohydrate terms (range)
         for j = 1:hist_len
            idx_carb = t_reg - w - (j-1);
             if idx_carb > 0 && idx_carb <= n_obs % Check upper limit n_obs
                X(row_idx, 2 + hist_len + j) = carbs(idx_carb);
             else
                 valid_row(row_idx) = false;
                 break; % Exit the inner loop j for this row
             end
         end
    end % End loop for t_loop
    
    % Parameter estimation with Least Squares
    try
        theta = X \ y;
        
        % Extract estimated parameters
        phi1_est = theta(1);
        phi2_est = theta(2);
        beta_est = theta(3 : 3 + hist_len - 1);
        gamma_est = theta(3 + hist_len : end);
        
        % Age Group Ranking
        if age < 13
            ageGroup = 'Child';
        elseif age >= 13 && age < 20
            ageGroup = 'Adolescent';
        else
            ageGroup = 'Adult';
        end
       
        % Save results
        estimatedParams(i).PatientID = patientID; % Save the original ID (numeric or string)
        estimatedParams(i).Age = age;
        estimatedParams(i).AgeGroup = ageGroup;
        estimatedParams(i).numObservations = n_obs;
        estimatedParams(i).phi1 = phi1_est;
        estimatedParams(i).phi2 = phi2_est;
        estimatedParams(i).beta = beta_est(:)';   % Ensures both row vector
        estimatedParams(i).gamma = gamma_est(:)'; % Ensures both row vector
        estimationSuccess(i) = true; % check as successful estimate
     end
end % End patient loop

%% Saving datasets with parameters
numSubjects = length(estimatedParams);
PatientID = cell(numSubjects,1);
Age = zeros(numSubjects,1);
AgeGroup = cell(numSubjects,1);
numObs = zeros(numSubjects,1);
phi1 = zeros(numSubjects,1);
phi2 = zeros(numSubjects,1);
beta_cols = zeros(numSubjects, hist_len);
gamma_cols = zeros(numSubjects, hist_len);

for i = 1:numSubjects
    PatientID{i} = string(estimatedParams(i).PatientID);
    Age(i) = estimatedParams(i).Age;
    AgeGroup{i} = estimatedParams(i).AgeGroup;
    numObs(i) = estimatedParams(i).numObservations;
    phi1(i) = estimatedParams(i).phi1;
    phi2(i) = estimatedParams(i).phi2;
    beta_cols(i,:) = estimatedParams(i).beta;
    gamma_cols(i,:) = estimatedParams(i).gamma;
end

% Creating the Final Table
T = table(PatientID, Age, AgeGroup, numObs, phi1, phi2, ...
    beta_cols(:,1), beta_cols(:,2), beta_cols(:,3), beta_cols(:,4), beta_cols(:,5), beta_cols(:,6), ...
    gamma_cols(:,1), gamma_cols(:,2), gamma_cols(:,3), gamma_cols(:,4), gamma_cols(:,5), gamma_cols(:,6), ...
    'VariableNames', {'PatientID','Age','AgeGroup','NumObservations', 'phi1', 'phi2', ...
                      'beta1','beta2','beta3','beta4','beta5','beta6', ...
                      'gamma1','gamma2','gamma3','gamma4','gamma5','gamma6'});

writetable(T, 'Subject_Parameters_and_Variables.csv');
disp('Saving completed in "Subject_Parameters_and_Variables.csv".');

%% Visualization of Results by Age Group (Beta and Gamma Parameters)
figure('Name', 'Estimation of Beta and Gamma Parameters by Age Group', 'NumberTitle', 'off', 'Position', [100, 100, 1200, 500]); % Finestra più larga

% Colors and markers to distinguish patients
numEstimated = length(estimatedParams);
colors = lines(max(numEstimated,1)); % Use 'lines' for distinct colors, minimum 1
markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '*'}; 
marker_indices = 1:length(markers);
% --- Subplot 1: Child ---
subplot(1, 3, 1);
hold on;
title('Child (< 13 years)');
xlabel('Coefficient Index (1-6)');
ylabel('Estimated Value');
legend_entries_child = {};
patient_idx_child = 0;
for i = 1:length(estimatedParams)
    if strcmp(estimatedParams(i).AgeGroup, 'Child')
        patient_idx_child = patient_idx_child + 1;
        plot(1:hist_len, estimatedParams(i).beta, '-', 'Color', colors(patient_idx_child,:), 'Marker', markers{mod(patient_idx_child-1, length(markers))+1}, 'DisplayName', sprintf('Beta - Paz. %s', string(estimatedParams(i).PatientID)));
        plot(1:hist_len, estimatedParams(i).gamma, '-.', 'Color', colors(patient_idx_child,:), 'Marker', markers{mod(patient_idx_child-1, length(markers))+1}, 'DisplayName', sprintf('Gamma - Paz. %s', string(estimatedParams(i).PatientID)));
        legend_entries_child{end+1} = sprintf('Beta - Paz. %s', string(estimatedParams(i).PatientID));
        legend_entries_child{end+1} = sprintf('Gamma - Paz. %s', string(estimatedParams(i).PatientID));
    end
end

legend('show', 'Location', 'best');
grid on;
hold off;
% --- Subplot 2: Adolescent ---
subplot(1, 3, 2);
hold on;
title('Adolescent (13-19 years)');
xlabel('Coefficient Index (1-6)');
ylabel('Estimated Value');
legend_entries_adol = {};
patient_idx_adol = 0;
for i = 1:length(estimatedParams)
    if strcmp(estimatedParams(i).AgeGroup, 'Adolescent')
        patient_idx_adol = patient_idx_adol + 1;
        plot(1:hist_len, estimatedParams(i).beta, '-', 'Color', colors(patient_idx_adol,:), 'Marker', markers{mod(patient_idx_adol-1, length(markers))+1}, 'DisplayName', sprintf('Beta - Paz. %s', string(estimatedParams(i).PatientID)));
        plot(1:hist_len, estimatedParams(i).gamma, '-.', 'Color', colors(patient_idx_adol,:), 'Marker', markers{mod(patient_idx_adol-1, length(markers))+1}, 'DisplayName', sprintf('Gamma - Paz. %s', string(estimatedParams(i).PatientID)));
        legend_entries_adol{end+1} = sprintf('Beta - Paz. %s', string(estimatedParams(i).PatientID));
        legend_entries_adol{end+1} = sprintf('Gamma - Paz. %s', string(estimatedParams(i).PatientID));
    end
end

legend('show', 'Location', 'best');
grid on;
hold off;
% --- Subplot 3: Adult ---
subplot(1, 3, 3);
hold on;
title('Adult (>= 20 years)');
xlabel('Coefficient Index(1-6)');
ylabel('Estimated Value');
legend_entries_adult = {};
patient_idx_adult = 0;
for i = 1:length(estimatedParams)
    if strcmp(estimatedParams(i).AgeGroup, 'Adult')
        patient_idx_adult = patient_idx_adult + 1;
        plot(1:hist_len, estimatedParams(i).beta, '-', 'Color', colors(patient_idx_adult,:), 'Marker', markers{mod(patient_idx_adult-1, length(markers))+1}, 'DisplayName', sprintf('Beta - Paz. %s', string(estimatedParams(i).PatientID)));
        plot(1:hist_len, estimatedParams(i).gamma, '-.', 'Color', colors(patient_idx_adult,:), 'Marker', markers{mod(patient_idx_adult-1, length(markers))+1}, 'DisplayName', sprintf('Gamma - Paz. %s', string(estimatedParams(i).PatientID)));
        legend_entries_adult{end+1} = sprintf('Beta - Paz. %s', string(estimatedParams(i).PatientID));
        legend_entries_adult{end+1} = sprintf('Gamma - Paz. %s', string(estimatedParams(i).PatientID));
    end
end

legend('show', 'Location', 'best');
grid on;
hold off;

sgtitle('Estimated Beta and Gamma Parameters by Age Group');


