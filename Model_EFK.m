%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of parameters with the method of least squares and use
% of estimated parameters for Extended Kalman Filter (EKF) simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Cleaning
clear all; close all; clc;

%% Data Upload
dataTable = readtable('Dataset_glucosio_finale.csv');

%% Model parameters and estimation
v = 5;         % Delay for insulin (minutes)
w = 5;         % Delay for carbohydrates (minutes)
hist_len = 6;  % Length of history for beta and gamma

% EKF Covariances
Q = 0.01;  % Process noise
R = 1;     % Measurement noise

%% Patient management
if iscell(dataTable.PatientID)
    uniquePatients = unique(string(dataTable.PatientID));
else
    uniquePatients = unique(dataTable.PatientID);
end
numPatients = length(uniquePatients);

% Structure for EKF results
simulationResults = struct('PatientID', [], 'xi', [], 'z', [], 'x_est_store', [], 'AgeGroup', [], 'MSE', [], 'R2', []);


for i = 1:numPatients
    patientID = uniquePatients(i);
    
    % Filtering data by patient
    if iscell(dataTable.PatientID)
        patientData = dataTable(strcmp(string(dataTable.PatientID), patientID), :);
    else
        patientData = dataTable(dataTable.PatientID == patientID, :);
    end
    
    patientData = sortrows(patientData, 'Minute');
    
    % Data extraction
    glucose_measured = patientData.Glucose; 
    insulin = patientData.Insulin;
    carbs = patientData.CarbsIntake;
    n_obs = height(patientData);
    
    % Determination of the age group
    age = patientData.Age(1); 
    if ~isnumeric(age)
        age = str2double(string(age));
    end

    if age < 13
        ageGroup = 'Child';
    elseif age >= 13 && age < 20
        ageGroup = 'Adolescent';
    else
        ageGroup = 'Adult';
    end
    
    % Calculation of parameters using the least squares method
    start_idx = max([2, v + hist_len, w + hist_len]);
    num_estimation_points = n_obs - start_idx + 1;
    X = zeros(num_estimation_points, 2 + 2*hist_len); 
    y_vec = zeros(num_estimation_points,1);
    
    for t = start_idx:n_obs
        row = t - start_idx + 1;
        X(row,1) = glucose_measured(t-1);
        X(row,2) = glucose_measured(t-2);
        for j = 1:hist_len
            idx_ins = t - v - (j-1);
            if idx_ins > 0 && idx_ins <= n_obs
                X(row,2+j) = insulin(idx_ins);
            else
                X(row,2+j) = 0;
            end
        end
        for j = 1:hist_len
            idx_carb = t - w - (j-1);
            if idx_carb > 0 && idx_carb <= n_obs
                X(row,2+hist_len+j) = carbs(idx_carb);
            else
                X(row,2+hist_len+j) = 0;
            end
        end
        y_vec(row) = glucose_measured(t);
    end
    theta = X \ y_vec;
    phi1 = theta(1);
    phi2 = theta(2);
    beta = theta(3:2+hist_len);
    gamma = theta(2+hist_len+1:end);
    
    % EKF Simulation
    P = 1;
    x_est = glucose_measured;
    for t = start_idx:n_obs
        insulin_sum = sum(beta .* insulin(max(1, t - v - (0:hist_len-1))));
        carbs_sum = sum(gamma .* carbs(max(1, t - w - (0:hist_len-1))));
        x_pred = phi1 * x_est(t-1) + phi2 * x_est(t-2) + insulin_sum + carbs_sum;
        P_pred = P + Q;
        K = P_pred / (P_pred + R);
        x_est(t) = x_pred + K * (glucose_measured(t) - x_pred);
        P = (1 - K) * P_pred;
    end
    
    % MSE and R² calculation
    y_true = glucose_measured(start_idx:end);
    y_pred = x_est(start_idx:end);
    MSE = mean((y_true - y_pred).^2);
    R2 = 1 - sum((y_true - y_pred).^2) / sum((y_true - mean(y_true)).^2);
    
    % Saving Results
    simulationResults(i).PatientID = patientID;
    simulationResults(i).xi = glucose_measured;
    simulationResults(i).z = glucose_measured;
    simulationResults(i).x_est_store = x_est;
    simulationResults(i).AgeGroup = ageGroup;
    simulationResults(i).MSE = MSE;
    simulationResults(i).R2 = R2;
    
    % Print results by patient
    fprintf('Patient %s: MSE = %.3f, R² = %.3f\n', string(patientID), MSE, R2);
end

%% Visualization Results by Age and Average MSE/R²
ageGroups = {'Child', 'Adolescent', 'Adult'};
titles = {'Child (< 13 years)', 'Adolescent (13-19 years)', 'Adult (>= 20 years)'};

fprintf("\n--- Average results by age group ---\n");

for g = 1:length(ageGroups)
    groupName = ageGroups{g};
    patientsInGroup = find(strcmp({simulationResults.AgeGroup}, groupName));
    
    mse_values = [simulationResults(patientsInGroup).MSE];
    r2_values = [simulationResults(patientsInGroup).R2];
    
    fprintf('Group %s - Average MSE: %.3f, Average R²: %.3f\n', groupName, mean(mse_values), mean(r2_values));
    
    figure('Name', titles{g}, 'NumberTitle', 'off', 'Position', [100, 100, 1200, 800]);
    sgtitle(titles{g});
    
    for p = 1:length(patientsInGroup)
        patientIdx = patientsInGroup(p);
        subplot(ceil(sqrt(length(patientsInGroup))), ceil(sqrt(length(patientsInGroup))), p);
        hold on;
        plot(simulationResults(patientIdx).xi, '.', 'DisplayName', 'Glucosio Vero');
        plot(simulationResults(patientIdx).x_est_store, 'r', 'DisplayName', 'Stima EKF');
        title(['Patient ', string(simulationResults(patientIdx).PatientID)]);
        xlabel('Time (minutes)');
        ylabel('Glucose (mg/dL)');
        grid on;
        hold off;
    end
     legend('Location', 'best');
end


