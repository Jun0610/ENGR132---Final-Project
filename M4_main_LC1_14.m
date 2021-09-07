function M4_main_LC1_14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENGR 132 
% Program Description 
% This function calculates the average v0, Vm and Km of the PGOX50
% enzyme, and the 5 enzymes from the large enzyme test dataset. 
% It then creates a Michaelis-Menten model using the calculated Km and Vmax
% values and plots the model and v0 values on a Michaelis-Menten graph. It
% also calculates the SSE between the model and the v0 values. 
%
% Function Call
% M4_main_LC1_14
%
% Input Arguments
% none
%
% Output Arguments
% none
%
% Assignment Information
%   Assignment:     M4
%   Team member:    Kevin Crowley, crowlek@purdue.edu
%                   Shravan Ranganathan, rangana2@purdue.edu
%                   Tyler Maslak, tmaslak@purdue.edu
%                   Jun Shern Lim, lim321@purdue.edu 
%   Team ID:        LC1-14
%
%   Academic Integrity:
%     [] We worked with one or more peers but our collaboration
%        maintained academic integrity.
%     Peers we worked with: Name, login@purdue [repeat for each]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INITIALIZATION
PGOX50_data = readmatrix("Data_PGOX50_enzyme.csv"); %importing PGOX50 enzyme data
enzy_data_input = PGOX50_data(4:1226, 2:11);    %parsing out only the product concentration data (uM)
sub_conc = PGOX50_data(2, 2:11);    %parsing out substrate concentration data (uM)

test_enzyme_data = readmatrix("Data_nextGen_KEtesting_allresults.xlsx");
enzymeA_org = test_enzyme_data(3:7486, 2:11);   %parsing out enzyme A original product concentration data (uM)
enzymeA_dup = test_enzyme_data(3:7486, 12:21);  %parsing out enzyme A duplicate product concentration data (uM)
enzymeB_org = test_enzyme_data(3:7486, 22:31);  %parsing out enzyme B original product concentration data (uM)
enzymeB_dup = test_enzyme_data(3:7486, 32:41);  %parsing out enzyme B duplicate product concentration data (uM)
enzymeC_org = test_enzyme_data(3:7486, 42:51);  %parsing out enzyme C original product concentration data (uM)
enzymeC_dup = test_enzyme_data(3:7486, 52:61);  %parsing out enzyme C duplicate product concentration data (uM)
enzymeD_org = test_enzyme_data(3:7486, 62:71);  %parsing out enzyme D original product concentration data (uM)
enzymeD_dup = test_enzyme_data(3:7486, 72:81);  %parsing out enzyme D duplicate product concentration data (uM)
enzymeE_org = test_enzyme_data(3:7486, 82:91);  %parsing out enzyme E original product concentration data (uM)
enzymeE_dup = test_enzyme_data(3:7486, 92:101); %parsing out enzyme E duplicate product concentration data (uM)

%% ____________________
%% CALCULATIONS
%calling the algorithm to calculate v0, Vmax and Km
[v0_ref, Vmax_ref, Km_ref] = M4_Algorithm_LC1_14(enzy_data_input, 0, sub_conc);

[v0_enzA, VmaxA, KmA] = M4_Algorithm_LC1_14(enzymeA_org,enzymeA_dup, sub_conc);  %calculating parameters for enzyme A
[v0_enzB, VmaxB, KmB] = M4_Algorithm_LC1_14(enzymeB_org,enzymeB_dup, sub_conc);  %calculating parameters for enzyme B
[v0_enzC, VmaxC, KmC] = M4_Algorithm_LC1_14(enzymeC_org,enzymeC_dup, sub_conc);  %calculating parameters for enzyme C
[v0_enzD, VmaxD, KmD] = M4_Algorithm_LC1_14(enzymeD_org,enzymeD_dup, sub_conc);  %calculating parameters for enzyme D
[v0_enzE, VmaxE, KmE] = M4_Algorithm_LC1_14(enzymeE_org,enzymeE_dup, sub_conc);  %calculating parameters for enzyme E

%calculating SSE and MAPE between calculated v0 values and model
MM_model_SSE = (Vmax_ref .* sub_conc) ./ (Km_ref + sub_conc);   %Reference Michaellis-Menten Model based on calculated Vmax and Km
SSE_ref = sum((v0_ref - MM_model_SSE).^2);   %Reference SSE between model and calculated v0 values
mean_abs_error_ref = 100/10 * sum(abs((v0_ref - MM_model_SSE) / v0_enzA)); 

MM_model_SSE_A = (VmaxA .* sub_conc) ./ (KmA + sub_conc);   %Enzyme A Michaellis-Menten Model based on calculated Vmax and Km
SSE_A = sum((v0_enzA - MM_model_SSE_A).^2);   %Enzyme A SSE between model and calculated v0 values
mean_abs_errorA = 100/10 * sum(abs((v0_enzA - MM_model_SSE_A) / v0_enzA)); 

MM_model_SSE_B = (VmaxB .* sub_conc) ./ (KmB + sub_conc);   %Enzyme B Michaellis-Menten Model based on calculated Vmax and Km
SSE_B = sum((v0_enzB - MM_model_SSE_B).^2);  %Enzyme B SSE between model and calculated v0 values
mean_abs_errorB = 100/10 * sum(abs((v0_enzB - MM_model_SSE_B) / v0_enzB)); 

MM_model_SSE_C = (VmaxC .* sub_conc) ./ (KmC + sub_conc);   %Enzyme C Michaellis-Menten Model based on calculated Vmax and Km
SSE_C = sum((v0_enzC - MM_model_SSE_C).^2);   %Enzyme C SSE between model and calculated v0 values
mean_abs_errorC = 100/10 * sum(abs((v0_enzC - MM_model_SSE_C) / v0_enzC));

MM_model_SSE_D = (VmaxD .* sub_conc) ./ (KmD + sub_conc);   %Enzyme D Michaellis-Menten Model based on calculated Vmax and Km
SSE_D = sum((v0_enzD - MM_model_SSE_D).^2);   %Enzyme D SSE between model and calculated v0 values
mean_abs_errorD = 100/10 * sum(abs((v0_enzD - MM_model_SSE_D) / v0_enzD)); 

MM_model_SSE_E = (VmaxE .* sub_conc) ./ (KmE + sub_conc);   %Enzyme E Michaellis-Menten Model based on calculated Vmax and Km
SSE_E = sum((v0_enzE - MM_model_SSE_E).^2);   %Enzyme E SSE between model and calculated v0 values
mean_abs_errorE = 100/10 * sum(abs((v0_enzE - MM_model_SSE_E) / v0_enzE));

%% ____________________
%% FORMATTED TEXT/FIGURE DISPLAYS
%plots the Michaelis-Menten model with the calculated v0 values for PGO-X50
%enzyem
figure(1)
modelx = 0:2000;      %x-axis data for model (uM) 
MM_model = (Vmax_ref .* modelx) ./ (Km_ref + modelx);   %michaellis-menten model based on Vmax and Km 
plot(sub_conc, v0_ref, "b*")
grid on
hold on
plot(modelx, MM_model, "r-")
legend("v0 Data", "Michaelis-Menten Model", 'Location', 'southeast')
xlabel("Substrate Concentration (uM)")
ylabel("Reaction Velocity (uM/s)")
title("Michaelis-Menten Plot for PGOX50 enzyme (Algorithm)")

%Michaelis-Menten plot with the calculated v0 values for enzyme A
figure(2)
MM_modelA = (VmaxA .* modelx) ./ (KmA + modelx);   %michaellis-menten model of enzyme A based on Vmax and Km 
plot(sub_conc, v0_enzA, "b*")
grid on
hold on
plot(modelx, MM_modelA, "r-")
legend("v0 Data", "Michaelis-Menten Model", 'Location', 'southeast')
xlabel("Substrate Concentration (uM)")
ylabel("Reaction Velocity (uM/s)")
title("Michaelis-Menten Plot for NextGen-A Enzyme")

%Michaelis-Menten plot with the calculated v0 values for enzyme B
figure(3)
MM_modelB = (VmaxB .* modelx) ./ (KmB + modelx);   %michaellis-menten model of enzyme B based on Vmax and Km 
plot(sub_conc, v0_enzB, "b*")
grid on
hold on
plot(modelx, MM_modelB, "r-")
legend("v0 Data", "Michaelis-Menten Model", 'Location', 'southeast')
xlabel("Substrate Concentration (uM)")
ylabel("Reaction Velocity (uM/s)")
title("Michaelis-Menten Plot for NextGen-B Enzyme")

%Michaelis-Menten plot with the calculated v0 values for enzyme C
figure(4)
MM_modelC = (VmaxC .* modelx) ./ (KmC + modelx);   %michaellis-menten model of enzyme C based on Vmax and Km 
plot(sub_conc, v0_enzC, "b*")
grid on
hold on
plot(modelx, MM_modelC, "r-")
legend("v0 Data", "Michaelis-Menten Model", 'Location', 'southeast')
xlabel("Substrate Concentration (uM)")
ylabel("Reaction Velocity (uM/s)")
title("Michaelis-Menten Plot for NextGen-C Enzyme")

%Michaelis-Menten plot with the calculated v0 values for enzyme D
figure(5)
MM_modelD = (VmaxD .* modelx) ./ (KmD + modelx);   %michaellis-menten model of enzyme D based on Vmax and Km 
plot(sub_conc, v0_enzD, "b*")
grid on
hold on
plot(modelx, MM_modelD, "r-")
legend("v0 Data", "Michaelis-Menten Model", 'Location', 'southeast')
xlabel("Substrate Concentration (uM)")
ylabel("Reaction Velocity (uM/s)")
title("Michaelis-Menten Plot for NextGen-D Enzyme")

%Michaelis-Menten plot with the calculated v0 values for enzyme E
figure(6)
MM_modelE = (VmaxE .* modelx) ./ (KmE + modelx);   %michaellis-menten model of enzyme E based on Vmax and Km 
plot(sub_conc, v0_enzE, "b*")
grid on
hold on
plot(modelx, MM_modelE, "r-")
legend("v0 Data", "Michaelis-Menten Model", 'Location', 'southeast')
xlabel("Substrate Concentration (uM)")
ylabel("Reaction Velocity (uM/s)")
title("Michaelis-Menten Plot for NextGen-E Enzyme")
%% ____________________
%% RESULTS
%displays v0, Vmax and Km for each enzyme
%PGO-X50 Enzyme
fprintf("PGO-X50 Enzyme\n")
fprintf("v0 PGO-X50 enzyme (uM/s) = ")
fprintf("\n")
disp(v0_ref)
fprintf("Vmax PGO-X50 Enzyme  = %.4fuM/s\n", Vmax_ref)
fprintf("Km PGO-X50 Enzyme  = %.4fuM\n", Km_ref)
fprintf("SSE PGO-X50 Enzyme  = %.4f(uM/s)^2\n", SSE_ref);
fprintf("Mean Absolute Percent Error of PGO-X50 Enzyme = %.4f%%\n", mean_abs_error_ref)

%enzyme A
fprintf("\nNextGen-A Enzyme\n")
fprintf("v0 enzymeA (uM/s) = ")
fprintf("\n")
disp(v0_enzA)
fprintf("Vmax Enzyme A = %.4fuM/s\n", VmaxA)
fprintf("Km Enzyme A = %.4fuM\n", KmA)
fprintf("SSE Enzyme A = %.4f(uM/s)^2\n", SSE_A);
fprintf("Mean Absolute Percent Error Enzyme A = %.4f%%\n", mean_abs_errorA)

%enzyme B
fprintf("\nNextGen-B Enzyme\n")
fprintf("v0 enzymeB (uM/s) = ")
fprintf("\n")
disp(v0_enzB)
fprintf("Vmax Enzyme B = %.4fuM/s\n", VmaxB)
fprintf("Km Enzyme B = %.4fuM\n", KmB)
fprintf("SSE Enzyme B = %.4f(uM/s)^2\n", SSE_B);
fprintf("Mean Absolute Percent Error Enzyme B = %.4f%%\n", mean_abs_errorB)

%enzyme C
fprintf("\nNextGen-C Enzyme\n")
fprintf("v0 enzymeC (uM/s) = ")
fprintf("\n")
disp(v0_enzC)
fprintf("Vmax Enzyme C = %.4fuM/s\n", VmaxC)
fprintf("Km Enzyme C = %.4fuM\n", KmC)
fprintf("SSE Enzyme C = %.4f(uM/s)^2\n", SSE_C);
fprintf("Mean Absolute Percent Error Enzyme C = %.4f%%\n", mean_abs_errorC)

%enzyme D
fprintf("\nNextGen-D Enzyme\n")
fprintf("v0 enzymeD (uM/s) = ")
fprintf("\n")
disp(v0_enzD)
fprintf("Vmax Enzyme D = %.4fuM/s\n", VmaxD)
fprintf("Km Enzyme D = %.4fuM\n", KmD)
fprintf("SSE Enzyme D = %.4f(uM/s)^2\n", SSE_D);
fprintf("Mean Absolute Percent Error Enzyme D = %.4f%%\n", mean_abs_errorD)

%enzyme E
fprintf("\nNextGen-E Enzyme\n")
fprintf("v0 enzymeE (uM/s) = ")
fprintf("\n")
disp(v0_enzE)
fprintf("Vmax Enzyme E = %.4fuM/s\n", VmaxE)
fprintf("Km Enzyme E = %.4fuM\n", KmE)
fprintf("SSE Enzyme E = %.4f(uM/s)^2\n", SSE_E);
fprintf("Mean Absolute Percent Error Enzyme E = %.4f%%\n", mean_abs_errorE)

%% ____________________
%% ACADEMIC INTEGRITY STATEMENT
% We have not used source code obtained from any other unauthorized
% source, either modified or unmodified. Neither have we provided
% access to my code to another. The program we are submitting
% is our own original work.
