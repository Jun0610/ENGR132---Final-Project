function [avg_v0, Vmax, Km] = M4_Algorithm_LC1_14(enzyData_org, enzyData_dup, sub_conc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENGR 132
% Program Description
% This function calculates the average v0, Vm and Km of a given
% an enzyme product concentration test dataset and its duplicate.
%
% Function Call
% [avg_v0, Vmax, Km] = M4_Algorithm_LC1_14(enzyData_org, enzyData_dup, sub_conc)
%
% Input Arguments
% 1. enzyme original product concentration test data: enzyData_org (uM)
% 2. enzyme duplicate product concentration test data: enzyData_dup (uM), IF
%  enzyData_dup = 0, it will be omitted from calculation, and only the
%  original data will be used to calculate v0, Vmax, and Km
% 3. substrate concentration: sub_conc (uM)
%
% Output Arguments
% 1. average v0 of orignal and duplicate tests: avg_v0 (uM/s)
% 2. maximum enzyme kinematic velocity: Vmax (uM/s)
% 3. constant, Km: Km (uM)
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

%% ____________________
%% INITIALIZATION
[numrows, numcols] = size(enzyData_org); %determines the size of the original enzyme test data, which wil also be the size of the duplicate test data
count = 1;
%slope = zeros(1, 50);   %creates a 1x50 zero matrix for slope
%% Improvement 1: takes the slope of the first 39 points instead of 50
slope = zeros(1, 39);   %creates a 1x39 zero matrix for slope 
v0_org = zeros(1, 10);  %creates a 1x10 zero matrix for v0 for the original test
v0_dup = zeros(1, 10);  %creates a 1x10 zero matrix for v0 for the duplicate test

%% ____________________
%% CALCULATIONS

%nested repetition and selection constructs to determine v0 for the
%original test cases
for cols = 1:numcols
    for rows = 1:numrows
        %parsing out product concentration data at each substrate
        %concentration
        if enzyData_org(rows, cols) >= 0
            prod_conc(count) = enzyData_org(rows, cols); %product concentration (uM)
            count = count + 1;
        end
    end
    [rows_pc, cols_pc] = size(prod_conc.');
    time = 0:rows_pc - 1;   %time (s)
    prod_conc_smooth = smoothdata(prod_conc, 'rloess'); %smoothens the product concentration data (uM)
    for i = 2:40
        %calculates each of the slopes that points 2 to 40 (39 in total)
        %make with the first point
        slope(i-1) = (prod_conc_smooth(i) -prod_conc_smooth(1)) / (time(i) - time(1));
    end
    %v0_org(cols) = mean(slope); %averages the 50 slopes to find v0 (uM/s)
    %% Improvement 1: use mode instead of mean to calculate the average of the fisrt 39 slopes
    v0_org(cols) = mode(slope); %averages the 39 slopes using mode to find v0 (uM/s)
    prod_conc = [];
    count = 1;
end

%selection structure to account for presence of duplicate data
if enzyData_dup ~= 0
    %nested repetition and selection constructs to determine v0 for the
    %duplicate test cases
    for cols = 1:numcols
        for rows = 1:numrows
            %parsing out product concentration data at each substrate
            %concentration
            if enzyData_dup(rows, cols) >= 0
                prod_conc(count) = enzyData_dup(rows, cols); %product concentration (uM)
                count = count + 1;
            end
        end
        [rows_pc, cols_pc] = size(prod_conc.');
        time = 0:rows_pc -1;    %time (s)
        prod_conc_smooth = smoothdata(prod_conc, 'rloess'); %smoothens the product concentration data (uM)
        for i = 2:40
            %calculates each of the slopes that points 2 to 40 (39 in total)
            %make with the first point
            slope(i-1) = (prod_conc_smooth(i) - prod_conc_smooth(1)) / (time(i) - time(1));
        end
        %v0_dup(cols) = mean(slope); %averages the 50 slopes to find v0 (uM/s)
        %% Improvement 1: use mode instead of mean to calculate the average of the fisrt 39 slopes
        v0_dup(cols) = mode(slope); %averages the 39 slopes using mode to find v0 (uM/s)
        prod_conc = [];
        count = 1;
    end
end

%finds the average v0 values of the original test data and duplicate data
if enzyData_dup ~= 0
    combined_v0 = [v0_org; v0_dup];
    avg_v0 = mean(combined_v0);
else
    avg_v0 = v0_org;
end

%Determining x values and y values for the Eadie-Hofstee plot
% EadieHofsteeX = avg_v0./sub_conc;
% EadieHofsteeY = avg_v0;
% 
% %Determining the coefficients of the Lineweaver-Burk plot
% coeff = polyfit(EadieHofsteeX, EadieHofsteeY, 1);
% 
% %Calculating Vmax and Km
% Vmax = coeff(2);    %Vmax (uM/s)
% Km = - coeff(1);   %Km (uM)

%% Improvement 2: Using the Hanes-Woolf plot instead of Eadie-Hofstee to
%determine Km and Vmax
HanesWoolfX = sub_conc;     % HanesWoolf X axis (uM)
HanesWoolfY = sub_conc./avg_v0; % HanesWoolf Y axis (s)

%Determining coefficients of Hanes-Woolf plot
coeff = polyfit(HanesWoolfX, HanesWoolfY, 1);

%Calculating Vmax and Km
Vmax = 1/coeff(1);
Km = coeff(2) * Vmax; 


%% ____________________
%% FORMATTED TEXT/FIGURE DISPLAYS


%% ____________________
%% RESULTS


%% ____________________
%% ACADEMIC INTEGRITY STATEMENT
% We have not used source code obtained from any other unauthorized
% source, either modified or unmodified. Neither have we provided
% access to my code to another. The program we are submitting
% is our own original work.



