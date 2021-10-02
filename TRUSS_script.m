% Coursework - Collapse of a Redundant Truss
% Ken Tan % 
% Script file

clear;
close all;
clc;

%% Formulation of truss geometry by specifying nodal coordinates, dofs and elemental connections
% Parameters
% E, A and L are constant throughout this truss geometry, but they can always
% be expressed in arrays based on each bar's properties
E = 2e5; % N/mm2 - modulus of elasticity of each bar element
R = 20; % mm2 - cross-sectional area of each bar element
L = 1000; %  mm - horizontal/vertical bar length
yieldS = 250; % MPA - yield stress value

P = 100000; % N - point end load

% Specifying nodal x-y coordinates
NODES.coords = [0   0 ; % node 1 - x,y coordinates
                0   L ; % node 2 - x,y coordinates
                L   0 ; % node 3 - x,y coordinates
                L   L ; % node 4 - x,y coordinates
                2*L 0 ; % node 5 - x,y coordinates
                2*L L ; % node 6 - x,y coordinates
                3*L 0 ; % node 7 - x,y coordinates
                3*L L ; % node 8 - x,y coordinates
                4*L 0 ; % node 9 - x,y coordinates
                4*L L]; % node 10 - x,y coordinates
            
NODES.dofs = [1 2 ;  % node 1 - dofs u1,v1 (1,2) - FIXED  
              3 4 ;  % node 2 - dofs u2,v2 (3,4) - FIXED
              5 6 ;  % node 3 - dofs u3,v3 (5,6) - FREE
              7 8 ;  % node 4 - dofs u4,v4 (7,8) - FREE
              9 10 ; % node 5 - dofs u5,v5 (9,10) - FREE
              11 12 ; % node 6 - dofs u6,v6 (11,12) - FREE
              13 14 ; % node 7 - dofs u7,v7 (13,14) - FREE
              15 16 ; % node 8 - dofs u8,v8 (15,16) - FREE
              17 18 ; % node 9 - dofs u9,v9 (17,18) - FIXED
              19 20]; % node 10 - dofs u10,v10 (19,20) - FIXED
          
% Note that node numbers are identified by their row number
     
% Specifying element nodal connectivity (order does not matter)
ELEMENTS = [1 3 ; % element 1 - element dofs u1,v1,u3,v3 (1,2,5,6)
            1 4 ; % element 2 - element dofs u1,v1,u4,v4 (1,2,7,8)
            2 3 ; % element 3 - element dofs u2,v2,u3,v3 (3,4,5,6)
            2 4 ; % element 4 - element dofs u2,v2,u4,v4 (3,4,7,8)
            4 3 ; % element 5 - element dofs u4,v4,u3,v3 (7,8,5,6)
            3 5 ; % element 6 - element dofs u3,v3,u5,v5 (5,6,9,10)
            3 6 ; % element 7 - element dofs u3,v3,u5,v5 (5,6,11,12)
            4 5 ; % element 8 - element dofs u4,v4,u5,v5 (7,8,9,10)
            4 6 ; % element 9 - element dofs u4,v4,u6,v6 (7,8,11,12)
            6 5 ; % element 10 - element dofs u6,v6,u5,v5 (11,12,9,10)
            5 7 ; % element 11 - element dofs u5,v5,u7,v7 (9,10,13,14)
            5 8 ; % element 12 - element dofs u5,v5,u8,v8 (9,10,15,16)
            6 7 ; % element 13 - element dofs u6,v6,u7,v7 (11,12,13,14)
            6 8 ; % element 14 - element dofs u6,v6,u7,v7 (11,12,15,16)
            8 7 ; % element 15 - element dofs u8,v8,u7,v7 (17,18,15,16)
            7 9 ; % element 16 - element dofs u7,v7,u9,v9 (15,16,19,20)
            7 10 ; % element 17 - element dofs u7,v7,u10,v10 (15,16,21,22)
            8 9 ; % element 18 - element dofs u8,v8,u9,v9 (17,18,19,20)
            8 10]; % element 19 - element dofs u8,v8,u10,v10 (17,18,21,22)

% Creating a global force vector based on external forces applied to the
% truss system
% Initialising an empty vector of size 2*nodes
F = zeros(2*length(NODES.coords), 1);

% Constructing the global nodal force vector (initial loads exerted on the
% system)
F(8) = -3*P;
F(12) = -1.5*P;
F(16) = -P;

% Stored degrees of freedom
dofsFree = 5:16; % unknown nodal x-y dofs (at nodes 3,4,5,6,7,8)
dofsRestrained = [1:4,17:20]; % known nodal x-y dofs due to BC (at nodes 1,2,9,10)

%% Question 2
% MNA by iLA
% Create the TRUSS object called myTruss
% For the second last input, type = 1 refers to MNA by iLA, type = 2 refers to
% GMNA by iLa. The last input is just to indicate whether or not force
% equilibrium checks should be printed out (please see Line 77 of Truss.m)
myTruss = TRUSS(E, R, yieldS, NODES.coords, NODES.dofs, ELEMENTS, F, dofsFree, dofsRestrained, 1, 1);

% Intialise the global increment counter
I = 1;
for i = 1:size(ELEMENTS,1)
    fprintf(['Initialising collapse analysis ... Event ', num2str(I)]);
    disp(' ');
    myTruss = myTruss.assembly();
    myTruss = myTruss.invertible();
    % Break out of the while loop if the KFF matrix is not invertible
    if myTruss.rcondValue < 1e-15
        fprintf('An additional member cannot be removed, otherwise the truss will collapse \n');
        fprintf('Event 7 will not happen \n');
        break
    end
    
    myTruss = myTruss.solver();
    disp(' ');
    myTruss = myTruss.postproc();
    myTruss = myTruss.plot();
    if myTruss.type == 1
        title(['MNA Event #', num2str(I), ', \lambda_{tot} = ', num2str(myTruss.lambdatot)]);
    elseif myTruss.type == 2
        title(['GMNA Event #', num2str(I), ', \lambda_{tot} = ', num2str(myTruss.lambdatot)]);
    end
    set(gca, 'FontSize', 30);
    myTruss = myTruss.update();
    
    I = I + 1;
end

fprintf('The total load proportionality factor is %g\n', myTruss.lambdatot);
fprintf('The number of bar members permitted to fail before global instability is %g\n', I-1);

%% Question 3
% GMNA by iLA
% Create the TRUSS object called myTruss
% For the second last input, type = 1 refers to MNA by iLA, type = 2 refers to
% GMNA by iLa. The last input is just to indicate whether or not force
% equilibrium checks should be printed out (please see Line 77 of Truss.m)
myTruss = TRUSS(E, R, yieldS, NODES.coords, NODES.dofs, ELEMENTS, F, dofsFree, dofsRestrained, 2, 1);

% Intialise the global increment counter
I = 0;
for i = 1:size(ELEMENTS,1)
    fprintf(['Initialising collapse analysis ... Event ', num2str(I)]);
    disp(' ');
    myTruss = myTruss.assembly();
    myTruss = myTruss.solver();
    disp(' ');
    % Break out of the while loop if the KFF matrix is not invertible
    myTruss = myTruss.invertible();
    if myTruss.rcondValue < 1e-15
        fprintf('An additional member cannot be removed, otherwise the truss will collapse \n');
        break
    end

    myTruss = myTruss.postproc();
    myTruss = myTruss.plot();
    if myTruss.type == 1
        title(['MNA Event #', num2str(I), ', \lambda_{tot} = ', num2str(myTruss.lambdatot)]);
    elseif myTruss.type == 2
        title(['GMNA Event #', num2str(I), ', \lambda_{tot} = ', num2str(myTruss.lambdatot)]);
    end
    
    myTruss = myTruss.update();
    I = I + 1;
end

fprintf('The total load proportionality factor is %g\n', myTruss.lambdatot);
fprintf('The number of bar members permitted to fail before global instability is %g\n', I);

%% Question 5
R_varied = 1:40;

% Initialise an empty array to log the number of members at failure for
% each analysis
failed_MNA = zeros(length(R_varied), 1);
failed_GMNA = zeros(length(R_varied), 1);

% Initialise an empty array to log the lambda_tot values for each analysis
lambdatot_MNA = zeros(length(R_varied), 1);
lambdatot_GMNA = zeros(length(R_varied), 1);

% Disable singular matrix/ nearly singular matrix warning (as we are
% looping through the loop many times we don't want to fill the command
% window with the same warning knowing that we will compute until near
% singular anyways
warning('off', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:singularMatrix');

% MNA by iLA
for i = 1:length(R_varied)
    % Create the truss object based on the radius
    myTrussMNA = TRUSS(E, R_varied(i), yieldS, NODES.coords, NODES.dofs, ELEMENTS, F, dofsFree, dofsRestrained, 1);
    myTrussGMNA = TRUSS(E, R_varied(i), yieldS, NODES.coords, NODES.dofs, ELEMENTS, F, dofsFree, dofsRestrained, 2);
    
    % Initialise the increment counters as 0
    I_MNA = 0;
    I_GMNA = 0;
    
    % MNA analysis
    for j = 1:size(ELEMENTS, 1)
       
        myTrussMNA = myTrussMNA.assembly();
        myTrussMNA = myTrussMNA.solver();

        % Break out of the while loop if the KFF matrix is not invertible
        myTrussMNA = myTrussMNA.invertible();
        if myTrussMNA.rcondValue < 1e-15
            failed_MNA(i) = I_MNA;
            lambdatot_MNA(i) = myTrussMNA.lambdatot;
            break
        end
        
        myTrussMNA = myTrussMNA.postproc();
        % We have omitted the 'plotting' method of the class as it will not
        % be required in this part. It was indeed the right choice to
        % seperate the 'plotting' method from the post processor!
        
        myTrussMNA = myTrussMNA.update();
        I_MNA = I_MNA + 1;
    end
    
    % GMNA analysis
    for k = 1:size(ELEMENTS, 1)
        
        myTrussGMNA = myTrussGMNA.assembly();
        myTrussGMNA = myTrussGMNA.solver();
        
        % Break out of the while loop if the KFF matrix is not invertible
        myTrussGMNA = myTrussGMNA.invertible();
        if myTrussGMNA.rcondValue < 1e-15
            failed_GMNA(i) = I_GMNA;
            lambdatot_GMNA(i) = myTrussGMNA.lambdatot;
            break
        end
        
        myTrussGMNA = myTrussGMNA.postproc();
        % Again, we do not require the 'plotting' method from the class
        % file
        
        myTrussGMNA = myTrussGMNA.update();
        I_GMNA = I_GMNA + 1;
    end
end

% Plotting the total load proportionality factor lambda_tot as a function
% of radius varied from 1mm to 40 mm in steps of 1 mm.
figure;
plot(R_varied, lambdatot_MNA, 'r', R_varied, lambdatot_GMNA, 'b--', 'LineWidth', 3); grid on
xlabel('Radius [mm]'); ylabel('\lambda_{tot}'); title('Plot of \lambda_{tot} against radius of member');
legend('MNA by iLA', 'GMNA by iLA');

% Plotting the number of failed members against radius
figure;
plot(R_varied, failed_MNA, 'b-', R_varied, failed_GMNA, 'r--', 'LineWidth', 3);
xlabel('Radius [mm]'); ylabel('Number of failed members'); title('Plot of Number of Failed Members Against Radius');
legend('MNA by iLA', 'GMNA by iLA'); ylim([3 8]);

% End of script file %


    

















