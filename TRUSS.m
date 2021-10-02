% Coursework - Collapse of a Redundant Truss
% Ken Tan 
% Class TRUSS

% Class file for analysing the collapse of a redundant truss

classdef TRUSS
    
    properties
        % Geometry and Loading
        % For the case of our trusses, E and A are constant throughout. E
        % and A can be inputted into the constructor method as an array if
        % we require different values of E and A for different bar
        % elements
        E % Elastic modulus of each member
        R % Radius of member
        Le % Element length
        A % Area of each member
        yieldS % Yield stress
        PsMNA % Yield Load
        PsGMNA % Euler Buckling Load
        type % To signify either MNA by iLA or GMNA by iLA (1 represnts MNA by iLA and 2 represents GMNA by iLA)
        
        % Nodal arrangements (coordinates)
        NODES
        
        % Dof arrangements
        DOFS
        
        % Element connectivity
        ELEMENTS
        
        % Useful parameters
        nodes % Number of nodes
        elements % Number of elements
        
        % Global force vector and global displacement vector 
        F % F here refers external forces on the system, not axial forces in bar elements
        U % U here refers to the displacements of the nodal dofs
        
        % Stored Dofs of Interest
        dofsFree % Global indices of free dofs
        dofsRestrained % Global indices of restrained dofs
        
        % Stiffness matrix
        K 
        KFF % Defined globally to be checked for singularity
        
        % Stored cumulative F, U & lambda
        Ftot % Cumulative axial force within each bar element 
        Utot % Cumulative axial displacement of each bar element
        lambdatot % Cumulative load proportionality factor of the system over all previous events
        
        % Miscellaneous properties defined so that they can be called
        % globally in other methods
        indexOfElement % Vector to store the elemental indices of each failed member
        FInd % Collapse force indicator (1 by tension, -1 by compression, -2 by buckling)
        rcondValue % Reciprocal number of the Kff matrix, to be called in the main script
        print % 1 for printing force equilibrium check, anything else for don't print
    end
    
    methods
        % METHOD 1 - CLASS CONSTRUCTOR
        function obj = TRUSS(E, R, yieldS, NODES, DOFS, ELEMENTS, F, dofsFree, dofsRestrained, type, print)
            obj.E = E;
            obj.R = R;
            obj.yieldS = yieldS;
            obj.NODES = NODES;
            obj.DOFS = DOFS;
            obj.ELEMENTS = ELEMENTS;
            obj.F = F;
            obj.dofsFree = dofsFree;
            obj.dofsRestrained = dofsRestrained;
            obj.type = type;
            
            % Optional argument to print force equilibrium checks (1 to
            % print, anything else or nothing to don't print)            
            if ~exist('print', 'var') || print ~= 1
                obj.print = 0;
            else
                obj.print = print;
            end
            
            % Some useful parameters
            obj.nodes = size(obj.NODES, 1); % Represents the number of nodes
            obj.elements = size(obj.ELEMENTS, 1); % Represents the number of elements
            
            % Preallocating the vectors for storing the cumulative Ftot, Utot and
            % lambdatot vectors. This is included in the constructor class
            % as the other functions will be called more than once by being
            % in the for loop and we don't want to overwrite the
            % initialisation of the vectors
      
            obj.Ftot = zeros(obj.elements, 1);
            obj.lambdatot = 0; 
            obj.Utot = zeros(2*obj.nodes, 1); 
            obj.indexOfElement = 0;
            obj.Le = zeros(obj.elements, 1);
            

        end
        
        % METHOD 2 - ASSEMBLY OF STIFFNESS MATRIX
        function obj = assembly(obj)
            
            % Initialising an empty matrix of size (2*nodes)
            obj.K = zeros(2*obj.nodes);
            
            % Computing the cross-sectional area of each member
            obj.A = pi*(obj.R)^2;
            
            for EL = 1:obj.elements % Loop through all elements & build stiffness matrix
                criticalMemberCheck = ismember(EL, obj.indexOfElement);
                if criticalMemberCheck == 1
                    continue
                else
                    n1 = obj.ELEMENTS(EL,1); n2 = obj.ELEMENTS(EL,2); % Identify element node numbers
                    x1 = obj.NODES(n1,1); y1 = obj.NODES(n1,2); % Element node 1 - x,y coordinates
                    x2 = obj.NODES(n2,1); y2 = obj.NODES(n2,2); % Element node 2 - x,y coordinates
                    dof11 = obj.DOFS(n1,1); dof12 = obj.DOFS(n1,2); % Element node 1 - dofs
                    dof21 = obj.DOFS(n2,1); dof22 = obj.DOFS(n2,2); % Element node 2 - dofs
                    alpha = atan2(y2-y1,x2-x1); % Angle of inclination relative to the POSITIVE direction of the x axis
                    c = cos(alpha); c2 = c*c; s = sin(alpha); s2 = s*s; cs = c*s; % Angle parameters
                    obj.Le(EL) = sqrt( (x2 - x1)^2 + (y2 - y1)^2 ); % Element length
                    ke = (obj.E*obj.A)/obj.Le(EL); % Element axial stiffness.
                    % Change obj.EA to obj.EA(EL) if EA is not constant
                    
                    % Updating global stiffness matrix [K] coefficients
                    % Row 1 - element dof11
                    obj.K(dof11,dof11) = obj.K(dof11,dof11) + ke*c2; % Col 1 - element dof11
                    obj.K(dof11,dof12) = obj.K(dof11,dof12) + ke*cs; % Col 2 - element dof12
                    obj.K(dof11,dof21) = obj.K(dof11,dof21) - ke*c2; % Col 3 - element dof21
                    obj.K(dof11,dof22) = obj.K(dof11,dof22) - ke*cs; % Col 4 - element dof22
                    
                    % Row 2 - element dof12
                    obj.K(dof12,dof11) = obj.K(dof12,dof11) + ke*cs; % Col 1 - element dof11
                    obj.K(dof12,dof12) = obj.K(dof12,dof12) + ke*s2; % Col 2 - element dof12
                    obj.K(dof12,dof21) = obj.K(dof12,dof21) - ke*cs; % Col 3 - element dof21
                    obj.K(dof12,dof22) = obj.K(dof12,dof22) - ke*s2; % Col 4 - element dof22
                    
                    % Row 3 - element dof21
                    obj.K(dof21,dof11) = obj.K(dof21,dof11) - ke*c2; % Col 1 - element dof11
                    obj.K(dof21,dof12) = obj.K(dof21,dof12) - ke*cs; % Col 2 - element dof12
                    obj.K(dof21,dof21) = obj.K(dof21,dof21) + ke*c2; % Col 3 - element dof21
                    obj.K(dof21,dof22) = obj.K(dof21,dof22) + ke*cs; % Col 4 - element dof22
                    
                    % Row 4 - element dof22
                    obj.K(dof22,dof11) = obj.K(dof22,dof11) - ke*cs; % Col 1 - element dof11
                    obj.K(dof22,dof12) = obj.K(dof22,dof12) - ke*s2; % Col 2 - element dof12
                    obj.K(dof22,dof21) = obj.K(dof22,dof21) + ke*cs; % Col 3 - element dof21
                    obj.K(dof22,dof22) = obj.K(dof22,dof22) + ke*s2; % Col 4 - element dof22
                end
            end
            
            obj.KFF = obj.K(obj.dofsFree, obj.dofsFree);
        end
        
        % METHOD 3 - SOLVER
        function obj = solver(obj)
            
            % Specification of submatrices directly, note there is no need to form
            % a rearranged K or F explicitly
            KRR = obj.K(obj.dofsRestrained, obj.dofsRestrained);
            KRF = obj.K(obj.dofsRestrained, obj.dofsFree);
            KFR = obj.K(obj.dofsFree, obj.dofsRestrained);
            fF = obj.F(obj.dofsFree); % fR will be computed later after solving for uF

            % Solution for the unknown nodal dofs
            uR = zeros(length(obj.dofsRestrained)); % BC - zero displacement at restrained nodes
            uF = obj.KFF\(fF - KFR*uR); % 1st matrix equation
            obj.U = zeros(2*obj.nodes,1);
            for i = 1:length(uF)
                obj.U(obj.dofsFree(i)) = uF(i); % Storing all the displacements (restrained and free dofs) into one global U vector
            end

            % Solution for the unknown reactions
            fR = KRF*uF + KRR*uR; % 2nd matrix equation
            % Note that there is no need to store the reaction forces at
            % the restrained dofs in the global obj.F vector because they
            % will not contribute to the analysis and hence we have
            % ommitted this step

            % Check for force equilibrium. This section is skipped if
            % obj.print is anything other than '1' when initialising the
            % truss object
            if obj.print == 1
                tol = 1e-6;
                % Initialise the inital values of each sum as zero
                sumOffR_vert = 0; sumOffR_hor = 0;
                sumOffF_vert = 0; sumOffF_hor = 0;
                
                % Compute sum of vertical and horizontal forces at
                % restrained dofs
                for i = 1:length(fR)
                    if rem(i, 2) == 0
                        sumOffR_vert = sumOffR_vert + fR(i);
                    elseif rem(i, 2) == 1
                        sumOffR_hor = sumOffR_hor + fR(i);
                    end
                end
                
                % Compute sum of vertical and horizontal forces at free
                % dofs
                for i = 1:length(fF)
                    if rem(i, 2) == 0
                        sumOffF_vert = sumOffF_vert + fF(i);
                    elseif rem(i, 2) == 1
                        sumOffF_hor = sumOffF_hor + fF(i);
                    end
                end
                
                % Displaying the results
                disp(' '); disp('Force equilibrium check:');
                disp(['Total vertical reaction forces = ', num2str(sumOffR_vert)]);
                disp(['Total vertical applied forces = ', num2str(sumOffF_vert)]);
                disp(['Total horizontal reaction forces = ', num2str(sumOffR_hor)]);
                disp(['Total horizontal applied forces = ', num2str(sumOffF_hor)]);
                if abs(abs(sumOffR_vert) - abs(sumOffF_vert)) >= tol
                    disp('Vertical force equilibrium is not satisfied');
                elseif abs(abs(sumOffR_hor) - abs(sumOffF_hor)) >= tol
                    disp('Horizontal force equilibrium is not satisfied');
                else
                    disp('Force equilibrium is satisfied');
                end
            end
        end

        % METHOD 4 - POST PROCESSOR (OBTAINING FTOT, UTOT AND LAMBDATOT)
        function obj = postproc(obj)
            
            % Updating nodal coordinates after deformation
            % Creating two vectors - one with real size deformed coordinates and one
            % with amplified deformed coordinate (for plotting purposes)
            nodesNewCoords = zeros(size(obj.NODES)); 

            for i = 1:size(obj.NODES,1)
                for j = 1:size(obj.NODES,2)   
                    nodesNewCoords(i,j) = obj.NODES(i,j) + obj.U(obj.DOFS(i,j));
                end
            end

            % Preallocating some arrays
            lambdaCurrent = zeros(obj.elements, 1); % Vector to store lambda values for each element before checking for the minimum positive amongst them
            Fcurrent = zeros(obj.elements, 1); % Vector to store the isolated axial forces in each element for the current event only, to be added on later to dlambda*dF
            dup = zeros(obj.elements, 1); % Elemental change in length
            obj.FInd = zeros(obj.elements, 1); % Vector to indicate whether or not a bar will fail by tension, compression or buckling
            obj.PsGMNA = zeros(obj.elements, 1); % Vector to store the euler buckling loads of each element (based on their lengths hence we require an array for this)
            
            % Computing the squash load Ps
            obj.PsMNA = obj.A*obj.yieldS; % Since all elements have the same area and steel grade, this does not have to be an array, althought it can be 
            
            for EL = 1:obj.elements
                criticalMemberCheck = ismember(EL, obj.indexOfElement); % Check whether EL is an index of a failed element member
                
                if criticalMemberCheck == 1 
                    continue % Skip this member as it has already failed and can be assumed to be non-existent for the rest of the analysis
                else
                    n1 = obj.ELEMENTS(EL,1); n2 = obj.ELEMENTS(EL,2); % Identify element node numbers
                    
                    % Computing the Elastic Euler buckling load
                    if obj.type == 2
                        I = (pi*(obj.R^4))/4;
                        obj.PsGMNA(EL) = ((pi^2)*(obj.E*I))/obj.Le(EL)^2;
                    end
                    % Obtaining the changes in bar element length
                    x1 = obj.NODES(n1,1); y1 = obj.NODES(n1,2); % Element node 1 - x,y original coordinates
                    x2 = obj.NODES(n2,1); y2 = obj.NODES(n2,2); % Element node 2 - x,y original coordinates
                    x1_new = nodesNewCoords(n1,1); y1_new = nodesNewCoords(n1,2); % Element node 1 - x,y actual deformed coordinates
                    x2_new = nodesNewCoords(n2,1); y2_new = nodesNewCoords(n2,2); % Element node 2 - x,y actual deformed coordinates
                    ke = (obj.E*obj.A)/obj.Le(EL); % Element axial stiffness
                    alpha = atan2(y2-y1,x2-x1); % Angle of inclination relative to the POSITIVE x axis direction
                    u1 = x1_new - x1; v1 = y1_new - y1; u2 = x2_new - x2; v2 = y2_new - y2; % Reconstruction of element global dofs
                    up1 = cos(alpha)*u1 + sin(alpha)*v1; up2 = cos(alpha)*u2 + sin(alpha)*v2; dup(EL) = up2 - up1; % Reconstruction of element local dofs
                    
                    % Calculate the axial forces in each element and store
                    % it in obj.Fcurrent
                    Fcurrent(EL) = ke*dup(EL);
                    
                    % Determine the critical member (smallest non-negative dlambda) 
                    % For each bar, store 2 values of lambda in tempLambda1 and
                    % tempLambda2 (assume both tensile and compressive/buckling failure)
                    tempLambda1 = (obj.PsMNA - obj.Ftot(EL))/Fcurrent(EL); % Tensile failure check
                    tempLambda2_MNA = (-obj.PsMNA - obj.Ftot(EL))/Fcurrent(EL); % Compressive failure check
                    flag = 0;
                    if obj.type == 1 % Proceed with MNA by iLA
                        tempLambda2 = tempLambda2_MNA; 
                    elseif obj.type == 2 % Proceed with GMNA by iLA
                        tempLambda2_GMNA = (-obj.PsGMNA(EL) - obj.Ftot(EL))/Fcurrent(EL); % Compressive failure by buckling check
                        tempLambda2 = min([tempLambda2_GMNA tempLambda2_MNA]);
                        if tempLambda2 == tempLambda2_GMNA
                            flag = 1; % Flagging this element if it undergoes buckling failure prior to compressive failure
                        end
                    end
                    
                    % Choose the smallest non-negative lambda value and store it in
                    % obj.lambdaCurrent
                    if tempLambda1 > 0 %&& tempLambda1 > tempLambda2
                        lambdaCurrent(EL) = tempLambda1;
                        obj.FInd(EL) = 1; % Tensile axial force failure indicator
                    elseif tempLambda2 > 0 % && tempLambda2 > tempLambda1
                        lambdaCurrent(EL) = tempLambda2;
                        if flag == 0
                            obj.FInd(EL) = -1; % Compressive axial force failure indicator
                        elseif flag == 1
                            obj.FInd(EL) = -2; % Buckling failure indicator
                        end
                    elseif tempLambda1 < 0 && tempLambda2 < 0
                        lambdaCurrent(EL) = 999; % Storing it as 999 so that, when computing the minimum lambda value, it won't pick up the negative lambda
                    end
                end
            end

            % Compute the minimum non-negative lambda value
            sortedLambda = sort(lambdaCurrent);
            lambdaIndex = find(sortedLambda > 0, 1);
            lambdaFinal = sortedLambda(lambdaIndex);
            obj.indexOfElement(end+1) = find(lambdaCurrent == lambdaFinal); % Flagging the indices of the critical member

            % Compute Ftot, Utot and lambdatot
            for EL = 1:obj.elements
                criticalMemberCheck = ismember(EL, obj.indexOfElement);
                if criticalMemberCheck == 1 % Check if EL is a critical member
                    continue
                else
                    obj.Ftot(EL) = obj.Ftot(EL) + Fcurrent(EL)*lambdaFinal;
                end
            end

            for i = 1:2*obj.nodes
                obj.Utot(i) = obj.Utot(i) + obj.U(i)*lambdaFinal;
            end

            obj.lambdatot = obj.lambdatot + lambdaFinal;
            end
        
        % METHOD 5 - PLOTTING DEFORMATION STATE AND VISUALISATION OF AXIAL
        % FORCES
        function obj = plot(obj)
            amp = 3; % Amplification factor
            nodesAmpCoords = zeros(size(obj.NODES));
            
            for i = 1:size(obj.NODES, 1)
                for j = 1:size(obj.NODES, 2)
                    nodesAmpCoords(i, j) = obj.NODES(i, j) + amp*obj.U(obj.DOFS(i, j)); % Amplified coordinates to visualise deformation
                end
            end 
            
            % Plotting
            figure('units', 'normalized', 'outerposition', [0 0 1 1]); hold all; grid on; tol = 1e-3;
            xmin = min(nodesAmpCoords(:,1)); xmax = max(nodesAmpCoords(:,1)); difx = xmax - xmin;
            ymin = min(nodesAmpCoords(:,2)); ymax = max(nodesAmpCoords(:,2)); dify = ymax - ymin; fac = 0.25;
            axis([xmin-difx*fac  xmax+difx*fac  ymin-dify*fac  ymax+dify*fac]);
            
            for EL = 1:obj.elements
                criticalMemberCheck = ismember(EL, obj.indexOfElement);

                n1 = obj.ELEMENTS(EL,1); n2 = obj.ELEMENTS(EL,2); % Identify element node numbers
                x1 = obj.NODES(n1,1); y1 = obj.NODES(n1,2); % Element node 1 - x,y original coordinates
                x2 = obj.NODES(n2,1); y2 = obj.NODES(n2,2); % Element node 2 - x,y original coordinates
                
                % Checking the changes in member lengths and plotting
                % amplified deformed structure
                x1_amp = nodesAmpCoords(n1,1); y1_amp = nodesAmpCoords(n1,2); % Element node 1 - x,y amplified deformed coordinates
                x2_amp = nodesAmpCoords(n2,1); y2_amp = nodesAmpCoords(n2,2); % Element node 2 - x,y amplified deformed coordinates
                
                if criticalMemberCheck == 1 % Check if EL is a critical member
                    plot([x1, x2], [y1, y2], 'Color', [0.5 0.5 0.5], 'Linewidth', 3); % Plotting original critical member shape
                    plot([x1_amp,x2_amp],[y1_amp,y2_amp],'k--','Linewidth', 3); % Plotting amplified deformation of critical member
                    continue
                else
                    if obj.Ftot(EL) < tol % Element length has decreased - member in compression
                        col = 'b'; % Blue colour
                    elseif obj.Ftot(EL) > tol % Element length as increased - member in tension
                        col = 'r'; % Red colour
                    else % No change in element length
                        col = 'k'; % Black colour
                    end
                    
                    plot([x1, x2], [y1, y2], 'Color', [0.5 0.5 0.5], 'Linewidth', 3);
                    plot([x1_amp,x2_amp],[y1_amp,y2_amp],col,'Linewidth',3);
                    
                    % Calculate the axial force in the element, and display it within a text box
                    % Feel free to comment or uncomment this section to
                    % choose to display or not display the text boxes with
                    % the axial forces in each bar element
                    %x_lims = get(gca,'xlim'); xmin = x_lims(1); xmax = x_lims(2);
                    %y_lims = get(gca,'ylim'); ymin = y_lims(1); ymax = y_lims(2);
                    %ax_pos = get(gca,'position'); ax_xpos = ax_pos(1); ax_ypos = ax_pos(2); ax_width = ax_pos(3); ax_height = ax_pos(4);
                    %txt_x_pos = (0.4*x1_amp + 0.6*x2_amp - xmin)/(xmax - xmin) * ax_width + ax_xpos;
                    %txt_y_pos = (0.4*y1_amp + 0.6*y2_amp - ymin)/(ymax - ymin) * ax_height + ax_ypos;
                    %txt_pos = [txt_x_pos txt_y_pos 0 0];
                    %annotation('textbox',txt_pos,'String',[num2str(obj.Ftot(EL)/1000), ' kN'],'FitBoxToText','on','color',col,'FontSize',20,'BackgroundColor','w');
                    
                    % Plotting nodes last
                    plot(x1,y1,'ko','Markersize',7,'MarkerFaceColor','w');
                    plot(x2,y2,'ko','Markersize',7,'MarkerFaceColor','w');
                    plot(x1_amp,y1_amp,'ko','Markersize',7,'MarkerFaceColor','y');
                    plot(x2_amp,y2_amp,'ko','Markersize',7,'MarkerFaceColor','y');
                    xlabel('x [m]'); ylabel('y [m]'); 
                end
            end
        end

        % METHOD 6 - UPDATING GLOBAL FORCE VECTOR
        function obj = update(obj)
            % Take note of the nodes of the failed bar element
            failedNodes = obj.ELEMENTS(obj.indexOfElement(end), :);
            
            % Convert these nodes into their respective dofs
            failedDofs = [obj.DOFS(failedNodes(1), :) obj.DOFS(failedNodes(2), :)];
            
            % Get nodal coordinates of the failed nodes of the failed bar
            % element
            x1 = obj.NODES(failedNodes(1), 1); y1 = obj.NODES(failedNodes(1), 2);
            x2 = obj.NODES(failedNodes(2), 1); y2 = obj.NODES(failedNodes(2), 2);
            alpha = atan2(y2-y1,x2-x1);
            
            % Update the global force vector by replacing the nodes of the
            % failed bar element with the yield force
            if obj.FInd(obj.indexOfElement(end)) == 1 % Tensile axial force failure
                if alpha >= 0 && alpha <= pi/2
                    obj.F(failedDofs(1)) = obj.F(failedDofs(1)) - obj.PsMNA*cos(abs(alpha));
                    obj.F(failedDofs(3)) = obj.F(failedDofs(3)) + obj.PsMNA*cos(abs(alpha));
                    obj.F(failedDofs(2)) = obj.F(failedDofs(2)) - obj.PsMNA*sin(abs(alpha));
                    obj.F(failedDofs(4)) = obj.F(failedDofs(4)) + obj.PsMNA*sin(abs(alpha));
                    
                elseif alpha < 0 && alpha >= -pi/2
                    obj.F(failedDofs(1)) = obj.F(failedDofs(1)) - obj.PsMNA*cos(abs(alpha));
                    obj.F(failedDofs(3)) = obj.F(failedDofs(3)) + obj.PsMNA*cos(abs(alpha));
                    obj.F(failedDofs(2)) = obj.F(failedDofs(2)) + obj.PsMNA*sin(abs(alpha));
                    obj.F(failedDofs(4)) = obj.F(failedDofs(4)) - obj.PsMNA*sin(abs(alpha));
                    
                elseif alpha <= pi && alpha > pi/2
                    obj.F(failedDofs(1)) = obj.F(failedDofs(1)) + obj.PsMNA*cos(abs(alpha));
                    obj.F(failedDofs(3)) = obj.F(failedDofs(3)) - obj.PsMNA*cos(abs(alpha));
                    obj.F(failedDofs(2)) = obj.F(failedDofs(2)) + obj.PsMNA*sin(abs(alpha));
                    obj.F(failedDofs(4)) = obj.F(failedDofs(4)) - obj.PsMNA*sin(abs(alpha));
                    
                elseif alpha < -pi/2 && alpha >= -pi
                    obj.F(failedDofs(1)) = obj.F(failedDofs(1)) + obj.PsMNA*cos(abs(alpha));
                    obj.F(failedDofs(3)) = obj.F(failedDofs(3)) - obj.PsMNA*cos(abs(alpha));
                    obj.F(failedDofs(2)) = obj.F(failedDofs(2)) + obj.PsMNA*sin(abs(alpha));
                    obj.F(failedDofs(4)) = obj.F(failedDofs(4)) - obj.PsMNA*sin(abs(alpha));
                end
                
            elseif obj.FInd(obj.indexOfElement(end)) == -1 % Compressive axial force failure
                if alpha >= 0 && alpha <= pi/2
                    obj.F(failedDofs(1)) = obj.F(failedDofs(1)) + obj.PsMNA*cos(abs(alpha)); 
                    obj.F(failedDofs(3)) = obj.F(failedDofs(3)) - obj.PsMNA*cos(abs(alpha));
                    obj.F(failedDofs(2)) = obj.F(failedDofs(2)) + obj.PsMNA*sin(abs(alpha));
                    obj.F(failedDofs(4)) = obj.F(failedDofs(4)) - obj.PsMNA*sin(abs(alpha));
                    
                elseif alpha < 0 && alpha >= -pi/2
                    obj.F(failedDofs(1)) = obj.F(failedDofs(1)) + obj.PsMNA*cos(abs(alpha));
                    obj.F(failedDofs(3)) = obj.F(failedDofs(3)) - obj.PsMNA*cos(abs(alpha));
                    obj.F(failedDofs(2)) = obj.F(failedDofs(2)) - obj.PsMNA*sin(abs(alpha));
                    obj.F(failedDofs(4)) = obj.F(failedDofs(4)) + obj.PsMNA*sin(abs(alpha));
                    
               elseif alpha <= pi && alpha > pi/2
                    obj.F(failedDofs(1)) = obj.F(failedDofs(1)) - obj.PsMNA*cos(abs(alpha));
                    obj.F(failedDofs(3)) = obj.F(failedDofs(3)) + obj.PsMNA*cos(abs(alpha));
                    obj.F(failedDofs(2)) = obj.F(failedDofs(2)) - obj.PsMNA*sin(abs(alpha));
                    obj.F(failedDofs(4)) = obj.F(failedDofs(4)) + obj.PsMNA*sin(abs(alpha));
                    
                elseif alpha < -pi/2 && alpha >= -pi
                    obj.F(failedDofs(1)) = obj.F(failedDofs(1)) - obj.PsMNA*cos(abs(alpha));
                    obj.F(failedDofs(3)) = obj.F(failedDofs(3)) + obj.PsMNA*cos(abs(alpha));
                    obj.F(failedDofs(2)) = obj.F(failedDofs(2)) - obj.PsMNA*sin(abs(alpha));
                    obj.F(failedDofs(4)) = obj.F(failedDofs(4)) + obj.PsMNA*sin(abs(alpha));
                end
                
            elseif obj.FInd(obj.indexOfElement(end)) == -2 % Buckling failure
                if alpha >= 0 && alpha <= pi/2
                    obj.F(failedDofs(1)) = obj.F(failedDofs(1)) + obj.PsGMNA(obj.indexOfElement(end))*cos(abs(alpha));
                    obj.F(failedDofs(3)) = obj.F(failedDofs(3)) - obj.PsGMNA(obj.indexOfElement(end))*cos(abs(alpha));
                    obj.F(failedDofs(2)) = obj.F(failedDofs(2)) + obj.PsGMNA(obj.indexOfElement(end))*sin(abs(alpha));
                    obj.F(failedDofs(4)) = obj.F(failedDofs(4)) - obj.PsGMNA(obj.indexOfElement(end))*sin(abs(alpha));
                    
                elseif alpha < 0 && alpha >= -pi/2
                    obj.F(failedDofs(1)) = obj.F(failedDofs(1)) + obj.PsGMNA(obj.indexOfElement(end))*cos(abs(alpha));
                    obj.F(failedDofs(3)) = obj.F(failedDofs(3)) - obj.PsGMNA(obj.indexOfElement(end))*cos(abs(alpha));
                    obj.F(failedDofs(2)) = obj.F(failedDofs(2)) - obj.PsGMNA(obj.indexOfElement(end))*sin(abs(alpha));
                    obj.F(failedDofs(4)) = obj.F(failedDofs(4)) + obj.PsGMNA(obj.indexOfElement(end))*sin(abs(alpha));
                    
                elseif alpha <= pi && alpha > pi/2
                    obj.F(failedDofs(1)) = obj.F(failedDofs(1)) - obj.PsGMNA(obj.indexOfElement(end))*cos(abs(alpha));
                    obj.F(failedDofs(3)) = obj.F(failedDofs(3)) + obj.PsGMNA(obj.indexOfElement(end))*cos(abs(alpha));
                    obj.F(failedDofs(2)) = obj.F(failedDofs(2)) - obj.PsGMNA(obj.indexOfElement(end))*sin(abs(alpha));
                    obj.F(failedDofs(4)) = obj.F(failedDofs(4)) - obj.PsGMNA(obj.indexOfElement(end))*sin(abs(alpha));
                    
                elseif alpha < -pi/2 && alpha >= -pi
                    obj.F(failedDofs(1)) = obj.F(failedDofs(1)) - obj.PsGMNA(obj.indexOfElement(end))*cos(abs(alpha));
                    obj.F(failedDofs(3)) = obj.F(failedDofs(3)) + obj.PsGMNA(obj.indexOfElement(end))*cos(abs(alpha));
                    obj.F(failedDofs(2)) = obj.F(failedDofs(2)) - obj.PsGMNA(obj.indexOfElement(end))*sin(abs(alpha));
                    obj.F(failedDofs(4)) = obj.F(failedDofs(4)) + obj.PsGMNA(obj.indexOfElement(end))*sin(abs(alpha));
                end
                
            end
        end
        
        % METHOD 7 - CHECK FOR KFF MATRIX'S STABILITY
        function obj = invertible(obj)
               obj.rcondValue = rcond(obj.KFF);
        end                 
    end          
end

% End of TRUSS.m class file %
    