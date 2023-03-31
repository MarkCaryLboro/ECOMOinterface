classdef ecomoInterface < handle
    % A class to generate the Model parameter structure demanded by the
    % ECOMO code...
    events
        PROCESS_NEW_QUERY                                                   % New point to query available.  
    end

    properties ( SetAccess = protected )
        Src    (1,1)                                                        % Event source
        Lh     (1,1)                                                        % Listener handle for RUN_EXPERIMENT event
        FM     (1,:)   FoulingModel                                         % ECOMO fouling model object array
        IDdata (1,1)   string                                               % Name of identification data file
        B      (1,:)   bayesOpt = bayesOpt.empty                            % bayesOpt object
        Data   (1,1)   struct                                               % Identification training data
    end

    methods
        function obj = loadIdentificationData( obj, Fname, Varname )
            %--------------------------------------------------------------
            % Load the identification data file
            %
            % obj = obj.loadIdentificationData( Fname, Varname );
            %
            % Input Arguments:
            %
            % Fname     --> (string) full name (including path and extension)
            %               of identification data file. If empty, or not found 
            %               a gui is opened allowing the user to select the 
            %               file manually.
            % Varname   --> (string) Name of structure containing
            %               identification data
            %--------------------------------------------------------------
            arguments
                obj     (1,1) ecomoInterface
                Fname   (1,:) string         = string.empty
                Varname (1,:) string         = "EGRVars"
            end
            if ( nargin < 2 ) || ( exist( Fname, "file" ) ~= 2 )
                [ Fname, Path ] = uigetfile( ".mat",...
                        "Select the identification data file",...
                        "MultiSelect", "off");
                Fname = fullfile( Path, Fname );
            end
            Ok = ( exist( Fname, "file" ) == 2 );
            assert( Ok, 'File "%s" not found', Fname );
            obj.IDdata = Fname;
            load( obj.IDdata,  Varname );
            obj.Data = eval( Varname );
        end % loadIdentificationData

        function obj = addRunExperimentListener( obj, Src )
            %--------------------------------------------------------------
            % Add the listener for the RUN_EXPERIMENT event
            % 
            % obj = addRunExperimentListener( Src );
            %
            % Input Arguments:
            %
            % Src --> DoEhook object
            %--------------------------------------------------------------
            arguments 
                obj (1,1) ecomoInterface
                Src (1,1) DoEhook
            end
            obj.Src = Src;
            obj.Lh = addlistener( Src, "RUN_EXPERIMENT",...
                @( SrcObj, Evnt )obj.eventCbRun( SrcObj, Evnt ) );           
        end % ecomoInterface

        function obj = eventCbRun( obj, SrcObj, Evnt )
            %--------------------------------------------------------------
            % Callback for the RUN_EXPERIMENT event
            %
            % 1. Run the ECOMO parameter model script file and generate an
            %    initial ModelPara structure. 
            % 2. Determine if all variables defined in the Parameter table
            %    in the DOE object are defined in the ModelPara structure.
            %    Add any missing fields.
            % 3. Loop through the DoE and run each simulation as requested.
            %
            % obj = obj.eventCbRun( SrcObj, Evnt );
            % 
            % Input Arguments:
            % 
            % SrcObj --> (DoEhook) object generating the RUN_EXPERIMENT
            %            event.
            % Evnt   --> (EventData)
            %--------------------------------------------------------------
            arguments
                obj     (1,1) ecomoInterface
                SrcObj  (1,1) DoEhook
                Evnt    (1,1) event.EventData
            end
            %--------------------------------------------------------------
            % Event check
            %--------------------------------------------------------------
            Ename = string( Evnt.EventName );
            Ok = contains( "RUN_EXPERIMENT", Ename );
            assert( Ok, 'Not processing the %s event supplied', Ename );
            %--------------------------------------------------------------
            % 1. Run the ECOMO model configuration script
            %--------------------------------------------------------------
            run( SrcObj.ConfigFile );
            BoundCond.L = SrcObj.TubeLength;                                % Define length of tube
            BoundCond.D0 = SrcObj.TubeIntDia;                               % Define diameter of tube
            %--------------------------------------------------------------
            % 2. Determine if all variables defined in the Parameter table
            %    in the DOE object are defined in the ModelPara structure.
            %    Add any missing fields.
            %
            % 3. Loop through the data and run the requested simulations
            %--------------------------------------------------------------
            N = height( SrcObj.ParTable );
            for Q = 1:N
                if ~SrcObj.ParTable.Simulated( Q )
                    %------------------------------------------------------
                    % Only run each condition once
                    %------------------------------------------------------
                    ModelPara = obj.parameterCheck( ModelPara,...
                        SrcObj.ParTable, Q );
                    obj = obj.runSimulation( ModelPara, BoundCond,...
                        Options, Q );
                    SrcObj.setSimulated( Q, true );
                end
            end % Q
            %--------------------------------------------------------------
            % Process data and export the results
            %--------------------------------------------------------------
            Res = obj.processResiduals();
            obj = obj.exportData( Res );
        end % eventCbRun

        function plotBestSimulation( obj )
            %--------------------------------------------------------------
            % Plot the best simulation results
            %
            % obj.plotBestSimulation();
            %--------------------------------------------------------------
            D = obj.Src.Design;
            Idx = all( obj.B.Xbest == D, 2 );
            Ptr = find( Idx == 1, 1, "last" );
            FSim = obj.FM( Ptr );                                           % Retrieve the best simulation
            %--------------------------------------------------------------
            % Now plot the results
            %--------------------------------------------------------------
            figure;
            Ax( 4 ) = subplot(2,2,4);
            for Q = 1:4
                Ax( Q ) = subplot(2,2,Q);  
                Ax( Q ).NextPlot = "add";
                grid on;
                Ax( Q ).GridAlpha = 0.5;
                Ax( Q ).GridLineStyle = ":";
                if Q == 3
                    %------------------------------------------------------
                    % Add the identification data
                    %------------------------------------------------------
                    Ax( Q ).NextPlot = "add";
                    plot( obj.Data.T_g_out, obj.Data.deltaP, 'g+');
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % To Do : Generalise the plotting
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            plot( Ax( 1 ), FSim.ModelPara.k_d(:,1), FSim.ModelPara.k_d(:,2), '-' );     % Thermal conductivity
            plot( Ax( 2 ), FSim.ModelPara.rho_d(:,1), FSim.ModelPara.rho_d(:,2), '-');  % Deposit density
            plot( Ax( 3 ), FSim.T_out_L_degC, FSim.deltaPre_L_kPa, 's' );               % Temperature out versus delta pressure
            plot( Ax( 4 ), 1000 * FSim.phi_soot_tn1, '-' );                             % Deposit layer thickness
            xlabel( Ax( 1 ), "Axial Distance [m]");
            ylabel( Ax( 1 ), "Deposit Thermal Conductivity [W/K/m]");
            xlabel( Ax( 2 ), "Axial Distance [m]");
            ylabel( Ax( 2 ), "Deposit Density [kg/m^3]");
            xlabel( Ax( 3 ), "T_{out} [^oC]");
            ylabel( Ax( 3 ), "\Deltap [kPa]")
            legend( Ax( 3 ), "Model", "Data" )
            xlabel( Ax( 4 ), "Control Volume [#]");
            ylabel( Ax( 4 ), "\phi [mm]");
        end % plotBestSimulation

        function obj = genNewQuery( obj )
            %--------------------------------------------------------------
            % Optimise the acquisition function and generate a new query.
            % Augment the design with the new point.
            %
            % obj = obj.genNewQuery();
            %--------------------------------------------------------------
            Lo = 0.5*obj.B.Xlo;
            Hi = 1.5*obj.B.Xhi;
            obj.B = obj.B.acqFcnMaxTemplate( "lb", Lo, "ub", Hi );
        end % genNewQuery

        function exportNewQuery( obj )
            %--------------------------------------------------------------
            % Export the new query and update the design
            %
            % obj.exportNewQuery();
            %--------------------------------------------------------------
            notify( obj, 'PROCESS_NEW_QUERY' );
        end % exportNewQuery

        function L = processResiduals( obj )
            %--------------------------------------------------------------
            % Process the Tout and delta pressure residuals
            %
            % Res = obj.processResiduals();
            %
            % Output Arguments:
            %
            % L --> (double) (Nx1) vector of loglikelihood values.
            %--------------------------------------------------------------
            Tout = [ obj.FM( : ).T_out_L_degC ];                            % Tout predictions from the simulation
            DeltaP = [ obj.FM( : ).deltaPre_L_kPa ];                        % Delta pressure predictions from the simulation
            Tres = ( obj.Data.T_g_out - Tout );                             % Temperature residual matrix
            Pres = ( obj.Data.deltaP - DeltaP );                            % pressure residual matrix
            [ N, C ] = size( Tres );
            %--------------------------------------------------------------
            % Calculate normal negative loglikelihood
            %--------------------------------------------------------------
            L = zeros( C, 1 );
            for Q = 1:C
                L( Q ) = 0.5 * N * ( log( Tres( :,Q ).' * Tres( :,Q ) ) +...
                                log( Pres( :,Q ).' * Pres( :,Q ) ) );
                L( Q ) = N * log( 2 * pi ) + N + L( Q );
            end
            L = -L;
        end % processResiduals        
    end % Ordinary methods

    methods ( Access = protected )
        function obj = exportData( obj, Res, SurModel, AcqFcn )
            %--------------------------------------------------------------
            % Define the bayesOpt object and input the data.
            %
            % obj.exportData();
            %--------------------------------------------------------------
            arguments
                obj      (1,1) ecomoInterface
                Res      (:,1) double
                SurModel (1,1) string = "gpr"
                AcqFcn   (1,1) string = "aei"
            end
            if isempty( obj.B )
                obj.B = bayesOpt( SurModel, AcqFcn );
            end
            X = obj.makeXmatrix();
            obj.B = obj.B.conDataCoding( min( X ),  max( X ) );
            obj.B = obj.B.setTrainingData( X, Res(:) );
        end % exportData

        function obj = runSimulation( obj, ModelPara, BoundCond, Options, Idx )
            %--------------------------------------------------------------
            % Run an ECOMO simulation
            %
            % obj.runSimulation( ModelPara, BoundCond, Options, Idx );
            %
            % Input Arguments:
            %
            % ModelPara --> (struct) Simulation model parameters
            % BoundCond --> (struct) Initial conditions
            % Options   --> (struct) Configuration options
            % Idx       --> (double) Pointer to position in FM array
            %--------------------------------------------------------------
            BoundCond.IN_TimeSeries = obj.Data;
            obj.FM( Idx ) = FoulingModel(BoundCond,ModelPara,Options);      % Define the simulation object
            obj.FM( Idx ).run();                                            % Run the simulation
        end % runSimulation
    end % protected methods

    methods ( Access = private )
        function X = makeXmatrix( obj )
            %--------------------------------------------------------------
            % Construct the X-matrix for the BO process
            %
            % X = obj.makeXmatrix();
            %--------------------------------------------------------------

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % To do: This code only supports identification at a single
            % operating point. Needs generalising for the characterisation
            % experimental data
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            X = obj.Src.Design;
        end  
    end % private methods

    methods ( Access = private, Static = true )
        function P = parameterCheck( M, T, R )
            %--------------------------------------------------------------
            % Return a structure containing all the identification
            % parameter fields. Add to the fields if necessary. Populate
            % fields for identification with row R of the source parameter
            % table.
            % 
            % P = obj.parameterCheck( M, T, R );
            %
            % Input Arguments:
            %
            % M --> (struct) ECOMO model parameter structure
            % T --> (table) List of identification parameter values
            % R --> (double) DoE run to load 
            %--------------------------------------------------------------
            if ( nargin < 3 )
                R = 1;                                                      % apply default
            end
            P = M;
            IdNames = string( T.Properties.VariableNames );
            Idx = ~contains( IdNames, "Simulated" );
            IdNames = IdNames( Idx );
            N = numel( IdNames );
            for Q = 1:N
                %----------------------------------------------------------
                % Add the field to the structure
                %----------------------------------------------------------
                Val = T{ R, IdNames( Q ) };
                if iscell( Val )
                    Val = cell2mat( Val ).';
                end
                %----------------------------------------------------------
                % Make sure the size of the data is correct
                %----------------------------------------------------------
                try
                    %------------------------------------------------------
                    % Match the size as a default
                    %------------------------------------------------------
                    Val = reshape( Val, size( P.( IdNames{ Q } ) ) );
                    P.( IdNames{ Q } ) = Val;
                catch
                    %------------------------------------------------------
                    % Overwrite if required
                    %------------------------------------------------------
                    P.( IdNames{ Q } ) = Val;
                end
            end % Q
        end % parameterCheck
    end % private and static methods
end % classdef