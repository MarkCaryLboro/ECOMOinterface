classdef ecomoInterface < handle
    % A class to generate the Model parameter structure demanded by the
    % ECOMO code...
    events
        PROCESS_NEW_QUERY                                                   % New point to query available.  
    end

    properties ( SetAccess = protected )
        Src         (1,1)                                                   % Event source
        Lh          (1,1)                                                   % Listener handle for RUN_EXPERIMENT event
        FM          (1,:) FoulingModel                                      % ECOMO fouling model object array
        IDdata      (1,1) string                                            % Name of identification data file
        B           (1,1) bayesOpt = bayesOpt( "gpr", "ucb" )               % bayesOpt object
        Data        (1,1) struct                                            % Identification training data
        PNdist      (1,1) struct                                            % PN distribution data
        NumTube     (1,1) double   = 35                                     % Number of heat exchanger tubes
        DuctGeo     (1,1) string   = "Circle"                               % Duct geometry
        TubeLen     (1,1) double   = 185                                    % Tube length in [mm]
        HydDiam     (1,1) double   = 4.5                                    % Hydraulic diameter of the wetted tube [mm]
    end

    properties 
        ConfigFile  (1,1) string                                            % ECOMO model configuration file
        UseParallel (1,1) logical = true
        ShowWaitbar (1,1) logical = true
    end

    properties ( Constant = true )
        X_H2O   (1,1)   double = 0.04508                                    % Molar fraction of water vapour present in exhaust gas
    end % Constant properties

    properties ( SetAccess = protected, Dependent = true )
        BestIdx        int64                                                % Pointer to best simulation
        BestFM         FoulingModel                                         % Best simulation object
        Problem        string                                               % Problem type
    end % dependent properties

    properties ( Access = private, Dependent = true )
    end

    methods
        function obj = setHydDiam( obj, Dia )
            %--------------------------------------------------------------
            % Set the internal hydraulic diameter of the tube [mm]
            %
            % obj = obj.setHydDiam( Dia );
            %
            % Input Arguments:
            %
            % Dia --> (doubleI Tube inner diameter [mm]
            %--------------------------------------------------------------
            arguments
                obj (1,1) ecomoInterface { mustBeNonempty( obj ) }
                Dia (1,1) double         { mustBePositive( Dia ) }  = 4.5
            end
            obj.HydDiam = Dia; 
        end % setHydDiam

        function obj = setTubeLength( obj, Len )
            %--------------------------------------------------------------
            % Set the length of the heat exchanger tubes
            %
            % obj = obj.setTubeLength( Len );
            %
            % Input Arguments:
            %
            % len   --> (double) lengthr of heat exchanger tubes [mm]
            %--------------------------------------------------------------
            arguments
                obj (1,1) ecomoInterface { mustBeNonempty( obj ) }
                Len (1,1) double         { mustBePositive( Len ) }
            end
            obj.TubeLen = Len;        
        end % setTubeLength

        function obj = setDuctShape( obj, Geo )
            %--------------------------------------------------------------
            % Set the geometric shape of the interior tube geometry. 
            %
            % obj = obj.setDuctShape( Geo );
            %
            % Input Arguments:
            %
            % Geo --> (string) Must be "Circle", "Square" or "Rectangular"
            %--------------------------------------------------------------
            arguments 
                obj (1,1) ecomoInterface { mustBeNonempty( obj ) }
                Geo (1,1) string { mustBeMember( Geo,...
                                   ["Circle", "Square", "Rectangular"]) } = "Circle"
            end
            obj.DuctGeo = ductGeometry( Geo );
        end % setDuctShape

        function obj = setNumberOfTubes( obj, Num )
            %--------------------------------------------------------------
            % Set the number of heat exchanger tubes
            %
            % obj = obj.setNumberOfTubes( Num );
            %
            % Input Arguments:
            %
            % Num   --> (int8) number of heat exchanger tubes
            %--------------------------------------------------------------
            arguments
                obj (1,1) ecomoInterface { mustBeNonempty( obj ) }
                Num (1,1) int8 = 35;
            end
            obj.NumTube = Num;
        end % setNumberOfTubes

        function [ BoundCond, ModelPara, Options ] = initialisation( obj )
            %--------------------------------------------------------------
            % Create the three structures required to configure and run an
            % ECOMO simulation:
            %
            % [ BoundCond, ModelPara, Options ] = obj.initialisation()
            %
            % Notes: the user will be prompted for a MATLAB script defing
            % the setup.
            %--------------------------------------------------------------
            arguments
                obj  (1,1)  ecomoInterface   { mustBeNonempty( obj ) }
            end
            %--------------------------------------------------------------
            % Pre-load the PN distribution data if required
            %--------------------------------------------------------------
            if isempty( obj.PNdist )
                obj = obj.loadPNdistDataFile();                              
            end
            %--------------------------------------------------------------
            % Load the configuration function if required
            %--------------------------------------------------------------
            if ( exist( obj.ConfigFile, 'File' ) ~= 2 )
                [ Fname, Path ] = uigetfile( "*.m",...
                "Select ECOMO Simulation Initialisation function",...
                "ECOMO_sim_initialization.m", "MultiSelect","off");
                obj.ConfigFile = fullfile( Path, Fname );
            end
            [ ~, Fname ] = fileparts( obj.ConfigFile );
            Cmd = strjoin( [Fname, "( obj )"], "" );
            %--------------------------------------------------------------
            % Execute the configuration function
            %--------------------------------------------------------------
            [ BoundCond, ModelPara, Options ] = eval( Cmd );
        end % initialisation

        function obj = defineBayesOpt( obj, Model, AcqFcn )
            %--------------------------------------------------------------
            % Set the bayesOpt object model type and acquisition function
            %
            % Input Arguments:
            %
            % Model     --> (string) surrogate model type. Must be either
            %               {"gpr"} or "rf".
            % AcqFcn    --> (string) Acquisition function name. Must be 
            %               either {"ucb"}, "aei" or "ei".
            %--------------------------------------------------------------
            arguments
                obj     (1,1) ecomoInterface { mustBeNonempty( obj ) }
                Model   (1,1) string = "gpr"
                AcqFcn  (1,1) string = "ucb"
            end
            obj.B = bayesOpt( Model, AcqFcn );
        end % defineBayesOpt

        function obj = loadPNdistDataFile( obj, Fname )
            %--------------------------------------------------------------
            % Load the data file containing the PN distribution particle
            % bin sizes and histogram relative frequencies. Store in
            % property "PNdist".
            %
            % obj = obj.loadPNdistDataFile( Fname );
            %
            % Input Arguments:
            %
            % Fname --> (string) Full file name of PN distirbution values.
            %           If not supplied, the code will prompt the user for
            %           the file via the standard windows file selection
            %           GUI
            %--------------------------------------------------------------
            arguments
                obj     (1,1) ecomoInterface { mustBeNonempty( obj ) }
                Fname   (1,:) string = string.empty( 1, 0)
            end
            if ( nargin < 2 ) || isempty( Fname )
                %----------------------------------------------------------
                % Prompt user to select file
                %----------------------------------------------------------
                [ Fname, Path ] = uigetfile( "*.mat",...
                        "Enter name of file defining PN distribution data",...
                        "PNdata.mat", "MultiSelect", "off");
                Fname = fullfile( Path, Fname );
            end
            Ok = ( exist( Fname, "file" ) == 2 );
            assert( Ok, 'File "%s" does cannot be found', Fname );
            obj.PNdist = load( Fname );
        end % loadPNdistDataFile 

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
        
        function obj = setProblemType( obj, M )
            %--------------------------------------------------------------
            % Set the optimisation problem type (maximisation or 
            % minimisation).
            %
            % obj = obj.setProblemType( M );
            %
            % Input Arguments:
            %
            % M --> (string) set to either "Maximisation" or "Minimisation" 
            %                as appropriate. Default is "Maximisation".
            %--------------------------------------------------------------
            arguments
                obj (1,1) ecomoInterface  { mustBeNonempty( obj ) }
                M   (1,1) string          { mustBeMember( M, ...
                                            ["Maximisation",...
                                             "Minimisation"] ) }            = "Maximisation"
            end
            obj.B = obj.B.setProblemTypeState( M );
        end % setProblemType

        function obj = setXbestAsXnext( obj )
            %--------------------------------------------------------------
            % Set the next query to be evaluated to the best encountered so
            % far.
            %
            % obj.setXbestAsXnext();
            %--------------------------------------------------------------
            obj.B = obj.B.setXbestAsXnext();
        end % setXbestAsXnext

        function [ Lo, Hi ] = setDataBounds( obj, Dx )
            %--------------------------------------------------------------
            % Return the data bounds for the parameters
            %
            % [ Lo, Hi ] = obj.setDataBounds( Dx );
            %
            % Input Arguments:
            %
            % Dx --> (double) Percentage delta to decrease (increase) the
            %        data limits beyond the mimimum (maximum) DoE levels.
            %        Note, 0.05 <= Dx <= 0.95;
            %
            % For example, if min( X ) = [ 1, -1 ] and if Dx = 0.1, then 
            % Lo = [ 0.9 -1.1 ]. Similarly, if max( X ) = [ 1, -1 ] and if 
            % Dx = 0.1, then Hi = [ 1.1, -0.9 ].
            %--------------------------------------------------------------
            Lo = obj.B.Xlo;
            Hi = obj.B.Xhi;
            [ PidxHi, PidxLo ] = obj.getPositive( Lo, Hi );
            %--------------------------------------------------------------
            % Now adjust the bounds
            %--------------------------------------------------------------
            Lo( PidxLo ) = ( 1 - Dx ) * Lo( PidxLo );
            Lo( ~PidxLo ) = ( 1 + Dx ) * Lo( ~PidxLo );
            Hi( PidxHi ) = ( 1 + Dx ) * Hi( PidxHi );
            Hi( ~PidxHi ) = ( 1 - Dx ) * Hi( ~PidxHi );
        end % setDataBounds
        
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
        end % addRunExperimentListener


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % NOTE: this method was modified by MathWorks to show an example of
        % how parallel computing can be leveraged to run simulations in
        % parallel in a parfor loop instead of a for loop.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
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
            % 1. Run the ECOMO model configuration function
            %--------------------------------------------------------------
            [ BoundCond, ModelPara, Options ] = obj.initialisation( );
            %--------------------------------------------------------------
            % 2. Determine if all variables defined in the Parameter table
            %    in the DOE object are defined in the ModelPara structure.
            %    Add any missing fields.
            % 3. Loop through the data and run the requested simulations
            %--------------------------------------------------------------
            N = find( ~SrcObj.ParTable.Simulated, height( SrcObj.ParTable));
            F = waitbar(0, 'DoE Simulation Progress');
            MaxN = max( N );
            Nx = numel( N );
            
            % Initialise array of Fouling model objects.
            FMarray = cell(1, Nx );

            % Extract properties of interest from the DoEhook, so that the
            % we do not need the full DoEhook object within the parfor
            % loop.
            parTable = SrcObj.ParTable;
            type = SrcObj.Type;

            % Simulate new conditions, serially or in parallel.
            if obj.UseParallel 

                % Create DataQueue and listener for waitbar.
                D = parallel.pool.DataQueue;
                afterEach(D, @parforWaitbar);
                parforWaitbar(F, MaxN) % Initialise waitbar.

                % Run simulations in parallel.
                parfor Q = 1:Nx
                    Idx = N( Q );

                    % Call the updated parameterCheck method.
                    [ModelParaTmp, BoundCondTmp] = ecomoInterface.parameterCheck_Updated( ...
                        ModelPara, BoundCond, Idx, parTable, type);

                    % Run the Fouling model simulation and assign the
                    % object in the correct index in the Fouling model
                    % array.
                    FMarray{Q}  = ecomoInterface.runSimulation(ModelParaTmp,...
                        BoundCondTmp, Options);
                    send(D, []); % update waitbar
                end

                % No need to do this within the parfor loop.
                for Q = 1:Nx
                    Idx = N( Q );
                    SrcObj.setSimulated(Idx, true);
                end

            else

                % Run simulations serially.
                for Q = 1:Nx
                    Idx = N( Q );
                    msg = sprintf('Simulation %4.0f out of %4.0f', Idx, MaxN);
                    waitbar(Idx / MaxN, F, msg);

                    [ModelParaTmp, BoundCondTmp] = ecomoInterface.parameterCheck( ...
                        ModelPara, BoundCond, SrcObj, Idx);
                    FMarray{Q}  = ecomoInterface.runSimulation(ModelParaTmp,...
                        BoundCondTmp, Options);
                    SrcObj.setSimulated(Idx, true);
                end

            end

            % Convert Fouling model cell array to a normal array.
            FMarray = [FMarray{:}];

            if isempty( obj.FM )
                obj.FM = FMarray;
            else
                obj.FM = [ obj.FM, FMarray ];
            end
            
            delete(F); % close waitbar.

            %--------------------------------------------------------------
            % Process data and export the results
            %--------------------------------------------------------------
            Res = obj.processResiduals();
            obj = obj.exportData( Res );
        end % eventCbRun

        function plotFitDiagnostics( obj, SimNum )
            %--------------------------------------------------------------
            % Plot the best simulation results
            %
            % obj.plotFitDiagnostics( SimNum );
            %
            % Input Arguments:
            %
            % SimNum --> (int64) Simulation number to plot. Default is the
            %                    best simulation
            %--------------------------------------------------------------
            arguments
                obj     (1,1) ecomoInterface
                SimNum  (1,1) int64                                         = obj.BestIdx
            end
            %--------------------------------------------------------------
            % Retrieve the desired simulation
            %--------------------------------------------------------------
            FSim = obj.FM( SimNum );                                              
            %--------------------------------------------------------------
            % Plot the simulation results
            %--------------------------------------------------------------
            figure;
            FSim.plot;
            %--------------------------------------------------------------
            % Plot the fit diagnostics
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
            %--------------------------------------------------------------
            % Allow for either single value or distributed parameters
            %--------------------------------------------------------------
            X = FSim.ModelPara.k_d(:,1);
            Xi = linspace( min( X ), max( X ), 1001 );
            Y = FSim.ModelPara.k_d(:,2);
            plot( Ax( 1 ), Xi, interp1( X, Y, Xi, 'spline' ), '-',...
                                        'LineWidth', 2.0 );                 % Thermal conductivity
            xlabel( Ax( 1 ), "Axial Distance [m]");
            ylabel( Ax( 1 ), "Deposit Thermal Conductivity [W/K/m]");
            X = FSim.ModelPara.rho_d(:,1);
            Y = FSim.ModelPara.rho_d(:,2);
            plot( Ax( 2 ), Xi, interp1( X, Y, Xi, 'spline' ), '-',...
                                        'LineWidth', 2.0 );                 % Deposit density
            xlabel( Ax( 2 ), "Axial Distance [m]");
            ylabel( Ax( 2 ), "Deposit Density [kg/m^3]");            
            plot( Ax( 3 ), FSim.T_out_L_degC, FSim.deltaPre_L_kPa, 's' );   % Temperature out versus delta pressure
            xlabel( Ax( 3 ), "T_{out} [^oC]");
            ylabel( Ax( 3 ), "\Deltap [kPa]")
            legend( Ax( 3 ), "Data", "Model", "Location", "NorthWest" )
            X = FSim.BoundCond.soot_phi0( :, 1 );
            Y = 1000*FSim.BoundCond.soot_phi0( :, 2 );
            Yend = 1000*FSim.phi_soot_tn1;
            yyaxis( Ax( 4 ), 'left' );
            plot( Ax( 4 ), Xi, interp1( X, Y, Xi, 'spline' ), 'b-',...
                           Xi, interp1( X, Yend, Xi, 'spline' ), 'b:',...
                           'LineWidth', 2.0 );                              % Deposit layer thickness
            xlabel( Ax( 4 ), "Axial Distance [m]");
            ylabel( Ax( 4 ), "\phi [mm]");
            yyaxis( Ax( 4 ), 'right' );
            X = linspace( 0, FSim.BoundCond.L, FSim.BoundCond.N_cell );
            Y = FSim.r_HCs_soot_tn1;
            ylabel( Ax( 4 ), 'HC/Soot ratio' );
            plot( Xi, interp1( X, Y, Xi, 'spline' ), 'LineWidth', 2.0 );
            legend( Ax(4), '\phi(x,0)', '\phi(x,t_{end})', '\Xi', 'Location',...
                'northoutside', 'Orientation', 'horizontal');
        end % plotBestSimulation
        
        function F = runIdentifiedSimulation( obj, R, Cln )
            %--------------------------------------------------------------
            % Run an ECOMO model simulation with the identified parameters.
            % If required the time series can be replicated multiple times.
            %
            % F = obj.runIdentifiedSimulation( R, Cln );
            %
            % Input Arguments:
            %
            % R   --> (int64) Number of times to replicate the time series.
            %                 {0}. For example, setting R = 1 doubles the
            %                 length of the time series. 
            % Cln --> (logical) Set to true to run from clean tube state
            %--------------------------------------------------------------
            arguments
                obj (1,1) ecomoInterface     { mustBeNonempty( obj ) }
                R   (1,1) int64              { mustBePositive( R )}    = 0
                Cln (1,1) logical   = false
            end
            %--------------------------------------------------------------
            % Load the identified structures
            %--------------------------------------------------------------
            [ Fname, Path ] = uigetfile( ".mat",...
                "Select the identification parameter file",...
                "MultiSelect", "off");
            Fname = fullfile( Path, Fname );
            load( Fname, "BoundCond", "ModelPara", "Options" );
            if Cln
                %----------------------------------------------------------
                % Run clean tube simulation if required
                %----------------------------------------------------------
                BoundCond.soot_phi0( :, 2 ) = 0;
            end
            %--------------------------------------------------------------
            % Replicate the time series
            %--------------------------------------------------------------
            if ~isfield( BoundCond, "IN_TimeSeries" )
                BoundCond.IN_TimeSeries = obj.Data;
            end
            for Q = 1:R
                %----------------------------------------------------------
                % Replicate the time series if necessary
                %----------------------------------------------------------
                BoundCond.IN_TimeSeries = ...
                            obj.repTimeSeries( BoundCond.IN_TimeSeries );
            end % /Q
            F = FoulingModel( BoundCond, ModelPara, Options );
            F.run;
        end

        function [ FM, ModelPara, BoundCond, Options ] = exportParameters( obj ) 
            %--------------------------------------------------------------
            % Export the identified parameters to the wotkspace. Store in
            % the Id property.
            %
            % [ FMbest, ModelPara, BoundCond ] = obj.exportParameters();
            %
            % Output Arguments:
            %
            % FMbest    - (FoulingModel) best simulation for the training
            % data
            % ModelPara - (struct) fouling model parameter structure
            % BoundCond - (struct) 
            %--------------------------------------------------------------
            H = obj.Lh.Source{ : };                                         % DoEhook object
            S = H.Lh.Source{ : };
            T = S.Factors.Type;
            NumFactors = S.NumFactors;
            Names = H.ParTable.Properties.VariableNames( 1:( end - 1) );
            P = H.ParTable( obj.B.Bidx, 1:( end - 1 ) );
            %--------------------------------------------------------------
            % Run the configuration file to define the parameters
            %--------------------------------------------------------------
            [ ~, Cmd ] = fileparts( obj.ConfigFile );
            Cmd = strjoin( [ Cmd,"( obj )"], "" );
            [ BoundCond, ModelPara, Options ] = eval( Cmd );
            for Q = 1:NumFactors
                %----------------------------------------------------------
                % Overwrite the parameters and boundary conditions as
                % required
                %----------------------------------------------------------
                X = P{ :, Names{ Q } };
                if iscell( X )
                    X = X{ : };
                end
                if contains( "Parameter", T( Q ) )
                    ModelPara.( Names{ Q } ) = X;
                else
                    BoundCond.( Names{ Q } ) = X;
                end
            end
            %--------------------------------------------------------------
            % Save the results
            %--------------------------------------------------------------
            [ Path, Fname, Ext ] = fileparts( obj.IDdata );
            Fname = strjoin( [ Fname, "Analysis"], "_" );
            Fname = strjoin( [ Fname, Ext ], "" );
            Fname = fullfile( Path, Fname );
            BoundCond.IN_TimeSeries = obj.Data;
            FM = obj.BestFM;
            save( Fname, "FM", "BoundCond", "ModelPara", "Options" );
        end % exportParameters

        function plotTimeSeries( obj, SimNum )
            %--------------------------------------------------------------
            % Plot the model predictions versus the identification data as
            % a time series.
            %
            % obj.plotTimeSeries( SimNum );
            %
            % Input Arguments:
            %
            % SimNum --> (int64) Simulation number to plot. Default is the
            %                    best simulation
            %--------------------------------------------------------------
            arguments
                obj    (1,1) ecomoInterface     { mustBeNonempty( obj ) }
                SimNum (1,1) int64              { mustBePositive( SimNum) } = obj.BestIdx
            end
            %--------------------------------------------------------------
            % Retrieve the simulation
            %--------------------------------------------------------------
            FSim = obj.FM( SimNum );                                        
            %--------------------------------------------------------------
            % Define signals to plot
            %--------------------------------------------------------------
            Torque = obj.Data.torque;
            RPM = obj.Data.rpm;
            Dtemp = FSim.deltaTemp_L_degC;
            Dpress = FSim.deltaPre_L_kPa;
            %--------------------------------------------------------------
            % Define plotting axes
            %--------------------------------------------------------------
            figure;
            for Q = 4:-1:1
                Ax( Q ) = subplot( 2, 2, Q);
                Ax( Q ).NextPlot = "add";
                grid on;
                Ax( Q ).GridAlpha = 0.5;
                Ax( Q ).GridLineStyle = ":";
            end
            %--------------------------------------------------------------
            % Make the plots
            %--------------------------------------------------------------
            Tres = obj.Data.T_g_in - obj.Data.T_g_out;
            plot( Ax( 1 ), obj.Data.t, Tres, 'g-', 'LineWidth', 2.0 );
            plot( Ax( 1 ), obj.Data.t, Dtemp, 'r-', 'LineWidth', 2.0 );
            xlabel( Ax( 1 ), "Time [s]", "FontSize", 14 );
            ylabel( Ax( 1 ), '\Delta Temp [^oc]', "FontSize", 14);
            legend( Ax( 1 ), "Data", "Model", "Location", "northoutside");
            plot( Ax( 2 ), obj.Data.t, obj.Data.deltaP, 'g-', 'LineWidth', 2.0 );
            plot( Ax( 2 ), obj.Data.t, Dpress, 'r-', 'LineWidth', 2.0 );
            xlabel( Ax( 2 ), "Time [s]", "FontSize", 14 );
            ylabel( Ax( 2 ), '\Delta Pressure [kPa]', "FontSize", 14);
            legend( Ax( 2 ), "Data", "Model", "Location", "northoutside");
            plot( Ax( 3 ), obj.Data.t, RPM,  'k-', 'LineWidth', 2.0 );
            xlabel( Ax( 3 ), "Time [s]", "FontSize", 14 );
            ylabel( Ax( 3 ), 'Engine Speed [RPM]', "FontSize", 14);
            plot( Ax( 4 ), obj.Data.t, Torque,  'k-', 'LineWidth', 2.0 );
            xlabel( Ax( 4 ), "Time [s]", "FontSize", 14 );
            ylabel( Ax( 4 ), 'Brake Torque [Nm]', "FontSize", 14);
        end % plotTimeSeries

        function obj = genNewQuery( obj )   
            %--------------------------------------------------------------
            % Optimise the acquisition function and generate a new query.
            % Augment the design with the new point.
            %
            % obj = obj.genNewQuery();
            %--------------------------------------------------------------
            [ Lo, Hi ] = obj.setDataBounds( 0 );
            SobSeqObj = obj.Src.DesObj;
            %--------------------------------------------------------------
            % Define the nonlinear constraints function
            %--------------------------------------------------------------
            NonLinCon = @(X)ecomoBsplineConstraintHandler( X, SobSeqObj );
            %--------------------------------------------------------------
            % Find the next point to query
            %--------------------------------------------------------------
            obj.B = obj.B.acqFcnMaxTemplate( "lb", Lo, "ub", Hi,...
                                             "nonlcon", NonLinCon );
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
%             L = -L;
        end % processResiduals        
    end % Ordinary methods

    methods
        function Idx = get.BestIdx( obj )
            % Return the pointer to the best simulation
            Idx = obj.B.Bidx;
        end % get.BestIdx

        function F = get.BestFM( obj )
            % Return best fouling model simulation
            Ptr = obj.B.Bidx;                                               % Point to the best simulation
            F = obj.FM( Ptr );                                              % Retrieve the best simulation
        end % get.BestFm

        function P = get.Problem( obj )
            % Return the problem type
            P = obj.B.Problem;
        end % get.Problem
    end % Get/Set methods

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
                AcqFcn   (1,1) string = "ucb"
            end
            if isempty( obj.B )
                obj.B = bayesOpt( SurModel, AcqFcn );
            end
            X = obj.makeXmatrix();
            [ Lo, Hi ] = fetchLimits( obj );
            obj.B = obj.B.conDataCoding( Lo, Hi );
            obj.B = obj.B.setTrainingData( X, Res(:) );
        end % exportData
    end % protected methods

    methods ( Access = private ) 
        function [ Lo, Hi ] = fetchLimits( obj )
            %--------------------------------------------------------------
            % Fetch the low and high limits for each factor
            %
            % [ Lo, Hi ] = obj.fetchLimits();
            %--------------------------------------------------------------
            Info = obj.Src.DesObj.DesignInfo;                               % Retrieve pointers to variables
            L = obj.Src.DesObj.Factors.Lo;                                  % Low limit for each factor 
            L = obj.ensureIsCell( L );
            H = obj.Src.DesObj.Factors.Hi;                                  % High limit for each factor
            H = obj.ensureIsCell( H );
            NumColsDesign = size( obj.Src.DesObj.Design, 2 );               % Number of parameters identified
            [ Lo, Hi ] = deal( zeros( 1, NumColsDesign ) );
            Name = string( obj.Src.DesObj.DesignInfo.Properties.RowNames );
            for Q = 1:obj.Src.NumFactors
                %----------------------------------------------------------
                % Generate the limit vectors one factor at a time
                %----------------------------------------------------------
                Coeff = Info{ Q, "Coefficients" };
                if iscell( Coeff )
                    Coeff = Coeff{ : };
                end
                Knots = Info{ Q, "Knots" };
                if iscell( Knots )
                    Knots = Knots{ : };
                end
                Tlo = reshape( L{ Q }, 1, numel( L{ Q } ) );                   
                Lo( Coeff ) = Tlo;
                Thi = reshape( H{ Q }, 1, numel( H{ Q } ) );             
                Hi( Coeff ) = Thi; 
                if ~isnan( Knots )
                    %------------------------------------------------------
                    % Determine knot parameter and limits
                    %------------------------------------------------------
                    S = obj.Src.DesObj.Factors{ Name( Q ), "Spline" };
                    if iscell( S )
                        S = S{ : };
                    end
                    if isempty(S.X) || matches( S.X, "x")
                        Lo( Knots ) = repmat( 0.05 * obj.TubeLen / 1000, size( Knots ) );
                        Hi( Knots ) = repmat( 0.95 * obj.TubeLen / 1000, size( Knots ) );
                    else
                        Lo( Knots ) = repmat( S.Xlo, size( Knots ) );
                        Hi( Knots ) = repmat( S.Xhi, size( Knots ) );
                    end
                end
            end % /Q
        end % fetchLimits

        function X = makeXmatrix( obj )
            %--------------------------------------------------------------
            % Construct the X-matrix for the BO process
            %
            % X = obj.makeXmatrix();
            %--------------------------------------------------------------
            X = obj.Src.Design;
        end  
    end % private methods

    methods ( Access = protected, Static = true )
        function S = repTimeSeries( S )
            %--------------------------------------------------------------
            % One replicate the identification data time series structure
            %
            % S = obj.repTimeSeries( S );
            %
            % Input Arguments:
            %
            % S --> (struct) Storage for replicated data
            %--------------------------------------------------------------
            Fnames = string( fieldnames( S ) );
            N = numel( Fnames );
            for Q = 1:N
                D = S.( Fnames( Q ) );
                D = repmat( D, 2, 1 );
                S.( Fnames( Q ) ) = D;
            end % /Q
        end

        function FM = runSimulation( ModelPara, BoundCond, Options )
            %--------------------------------------------------------------
            % Run an ECOMO simulation
            %
            % FM = ecomoInterface.runSimulation( ModelPara, BoundCond,... 
            %                                    Options );
            %
            % Input Arguments:
            %
            % ModelPara --> (struct) Simulation model 

            % BoundCond --> (struct) Initial conditions
            % Options   --> (struct) Configuration options
            %--------------------------------------------------------------
            FM = FoulingModel(BoundCond,ModelPara,Options);                 % Define the simulation object
            FM.run();                                                       % Run the simulation
        end % runSimulation

        function [ P, Bc ] = parameterCheck( M, B, S, R )
            %--------------------------------------------------------------
            % Return 2 structures containing all the identification
            % parameter and boundary condition fields. Add to the fields if 
            % necessary. Populate fields for identification with row R of 
            % the source parameter table.
            % 
            % [ P, B ] = obj.parameterCheck( M, B, S, R );
            %
            % Input Arguments:
            %
            % M  --> (struct)  ECOMO model parameter structure
            % B  --> (struct)  ECOMO model boundary conditions structure
            % S  --> (DoEhook) Event source object
            % R  --> (double)  DoE run to load 
            %
            % Output Arguments:
            %
            % P  --> (struct) Idenification parameter structure
            % Bc --> (struct) Boundary condition structure
            %--------------------------------------------------------------
            if ( nargin < 3 )
                R = 1;                                                      % apply default
            end
            %--------------------------------------------------------------
            % Define logical parameter vector
            %--------------------------------------------------------------
            Pidx = matches( S.Type, "Parameter");
            %--------------------------------------------------------------
            % Parse the parameter and boundary condition structures
            %--------------------------------------------------------------
            T = S.ParTable;
            P = M;
            Bc = B;
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
                    Val = cell2mat( Val );
                end
                %----------------------------------------------------------
                % Overwrite all parameters
                %----------------------------------------------------------
                if Pidx( Q )
                    %------------------------------------------------------
                    % Identification Parameter
                    %------------------------------------------------------
                    P.( IdNames{ Q } ) = Val;
                else
                    %------------------------------------------------------
                    % Identification Parameter
                    %------------------------------------------------------                    
                    Bc.( IdNames{ Q } ) = Val;
                end
            end % Q
        end % parameterCheck

        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % NOTE: this method was introduced by MathWorks as part of the
        % review, as an example of how the "parameterCheck" method could be
        % modified so that the DoEhook object is no longer a broadcast
        % variable in the parfor loop.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [ P, Bc ] = parameterCheck_Updated( M, B, R, T, type)
            %--------------------------------------------------------------
            % Return 2 structures containing all the identification
            % parameter and boundary condition fields. Add to the fields if 
            % necessary. Populate fields for identification with row R of 
            % the source parameter table.
            % 
            % [ P, B ] = obj.parameterCheck( M, B, S, R );
            %
            % Input Arguments:
            %
            % M  --> (struct)  ECOMO model parameter structure
            % B  --> (struct)  ECOMO model boundary conditions structure
            % S  --> (DoEhook) Event source object
            % R  --> (double)  DoE run to load 
            %
            % Output Arguments:
            %
            % P  --> (struct) Idenification parameter structure
            % Bc --> (struct) Boundary condition structure
            %--------------------------------------------------------------
            if ( nargin < 3 )
                R = 1;                                                      % apply default
            end
            %--------------------------------------------------------------
            % Define logical parameter vector
            %--------------------------------------------------------------
            Pidx = matches( type, "Parameter");
            %--------------------------------------------------------------
            % Parse the parameter and boundary condition structures
            %--------------------------------------------------------------
            P = M;
            Bc = B;
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
                    Val = cell2mat( Val );
                end
                %----------------------------------------------------------
                % Overwrite all parameters
                %----------------------------------------------------------
                if Pidx( Q )
                    %------------------------------------------------------
                    % Identification Parameter
                    %------------------------------------------------------
                    P.( IdNames{ Q } ) = Val;
                else
                    %------------------------------------------------------
                    % Identification Parameter
                    %------------------------------------------------------                    
                    Bc.( IdNames{ Q } ) = Val;
                end
            end % Q
        end % parameterCheck
        
        function [ PidxLo, PidxHi ] = getPositive( Lo, Hi )
            %--------------------------------------------------------------
            % Return logical pointers to establish sign of low and high
            % parameter limits
            %
            % Input Arguments:
            %
            % Lo    --> (double) low parameter limit
            % Hi    --> (double) high parameter limit
            %
            % Output Arguments:
            %
            % PidxLo --> (logical) is true if element of Lo is >= 0
            % PidxHi --> (logical) is true if element of Hi is >= 0
            %--------------------------------------------------------------
            PidxLo = ( Lo >= 0 );
            PidxHi = ( Hi >= 0);
        end % getPositive

        function Lim = ensureIsCell( L )
            %--------------------------------------------------------------
            % Ensure theoutput is a cell array
            %
            % Lim = obj.ensureIsCell( L );
            %
            % Input Arguments:
            %
            % L --> (double) input array
            %--------------------------------------------------------------
            if iscell( L )
                Lim = L;
            else
                Lim = num2cell( L );
            end
        end % ensureIsCell
    end % protected and static methods
end % classdef


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: this helper function was introduced by MathWorks as an example of a
% simple function that displays a waitbar of progress of parallel
% simulations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function parforWaitbar(waitbarHandle, iterations)
    persistent count h N
    
    msg = sprintf('Simulation %4.0f out of %4.0f', count, N);
    
    if nargin == 2
        % Initialize
    
        count = 0;
        h = waitbarHandle;
        N = iterations;
    else
        % Update the waitbar
    
        % Check whether the handle is a reference to a deleted object
        if isvalid(h)
            count = count + 1;
            waitbar(count / N, h, msg);
        end
    end
end