classdef ecomoInterface < handle
    % A class to generate the Model parameter structure demanded by the
    % ECOMO code...

    properties ( SetAccess = protected )
        Src         (1,1) DoEhook                                           % DoE Hook
        FM          (:,1) FoulingModel                                      % ECOMO fouling model object array
        IDdata      (1,1) string                                            % Name of identification data file
        B           (1,1) bayesOpt = bayesOpt( "gpr", "ucb" )               % bayesOpt object
        Data        (1,1) struct                                            % Identification training data
        PNdist      (1,1) struct                                            % PN distribution data
        NumTube     (1,1) double   = 35                                     % Number of heat exchanger tubes
        DuctGeo     (1,1) string   = "Circle"                               % Duct geometry
        TubeLen     (1,1) double   = 185                                    % Tube length in [mm]
        HydDiam     (1,1) double   = 4.5                                    % Hydraulic diameter of the wetted tube [mm]
        CostFcn     (1,1) costFcnType = "Combined"                          % Configure cost function
        NumCells    (1,1) double       = 10                                 % Number of computational cells
        Clean       (1,1) logical      = true                               % Set to logical to indicate clean cooler
    end

    properties ( Constant = true )
        UseParallel (1,1) logical = true
    end % Constants

    properties 
        ConfigFile  (1,1) string                                            % ECOMO model configuration file
        ShowWaitbar (1,1) logical = true
        ExpMult     (1,1) double { mustBeGreaterThanOrEqual( ExpMult, 1 ) } = 1
        SVdTgout    (1,1)   double = 0.0001                                 % Convergence criterion for stopping simulation early
    end

    properties ( Access = public )
        BoundCond   (1,1) struct
        ModelPara   (1,1) struct
    end % private properties

    properties ( Constant = true )
        X_H2O    (1,1)   double = 0.04508                                   % Molar fraction of water vapour present in exhaust gas
    end % Constant properties

    properties ( SetAccess = protected, Dependent = true )
        BestIdx        int64                                                % Pointer to best simulation
        BestFM         FoulingModel                                         % Best simulation object
        Problem        string                                               % Problem type
        InitialSize    double                                               % Size of initial surrogate model training DoE
        NumTsteps      int64                                                % Number of time time steps
        NumPoints      int64                                                % Number of simulation runs
    end % dependent properties

    properties ( Access = private, Dependent = true )
    end

    methods
        function obj = setCleanCooler( obj, State )
            %--------------------------------------------------------------
            % Define the initial state of the cooler as either clean or
            % pre-fouled
            %
            % obj = obj.setCleanCooler( State );
            %
            % Input Arguments:
            %
            % State --> (logical) Set to true (false) to indicate clean
            %                     (pre-fouled) cooler.
            %--------------------------------------------------------------
            arguments
                obj      (1,1) ecomoInterface { mustBeNonempty( obj ) }
                State    (1,1) logical = false
            end
            obj.Clean = State;
        end % setCleanCooler

        function obj = setNumCells( obj, NumCells )
            %--------------------------------------------------------------
            % Set the number of computational cellsused for the simulation
            % runs
            %
            % obj = obj.setNumCells( NumCells );
            %
            % Input Arguments:
            %
            % NumCells --> (int64) Number of computational cells
            %--------------------------------------------------------------
            arguments
                obj      (1,1) ecomoInterface { mustBeNonempty( obj ) }
                NumCells (1,1) double         { mustBePositive( NumCells ) }  = 4.5
            end
            obj.NumCells = NumCells;
        end % setNumCells

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

        function obj = setCostFcn( obj, Fcn )
            %--------------------------------------------------------------
            % Set the cost function type to either "Combined",
            % "Temperature" or "Pressure"
            %
            % obj = obj.setCostFcn( Fcn )
            %
            % Input Arguments:
            %
            % Fcn --> (string) Cost function type: "Combined",
            %                  "Temperature" or "Pressure"
            %--------------------------------------------------------------
            arguments
                obj (1,1) ecomoInterface { mustBeNonempty( obj ) }
                Fcn (1,1) string = "Combined"
            end
            Fcn = costFcnType( Fcn );
            obj.CostFcn = Fcn;
        end % setCostFcn

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

        function obj = initialisation( obj )
            %--------------------------------------------------------------
            % Create the two structures required to configure and run an
            % ECOMO simulation:
            %
            % obj = obj.initialisation()
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
            % Load the configuration function
            %--------------------------------------------------------------
            [ Fname, Path ] = uigetfile( "*.m",...
                "Select ECOMO Simulation Initialisation function",...
                "ECOMO_sim_initialization.m", "MultiSelect","off");
            obj.ConfigFile = fullfile( Path, Fname );
            [ ~, Fname ] = fileparts( obj.ConfigFile );
            Cmd = strjoin( [Fname, "( obj )"], "" );
            %--------------------------------------------------------------
            % Execute the configuration function
            %--------------------------------------------------------------
            [ obj.BoundCond, obj.ModelPara ] = eval( Cmd );
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

        function obj = configureBOPTcoding( obj, Lo, Hi )
            %--------------------------------------------------------------
            % Configure the coding scheme for the Bayesian Optimisation
            % object
            %
            % obj = obj.configureBOPTcoding( Lo, Hi );
            %
            % Input Arguments:
            %
            % Lo  --> (double) Low limits for parameters
            % Hi  --> (double) High limits for parameters
            %--------------------------------------------------------------
            arguments
                obj (1,1) ecomoInterface { mustBeNonempty( obj ) }
                Lo  (1,:) double         { mustBeNonempty( Lo ) }
                Hi  (1,:) double         { mustBeNonempty( Hi ) }
            end
            obj.B = obj.B.conDataCoding( Lo, Hi );
        end % configureBOPTcoding

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
                Varname (1,:) string         = string.empty
            end
            if ( nargin < 2 ) || isempty( Fname ) || ( exist( Fname, "file" ) ~= 2 )
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
        
        function obj = addDoEhook( obj, Src )
            %--------------------------------------------------------------
            % Add the DoEhook object
            % 
            % obj = addDoEhook( Src );
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
        end % addDoEhook

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % NOTE: this method was modified by MathWorks to show an example of
        % how parallel computing can be leveraged to run simulations in
        % parallel in a parfor loop instead of a for loop.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        function obj = doSimulations( obj, SrcObj )
            %--------------------------------------------------------------
            % Perform all outstanding simulations
            %
            % 1. Run the ECOMO parameter model script file and generate an
            %    initial ModelPara structure if required. 
            % 2. Determine if all variables defined in the Parameter table
            %    in the DOE object are defined in the ModelPara structure.
            %    Add any missing fields.
            % 3. Loop through the DoE and run each simulation as requested.
            %
            % obj = obj.doSimulations( SrcObj );
            % 
            % Input Arguments:
            % 
            % SrcObj --> (DoEhook) object. Used to translate design into
            %            corresponding simulation variables
            %--------------------------------------------------------------
            arguments
                obj     (1,1) ecomoInterface
                SrcObj  (1,1) DoEhook                                       = obj.Src
            end
            %--------------------------------------------------------------
            % Composite the DoEhook object
            %--------------------------------------------------------------
            obj = obj.addDoEhook( SrcObj );
            %--------------------------------------------------------------
            % 1. Run the ECOMO model configuration function
            %--------------------------------------------------------------
            if ( exist( obj.ConfigFile, "file" ) ~= 2 )
                obj = obj.initialisation( );
            end
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

            % Modification by M. Cary 11/01/2024
            % Define default configuration structures
            ModPara = obj.ModelPara;                                        % Parameters
            BndCond = obj.BoundCond;                                        % Boundary conditions

            % Simulate new conditions, serially or in parallel.
            if obj.UseParallel 

                % Create DataQueue and listener for waitbar.
                D = parallel.pool.DataQueue;
                afterEach(D, @parforWaitbar);
                parforWaitbar(F, MaxN) % Initialise waitbar.
                %----------------------------------------------------------
                % Run simulations in parallel.
                %----------------------------------------------------------
                parfor Q = 1:Nx
%                 for Q = 1:Nx    % for debugging

                    Idx = N( Q );

                    % Call the updated parameterCheck method.
                    [ModelParaTmp, BoundCondTmp] = ecomoInterface.parameterCheck_Updated( ...
                        ModPara, BndCond, Idx, parTable, type);

                    % Run the Fouling model simulation and assign the
                    % object in the correct index in the Fouling model
                    % array.
                    FMarray{ Q }  = ecomoInterface.runSimulation(ModelParaTmp,...
                        BoundCondTmp );
                    send(D, [] ); % update waitbar
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
                        ModPara, BndCond, SrcObj, Idx);
                    FMarray{Q}  = ecomoInterface.runSimulation(ModelParaTmp,...
                        BoundCondTmp );
                    SrcObj.setSimulated(Idx, true);
                end

            end

            % Convert Fouling model cell array to a normal array.
            FMarray = [FMarray{:}];

            if isempty( obj.FM )
                obj.FM = FMarray;
            else
                obj.FM = [ obj.FM; FMarray ];
            end
            
            delete(F); % close waitbar.

            %--------------------------------------------------------------
            % Process data and export the results
            %--------------------------------------------------------------
            Res = obj.processResiduals();
            obj = obj.exportData( Res );
        end % doSimulations

        function resetFMArray( obj )
            %--------------------------------------------------------------
            % Reset the Fouling Model array in the Ecomointerface object
            %
            % obj.resetFMarray();
            %--------------------------------------------------------------
            obj.FM = FoulingModel.empty;
        end % resetFMArray

        function resetBayes( obj )
            %--------------------------------------------------------------
            % Reset the Bayesian optimisation object in the Ecomointerface 
            % object
            %
            % obj.resetBayes();
            %--------------------------------------------------------------
            AcqFun = obj.B.AcqFcn;
            obj.B = bayesOpt( "gpr", AcqFun );
        end % resetBayes

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
            FSim.plot;
        end % plotFitDiagnostics

        function Ax = plotResiduals( obj, SimNum )
            %--------------------------------------------------------------
            % Residual plots for the selected simulation
            %
            % obj.plotResiduals( SimNum );
            %
            % Input Arguments:
            %
            % SimNum --> (int64) Simulation number to plot. Default is the
            %                    best simulation
            %
            % Output Arguments:
            %
            % Ax     --> (axes) Axes handles to plots
            %--------------------------------------------------------------
            arguments
                obj     (1,1) ecomoInterface
                SimNum  (1,1) int64                                         = obj.BestIdx
            end
            %--------------------------------------------------------------
            % Retrieve the desired simulation
            %--------------------------------------------------------------
            FSim = obj.FM( SimNum );           
            Ax = makeResidualPlots( obj, FSim );
        end % plotResiduals
        
        function plotBayesOpt( obj )
            %--------------------------------------------------------------
            % Plot the Bayesian Optimisation acquisition function results
            %
            % obj.plotBayesOpt();
            %--------------------------------------------------------------
            figure;
            X = 1:obj.InitialSize;
            Ax = axes;
            Ax.NextPlot = "add";
            plot( Ax, X, obj.B.Y( X ), 'bo:', 'LineWidth', 2,...
                                              'MarkerFaceColor', 'blue' );
            X = ( obj.InitialSize + 1 ):obj.Src.NumPoints;
            plot( Ax, X, obj.B.Y( X ), 'rs-', 'LineWidth', 2,...
                                              'MarkerFaceColor', 'red' );
            plot( Ax, obj.BestIdx, obj.B.Y( obj.BestIdx ), 'gh',... 
                                              'MarkerFaceColor', 'green', ...
                                              'MarkerSize', 10 );
            grid on;
            Ax.GridAlpha = 0.5;
            Ax.GridAlphaMode = "Manual";
            Ax.GridLineStyle = "--";
            legend( "Training", "BOpt", "Best", "location", "best",...
                                                "FontSize", 12 );
            if matches( string( obj.Problem ), "Maximisation" )
                Tstr = sprintf( "Bayesian Optimisation Solution Summary - Maximisation" );
            else
                Tstr = sprintf( "Bayesian Optimisation Solution Summary - Minimisation" );
            end
            title( Ax, Tstr, "FontSize", 16);
            xlabel( Ax, "Solution Index", "FontSize", 14);
            ylabel( Ax, "Cost Function", "FontSize", 14);
        end % plotBayesOpt

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
            load( Fname, "BoundCond", "ModelPara" );
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
            F = FoulingModel( BoundCond, ModelPara );
            F.run;
        end

        function FM = exportParameters( obj ) 
            %--------------------------------------------------------------
            % Export the identified parameters to the wotkspace. Store in
            % the Id property.
            %
            % [ FMbest, ModelPara, BoundCond ] = obj.exportParameters();
            %
            % Output Arguments:
            %
            % FM        - (FoulingModel) best simulation for the training
            %--------------------------------------------------------------
            FM = obj.BestFM;
            H = obj.Src;                                                    % DoEhook object
            S = H.DesObj;
            T = S.Factors.Type;
            NumFactors = S.NumFactors;
            Names = H.ParTable.Properties.VariableNames( 1:( end - 1) );
            P = H.ParTable( obj.B.Bidx, 1:( end - 1 ) );
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
                    obj.ModelPara.( Names{ Q } ) = X;
                else
                    obj.BoundCond.( Names{ Q } ) = X;
                end
            end
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
            [ Lo, Hi ] = obj.fetchLimits;
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
            %--------------------------------------------------------------
            % Now run a few iterations of a gradient optimiser using the
            % likelihood as a way of guranteeing improved performance
            %--------------------------------------------------------------
            Xnext = obj.minLogLikelihood();
            %--------------------------------------------------------------
            % Replace the value in the current design and resimulate
            %--------------------------------------------------------------
            obj.B = obj.B.setXnext( Xnext );
        end % genNewQuery

        function exportNewQuery( obj )
            %--------------------------------------------------------------
            % Export the new query and update the design
            %
            % obj.exportNewQuery();
            %--------------------------------------------------------------
            obj.Src.augmentParTable( obj.B.Xnext );
            obj.doSimulations();
        end % exportNewQuery

        function obj = overWrite( obj, Sim )
            %--------------------------------------------------------------
            % Overwrite parameter values in the default ModelPar and
            % BoundCond structures. This permits the results from a
            % previous identification step to be utilised in a subsequent
            % one.
            %
            % obj = obj.overWrite( Sim );
            %
            % Input Arguments:
            %
            % Sim  --> (int64) Simulation number to copy parameters from.
            %          Default is the best simulation.
            %--------------------------------------------------------------
            arguments
                obj  (1,1) ecomoInterface { mustBeNonempty( obj ) }
                Sim  (1,1) int64                                            = obj.BestIdx
            end
            %--------------------------------------------------------------
            % Retrieve the required simulation object
            %--------------------------------------------------------------
            FMsim = obj.FM( Sim );
            obj.BoundCond = FMsim.BoundCond;
            obj.ModelPara = FMsim.ModelPara;
        end % overWrite

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
            N = obj.NumPoints;
            %--------------------------------------------------------------
            % Calculate the pressure and temperature residuals
            %--------------------------------------------------------------
            [DeltaP, Tout ] = deal( zeros( obj.NumTsteps, N ) );
            for Q = 1:N
                Tout( :, Q ) = obj.FM( Q ).T_out_L_degC;                    % Tout predictions from the simulation
                DeltaP( :, Q ) = obj.FM( Q ).deltaPre_L_kPa;                % Delta pressure predictions from the simulation
            end
            Tres = ( obj.Data.T_g_out - Tout );                             % Temperature residual matrix
            Pres = ( obj.Data.deltaP - DeltaP );                            % pressure residual matrix

            [ N, C ] = size( Tres );
            %--------------------------------------------------------------
            % Calculate normal negative loglikelihood
            %--------------------------------------------------------------
            L = zeros( C, 1 );
            for Q = 1:C
                switch obj.CostFcn
                    case "Combined"
                        %--------------------------------------------------
                        % Identify pressure and temperature parameters
                        % together
                        %--------------------------------------------------
                        L( Q ) = 0.5 * N * ( log( Tres( :,Q ).' * Tres( :,Q ) ) +...
                            log( Pres( :,Q ).' * Pres( :,Q ) ) );
                    case "Temperature"
                        %--------------------------------------------------
                        % Identify temerpature parameters only
                        %--------------------------------------------------
                        L( Q ) = 0.5 * N * ( log( Tres( :,Q ).' * Tres( :,Q ) ));
                    case "Pressure"
                        L( Q ) = 0.5 * N * ( log( Pres( :,Q ).' * Pres( :,Q ) ));
                end
                L( Q ) = N * log( 2 * pi ) + N + L( Q );
            end
        end % processResiduals        
    end % Ordinary methods

    methods ( Hidden = true )
        function Ax = makeResidualPlots( obj, FSim )
            %--------------------------------------------------------------
            % Make residual plots for the simulation object provided.
            %
            % Ax = obj.makeResidualPlots( FSim );
            %
            % Input Arguments:
            %
            % FSim --> (FoulingModel) ECOMO simulation object
            %
            % Output Arguments:
            %
            % Ax     --> (axes) Axes handles to plots
            %--------------------------------------------------------------
            figure;
            for  Q = 4:-1:1
                Ax( Q ) = subplot( 2, 2, Q);
                switch Q
                    case 1
                        %--------------------------------------------------
                        % Temperature residual versus time
                        %--------------------------------------------------
                        Res = obj.Data.T_g_out - FSim.T_out_L_degC.';
                        obj.plotResidualsVstime( obj.Data.t, Res, Ax( Q ) );
                        xlabel( "Time [s]" );
                        ylabel( "Outlet Temperature Residual [^oC]" );
                    case 2
                        %--------------------------------------------------
                        % Temperature residual versus predicted
                        %--------------------------------------------------
                        Res = obj.Data.T_g_out - FSim.T_out_L_degC.';
                        Yhat = FSim.T_out_L_degC;
                        obj.plotResidualVsPredicted( Yhat, Res, Ax( Q ));
                        xlabel( "Predicted Outlet Temperature [^oC]");
                        ylabel( "Outlet Temperature Residual [^oC]" );
                    case 3
                        %--------------------------------------------------
                        % Pressure residuals versus time
                        %--------------------------------------------------
                        Res = obj.Data.deltaP - FSim.deltaPre_L_kPa.';
                        obj.plotResidualsVstime( obj.Data.t, Res, Ax( Q ) );
                        xlabel( "Time [s]" );
                        ylabel( "Residual \Deltapressure [kPa]");
                    otherwise
                        %--------------------------------------------------
                        % Pressure residuals versus predicted
                        %--------------------------------------------------
                        Res = obj.Data.deltaP - FSim.deltaPre_L_kPa.';
                        Yhat = FSim.deltaPre_L_kPa;
                        obj.plotResidualVsPredicted( Yhat, Res, Ax( Q ));
                        xlabel( "Predicted \Deltapressure [kPa]");
                        ylabel( "Residual \Deltapressure [kPa]");
                end
                Ax( Q ).GridAlphaMode = "manual";
                Ax( Q ).GridAlpha = 0.75;
                Ax( Q ).GridLineStyle = "-.";
                Ax( Q ).XGrid = "on";
                Ax( Q ).YGrid = "on";
            end
        end % makeResidualPlots
    end % Hidden methods

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

        function S = get.InitialSize( obj )
            % Return the size of the training data pool
            S = obj.Src.DesObj.InitialSize;
        end % get.InitialSize

        function N = get.NumTsteps( obj )
            % Return the number of time steps
            N = numel( obj.Data.pedal );
        end % get.NumTsteps

        function N = get.NumPoints( obj )
            % Return number ofpoints
            N = obj.Src.NumPoints;
        end % get.NumPoints
    end % Get/Set methods

    methods ( Access = protected )
        function Xnext = minLogLikelihood( obj, NumIter )
            %--------------------------------------------------------------
            % A method to directly reduce the likelihood in an attempt to
            % mitigate the concern regarding Bayesian Optimisation not
            % guaranteed to yield an improvement every iteration.
            %
            % Xnext = obj.minLogLikelihood( NumIter )
            %
            % Input Arguments:
            %
            % NumIter --> (int8) Number of iterations to be performed.
            %             NumIter must be in the interval [1, 10]. Default
            %             is 5.
            %--------------------------------------------------------------
            arguments
                obj      (1,1) ecomoInterface     { mustBeNonempty( obj ) }
                NumIter  (1,1) int8 { mustBeGreaterThanOrEqual( NumIter, 1 ), ...
                                      mustBeLessThanOrEqual( NumIter, 20 ) } = 15
            end
            %--------------------------------------------------------------
            % Parse the optional arguments
            %--------------------------------------------------------------
            Names = [ "lb", "ub", "nonlcon", "Aineq", "binq", "Aeq",...
                       "beq", "options" ];
            [ Xlo, Xhi ] = obj.fetchLimits();
            for Q = 1:numel( Names )
                if matches( Names( Q ), "lb" )
                   PROBLEM.( Names( Q ) ) = Xlo;
                elseif matches( Names( Q ), "ub" )
                   PROBLEM.( Names( Q ) ) = Xhi;
                else
                   PROBLEM.( Names( Q ) ) = [];
               end
            end
            %--------------------------------------------------------------
            % Set up the optimisation problem
            %--------------------------------------------------------------
            PROBLEM.options = optimoptions( "fmincon" );
            PROBLEM.options.Display = "none";
            PROBLEM.options.PlotFcn = "optimplotfval";
            PROBLEM.options.MaxIterations = double( NumIter );
            PROBLEM.solver = "fmincon"; 
            PROBLEM.objective = @(X)obj.evalLikelihood( X );
            PROBLEM.x0 = obj.B.Xnext;
            Xnext = fmincon( PROBLEM );
        end % minLogLikelihood

        function L = evalLikelihood( obj, X )
            %--------------------------------------------------------------
            % Evaluate the likelihood function
            %
            % L = obj.evalLikelihood( X);
            %
            % Input Arguments:
            %
            % X --> (double) Parameter vector
            %--------------------------------------------------------------
            L = obj.B.ModelObj.predict( X );
        end % evalLikelihood

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
            [ Lo, Hi ] = obj.fetchLimits();
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
                    TensorFlag = obj.Src.DesObj.Bspline{ Name( Q ), "Tensor" };
                    if TensorFlag
                        %--------------------------------------------------
                        % Tensor product B-spline
                        %--------------------------------------------------
                        Finish = 0;
                        for Kk = 1:numel( S.K )
                            Start = Finish + 1;
                            Finish = Start + S.K( Kk ) - 1;
                            K = Start:Finish;
                            if matches( "x", S.X( Kk ) )
                                %------------------------------------------
                                % Knot limits for tube length are hard
                                % coded for now
                                %------------------------------------------
                                Vlo = 0.10 * obj.TubeLen / 1000;
                                Vhi = 0.90 * obj.TubeLen / 1000;
                            else
                                Vlo = S.Klo( Kk );
                                Vhi = S.Khi( Kk );
                            end
                            Lo( Knots( K ) ) = Vlo;
                            Hi( Knots( K ) ) = Vhi;
                        end
                    else
                        %--------------------------------------------------
                        % One-dimensional B-spline
                        %--------------------------------------------------
                        if matches( "x", S.X)
                            %----------------------------------------------
                            % These are hard coded for now
                            %----------------------------------------------
                            Lo( Knots ) = repmat( 0.10 * obj.TubeLen ...
                                          / 1000, size( Knots ) );
                            Hi( Knots ) = repmat( 0.90 * obj.TubeLen ...
                                          / 1000, size( Knots ) );
                        else
                            Lo( Knots ) = repmat( S.Klo, size( Knots ) );
                            Hi( Knots ) = repmat( S.Khi, size( Knots ) );
                        end
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
        function plotResidualVsPredicted( Yhat, Res, Ax)
            %--------------------------------------------------------------
            % Plot residuals versus the predicted value
            %
            % obj.plotResidualVsPredicted( Yhat, Res, Ax);
            % 
            % Input Arguments:
            %
            % Yhat  --> (double) Predicted vector
            % Res   --> (double) residual vector
            % Ax    --> (axes) Axes to plot on
            %--------------------------------------------------------------
            if ( nargin < 3 )
                figure;
                Ax = axes;
            end
            plot( Ax, Yhat, Res, 'bo', "MarkerFaceColor", "blue" );
        end % plotResidualVsPredicted

        function plotResidualsVstime( T, Res, Ax )
            %--------------------------------------------------------------
            % Plot residuals as a time series
            %
            % obj.plotResidualsVstime( T, Res, Ax );
            %
            % Input Arguments:
            %
            % T     --> (double) Time vector
            % Res   --> (double) residual vector
            % Ax    --> (axes) Axes to plot on
            %--------------------------------------------------------------
            if ( nargin < 3 )
                figure;
                Ax = axes;
            end
            plot( Ax, T, Res, 'bo', "MarkerFaceColor", "blue" );
        end % plotResidualsVstime

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

        function FM = runSimulation( ModelPara, BoundCond )
            %--------------------------------------------------------------
            % Run an ECOMO simulation
            %
            % FM = ecomoInterface.runSimulation( ModelPara, BoundCond );
            %
            % Input Arguments:
            %
            % ModelPara --> (struct) Simulation model 
            % BoundCond --> (struct) Initial conditions
            %--------------------------------------------------------------
            FM = FoulingModel( BoundCond, ModelPara );                      % Define the simulation object
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
            % [ P, B ] = obj.parameterCheck( M, B, R, T, type );
            %
            % Input Arguments:
            %
            % M     --> (struct)  ECOMO model parameter structure
            % B     --> (struct)  ECOMO model boundary conditions structure
            % R     --> (double)  DoE run to load 
            % T     --> (table)   Parameter table containing DoE values
            % type  --> (string)  Either Parameter or Boundary
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
                if matches( IdNames( Q ), "fRegMP" )
                    %------------------------------------------------------
                    % Insert the parameters into the correct place in the
                    % (3x2) array
                    %------------------------------------------------------
                    Freg = P.( IdNames{ Q } );
                    Freg( 2, : ) = Val( 1:2 );
                    Freg( 3, 1 ) = Val( 3 );
                    Val = Freg;
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
%                 Lim = num2cell( L );
                Lim = mat2cell( L, ones( 1, size( L, 1 ) ) );
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