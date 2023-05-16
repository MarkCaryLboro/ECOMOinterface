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

    properties ( SetAccess = protected, Dependent = true )
        BestFM          FoulingModel
    end % dependent properties

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
            BoundCond.IN_TimeSeries = obj.Data;                             % Load the identification data
            %--------------------------------------------------------------
            % 2. Determine if all variables defined in the Parameter table
            %    in the DOE object are defined in the ModelPara structure.
            %    Add any missing fields.
            % 3. Loop through the data and run the requested simulations
            %--------------------------------------------------------------
            N = find( ~SrcObj.ParTable.Simulated, height( SrcObj.ParTable));
            F = waitbar( 0, 'DoE Simulation Progress' );
            MaxN = max( N );
            Nx = numel( N );
            FMarray = repmat( FoulingModel.empty, 1, Nx );
            for Q = 1:Nx
                %----------------------------------------------------------
                % Simulate new conditions
                %----------------------------------------------------------
                Idx = N( Q );
                Msg = sprintf( 'Simulation %4.0f out of %4.0f', Idx,...
                                                                MaxN );
                waitbar( Idx / MaxN, F, Msg );
                [ ModelParaTmp, BoundCondTmp ] = ecomoInterface.parameterCheck( ...
                    ModelPara, BoundCond, SrcObj, Idx );
                FMarray( Q ) = ecomoInterface.runSimulation( ModelParaTmp,...
                                    BoundCondTmp, Options );
                SrcObj.setSimulated( Idx, true );
            end % Q     
            if isempty( obj.FM )
                obj.FM = FMarray;
            else
                obj.FM = [ obj.FM, FMarray ];
            end
            delete( F );
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
            FSim = obj.BestFM;                                              % Retrieve the best simulation
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
            X = FSim.BoundCond.soot_phi0( 1, : );
            Y = FSim.BoundCond.soot_phi0( 2, : );
            yyaxis( Ax( 4 ), 'left' );
            plot( Ax( 4 ), Xi, interp1( X, Y, Xi, 'spline' ),...
                'LineWidth', 2.0 );                                         % Deposit layer thickness
            xlabel( Ax( 4 ), "Axial Distance [m]");
            ylabel( Ax( 4 ), "\phi [mm]");
            yyaxis( Ax( 4 ), 'right' );
            X = linspace( 0, FSim.BoundCond.L, FSim.BoundCond.N_cell );
            Y = FSim.r_HCs_soot_tn1;
            ylabel( Ax( 4 ), 'HC/Soot ratio' );
            plot( Xi, interp1( X, Y, Xi, 'spline' ), 'LineWidth', 2.0 );
        end % plotBestSimulation

        function [ ModelPara, BoundCond, Options ] = exportParameters( obj ) %#ok<STOUT> 
            %--------------------------------------------------------------
            % Export the identified parameters to the wotkspace. Store in
            % the Id property.
            %
            % [ ModelPara, BoundCond ] = obj.exportParameters();
            %
            % Output Arguments:
            %
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
            run( H.ConfigFile );
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
            save( Fname, "BoundCond", "ModelPara", "Options" );
        end % exportParameters

        function plotTimeSeries( obj )
            %--------------------------------------------------------------
            % Plot the model predictions versus the identification data as
            % a time series.
            %
            % obj.plotTimeSeries();
            %--------------------------------------------------------------
            FSim = obj.BestFM;                                              % Retrieve the best simulation
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
            [ Lo, Hi ] = obj.setDataBounds( 0.1 );
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

    methods
        function F = get.BestFM( obj )
            % Return best fouling model simulation
            Ptr = obj.B.Bidx;                                               % Point to the best simulation
            F = obj.FM( Ptr );                                              % Retrieve the best simulation
        end % get.BestFm
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
            obj.B = obj.B.conDataCoding( min( X ),  max( X ) );
            obj.B = obj.B.setTrainingData( X, Res(:) );
        end % exportData
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

    methods ( Access = protected, Static = true )
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
            Pidx = contains( S.Type, "Parameter");
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
    end % protected and static methods
end % classdef