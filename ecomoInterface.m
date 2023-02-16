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
        B      (1,1)   bayesOpt                                             % bayesOpt object
        Data   (1,1)   struct                                               % Identification training data
    end

    methods
        function obj = loadIdentificationData( obj, Fname )
            %--------------------------------------------------------------
            % Load the identification data file
            %
            % obj = obj.loadIdentificationData( Fname );
            %
            % Input Arguments:
            %
            % Fname --> (string) full name (including path and extension)
            %           of identification data file. If empty, or not found 
            %           a gui is opened allowing the user to select the 
            %           file manually.
            %--------------------------------------------------------------
            arguments
                obj   (1,1) ecomoInterface
                Fname (1,1) string         = ""
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
            load( obj.IDdata, "EGRVars" );
            obj.Data = EGRVars;
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

        function obj = genNewQuery( obj )
            %--------------------------------------------------------------
            % Optimise the acquisition function and generate a new query.
            % Augment the design with the new point.
            %
            % obj = obj.genNewQuery();
            %--------------------------------------------------------------
            Lo = 0.9 * obj.B.Xlo;
            Hi = 1.1 * obj.B.Xhi;
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
            BO_obj = bayesOpt( SurModel, AcqFcn );
            X = obj.makeXmatrix();
            BO_obj = BO_obj.conDataCoding( min( X ),  max( X ) );
            obj.B = BO_obj.setTrainingData( X, Res(:) );
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
         
        function Res = processResiduals( obj )
            %--------------------------------------------------------------
            % Process the Tout and delta pressure residuals
            %
            % Res = obj.processResiduals();
            %
            % Output Arguments:
            %
            % Res --> (double) Nx2 matrix of residual data. Column 1 is the
            %         temperature residuals and column2 the corresponding 
            %         delta pressure residuals.
            %--------------------------------------------------------------
            Tout = [ obj.FM( : ).T_out_L_degC ];                            % Tout predictions from the simulation
            DeltaP = [ obj.FM( : ).deltaPre_L_kPa ];                        % Delta pressure predictions from the simulation
            Tres = mean( obj.Data.T_g_out - Tout ).';                       % Temperature residual
            Pres = mean( obj.Data.deltaP - DeltaP ).';                      % pressure residual
            Res = ( Tres.^2 + Pres.^2 );                                    % Sum the squared residual. This is the BOpt function to maximise
        end % processResiduals
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
                P.( IdNames{ Q } ) = Val;
            end % Q
        end % parameterCheck
    end % private and static methods
end % classdef