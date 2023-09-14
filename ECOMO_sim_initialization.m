function [ BoundCond, ModelPara, Options ] = ECOMO_sim_initialization( E )
    %----------------------------------------------------------------------
    % Define the three structures required to execute an ECOMO model
    % simulation
    %
    % [ BoundCond, ModelPara, Options ] = ECOMO_sim_initialization( E );
    %
    % Input Arguments:
    %
    % E             --> (ecomoInterface) object
    % 
    % Output Arguments:
    %
    % BoundCond     --> (struct) ECOMO simulation boundary conditions
    % ModelPara     --> (struct) ECOMO simulation model parameters
    % Options       --> (struct) ECOMO simulation configuration options
    %----------------------------------------------------------------------
    arguments
        E (1,1) ecomoInterface { mustBeNonempty( E ) }
    end
    %----------------------------------------------------------------------
    % Initialise the output structures
    %----------------------------------------------------------------------
    BoundCond = []; 
    Options = []; 
    ModelPara = [];
    %======================================================================
    % Define the boundary conditions
    %======================================================================
    BoundCond.D0 = E.HydDiam * 0.001;                                       % Convert to [m] for simulation
    BoundCond.L = E.TubeLen * 0.001;                                        % Convert to [m] for simulatiom
    BoundCond.N_tube = E.NumTube;                                           % Number of tubes
    BoundCond.ductType = E.DuctGeo;                                         % Duct interior shape
    BoundCond.soot_phi0 = [ 0   0.0; BoundCond.L  0.0 ];                    % Default soot layer thickness [m]
    BoundCond.mRatio_HC_soot = 0;                                           % Initial percentage mass of HCs to soot; range:[0,1]
    BoundCond.N_cell = 10;                                                  % Number of computational cells
    %----------------------------------------------------------------------
    % DEfine HC specie and corresponding ppm levels.
    %----------------------------------------------------------------------
    BoundCond.HCsName_cell = {'C15H32','C17H36','C19H40','C20H42',...
        'C22H46','C24H50'}';
    HCs_ppm = [ 10, 5, 5, 5, 5, 2 ]';
    BoundCond.HCs_frac = HCs_ppm / sum ( HCs_ppm );
    %----------------------------------------------------------------------
    % Define molar fraction, by mass, of water vapour present in the
    % exhaust gas
    %----------------------------------------------------------------------
    BoundCond.X_H2O = obj.X_H2O;                                            % This is a reasonable assumption
    %----------------------------------------------------------------------
    % Define training data
    %----------------------------------------------------------------------
    BoundCond.N_steps = 200;
    BoundCond.IN_TimeSeries = obj.Data;
    BoundCond.dt = median( diff( obj.Data.t ) );                            % Compute time step
    %======================================================================
    % Define configuration options.
    %======================================================================
    SelectMechanism = [                                                     % Comment out any line to disable the mechanism.
        FoulingMechnism.ThermopheresisBT                    
        FoulingMechnism.Diffusiophoresis
        FoulingMechnism.TurbulentImpaction
        FoulingMechnism.GravitationalDrift
        FoulingMechnism.ShearStress
        %FoulingMechnism.SootEmpiricalRemoval                               % Empirical soot removal mechanism disabled
        FoulingMechnism.HCCondEvap
        FoulingMechnism.WaterCondEvap ];
    Options.SelectMechanism = SelectMechanism;
    Options.frictionRegType = 'frictionReg_2';                              % Alternative 'frictionReg_1'
    Options.NuRegType = 'DittusBoelter';                                    % Alternative is 'Gnielinski';
    Options.DeltaType = 'Warey2012';                                        % Alternative is 'He1998';
    Options.SootRemEmpEqType = 'Kuan2017';                                  % Alternative is 'Sul2016'
    Options.TempGradType = 'Warey2012_He1998';                              % Alternative is 'Abarham2013';
    %----------------------------------------------------------------------
    % Disable the sensitivity study perturbations
    %----------------------------------------------------------------------
    Options.SenStudy_1 = false;                                             % Do not edit
    Options.deltaPer = 0.1;                                                 % Do not edit
    %----------------------------------------------------------------------
    % Disable the Sc sweep mechanism
    %----------------------------------------------------------------------
    Options.SweepSc = false;                                                % Do not edit
    BoundCond.DeltaScPer = (-0.8:0.1:0.8);                                  % Do not edit
    %======================================================================
    % Define default settings for the eCoMO model parameters. Note some of
    % these may be overwritten by corresponding DoE values.
    %======================================================================
    c2 = 1; c3 = 1; c4 = 1; c5 = 1;                                         % Model paramters for soot removal empirical regression (assumes mechanism is disabled)
    SootRemEmp = [ c2, c3, c4, c5 ] * 10^-20;
    ModelPara.SootRemEmp = SootRemEmp;
    ModelPara.K_rem = 10^-20;
    ModelPara.NuRegMP =  [  0.5144, 0.5716, 1.1563];                        % Example Dittus-Boelter regression coeffcients
    ModelPara.fRegMP = [ 0.184  -0.2;    0.6488   -0.4757 ; 37.6372   0];   % Example rictionReg_2 regression coeffcients
    %----------------------------------------------------------------------
    % Default lookup table removal gain as a function of EGR mass flow.
    %
    % Column 1 is EGR mass flow.
    % Column 2 is the corresponding removal gain
    %----------------------------------------------------------------------
    ModelPara.K = 0.00001;                                                  % Scale factor use in removal by shear stress, smaller K,smaller removal, more deposited mass
    P_K = [12.66    0.00001
         53.78    0.00001
         100.67   0.008];
    ModelPara.P_K = P_K;
    ModelPara.psi = 1;                                                      % Strength of scale factor, used in shear stress soot removal model
    %----------------------------------------------------------------------
    % Soot particle sticking probability ration threshold
    %----------------------------------------------------------------------
    ModelPara.Ks = 0.8;
    %----------------------------------------------------------------------
    % Scale factor in computing deltaP, inversely proportional so larger
    % value results in a smaller scaling factor
    %---------------------------------------------------------------------- 
    ModelPara.Kp = 50;
    %----------------------------------------------------------------------
    % Scale factor in compute HCs deposit mass, directly proportional so
    % larger value results in a larger scaling factor
    %----------------------------------------------------------------------
    ModelPara.K_HCs = 10^-3;
    %----------------------------------------------------------------------
    % Scale factor used in computing F_V
    %----------------------------------------------------------------------
    ModelPara.K_H = 10;                                                     % Do not edit
    %----------------------------------------------------------------------
    % Scale factor for friction forces (friction factor)
    %----------------------------------------------------------------------
    ModelPara.Kf = 1;                                                       % Do not edit
    %----------------------------------------------------------------------
    % Default spline for distributed k_d (thermal conductivity)
    %----------------------------------------------------------------------
    ModelPara.k_d = [ 0   0.01; BoundCond.L  0.01 ];
    %----------------------------------------------------------------------
    % Default spline for distributed rho_d (deposit layer density)
    %----------------------------------------------------------------------
    ModelPara.rho_d = [ 0   35; BoundCond.L  35];
    %----------------------------------------------------------------------
    % Control the HC evaporation rate, 0 is no HC evaperation
    %----------------------------------------------------------------------
    ModelPara.K_HCevp = 0;
    %----------------------------------------------------------------------
    % Add percentage of residual mSoot, added on 2023-07-06
    %----------------------------------------------------------------------
    ModelPara.PerMsoot = 0.1;
    %----------------------------------------------------------------------
    % Define the deposit soot layer porosity, [0,100) (0 is solid)
    %----------------------------------------------------------------------
    ModelPara.fai = 98;
end % ECOMO_sim_initialization