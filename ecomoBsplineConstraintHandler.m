function [ C, Ceq ] = ecomoBsplineConstraintHandler( Theta, D )
    %----------------------------------------------------------------------
    % A function to parse and implement the ECOMO distributed parameter
    % constraints.
    %
    % [ C, Ceq ] = ecomoBsplineConstraintHandler( Theta, D );
    %
    % Input Arguments:
    %
    % Theta --> (double) coded decision variables
    % D     --> (SobolSequence) object containing DoE constraint data
    %
    % Output Arguments:
    %
    % C     --> Nonlinear inequality constraints
    % Ceq   --> Nonlinear equality constraints
    %----------------------------------------------------------------------
    arguments
        Theta (1,:)     double          { mustBeNonempty( Theta ) }
        D     (1,1)     SobolSequence   { mustBeNonempty( D ) }
    end
    %----------------------------------------------------------------------
    % Fetch the constraints
    %----------------------------------------------------------------------
    Ceq = [];
    C = [];
    if D.Constrained
        Con = D.Bspline.Constraint;
        Names = string( D.Bspline.Properties.RowNames );
        N = numel( Names );
        %------------------------------------------------------------------
        % Retain only splines with active constraints
        %------------------------------------------------------------------
        ConIdx = ~cellfun( @isempty, Con );
        for Q = 1:N
            if ConIdx( Q )
                %----------------------------------------------------------
                % Retrieve the corresponding B-spline object
                %----------------------------------------------------------
                B = D.Bspline{ Names( Q ), "Object" };
                if iscell( B )
                    B = B{ : };
                end
                %----------------------------------------------------------
                % Set the x-points to evaluate the constraint at
                %----------------------------------------------------------
                if isa( B, "tensorProductBspline" )
                    %------------------------------------------------------
                    % Two-dimensional tensor product B-spline
                    %------------------------------------------------------
                    X = arrayfun( @linspace, B.A, B.B, [101 101] ,...
                                                   'UniformOutput', false);
                    X = [ X{ 1 } ; X{ 2 } ].';
                else
                    %------------------------------------------------------
                    % One-dimensional B-spline
                    %------------------------------------------------------
                    X = linspace( B.a, B.b, 101 ).';
                end
                %----------------------------------------------------------
                % Evaluate the nonlinear constraint
                %----------------------------------------------------------
                Kidx = D.DesignInfo{ Names( Q ),  "Knots" };
                if iscell( Kidx )
                    Kidx = Kidx{ : };
                end
                Knot = Theta( Kidx );
                Cidx = D.DesignInfo{ Names( Q ),  "Coefficients" };
                if iscell( Cidx )
                    Cidx = Cidx{ : };
                end
                Coef = Theta( Cidx );
                if isa( B, "tensorProductBspline" )
                    %------------------------------------------------------
                    % Two-dimensional tensor product B-spline
                    %------------------------------------------------------
                    B = B.setAlpha( Coef );
                    Seq = B.convertKnotSequences( Knot );
                    B = B.setKnotSequences( Seq );
                else
                    %------------------------------------------------------
                    % One-dimensional B-spline
                    %------------------------------------------------------
                    B.n = Knot;
                    B.alpha = Coef;
                end
                %----------------------------------------------------------
                % Constrained if multiple dimensional structure
                %----------------------------------------------------------
                ApplyConstraint = ( max( size( Con{ Q } ) ) > 1 );
                if ~ApplyConstraint
                    ApplyConstraint = ~( isempty(Con{ Q }.type) &...
                        isempty(Con{ Q }.derivative) );
                end
                if ApplyConstraint
                    [ C, Ceq ] = B.evalNonlinConstraints( X, Con{ Q } );
                end
            end
        end % /Q
    end
end
