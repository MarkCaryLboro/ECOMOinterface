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
        %----------------------------------------------------------------------
        % Retain only splines with active constraints
        %----------------------------------------------------------------------
        ConIdx = ~cellfun( @isempty, Con );
        %----------------------------------------------------------------------
        % Decode the design
        %----------------------------------------------------------------------
        Theta = D.decodeDesign( Theta );
        for Q = 1:N
            if ConIdx( Q )
                %--------------------------------------------------------------
                % Set the x-points to evaluate the constraint at
                %--------------------------------------------------------------
                Idx = contains( D.Factors.Properties.RowNames, Con{ Q }.name );
                Inc = max( D.Factors.Sz( Idx, : ) );
                X = linspace( 0, D.TubeLength, Inc ).';
                %--------------------------------------------------------------
                % Evaluate the nonlinear constraint
                %--------------------------------------------------------------
                B = D.Bspline{ Con{ Q }.name, "Object" };
                Kidx = D.DesignInfo{ Con{ Q }.name,  "Knots" };
                if iscell( Kidx )
                    Kidx = Kidx{ : };
                end
                Knot = Theta( Kidx );
                Cidx = D.DesignInfo{ Con{ Q }.name,  "Coefficients" };
                if iscell( Cidx )
                    Cidx = Cidx{ : };
                end
                Coef = Theta( Cidx );
                B.n = Knot;
                B.alpha = Coef;
                ApplyConstraint = ~( isempty(Con{ Q }.type) &...
                    isempty(Con{ Q }.derivative) );
                if ApplyConstraint
                    [ C, Ceq ] = B.evalNonlinConstraints( X, Con{ Q } );
                end
            end
        end % /Q
    end
end
