function [ C, Ceq ] = ecomoBsplineConstraintHandler( E, D )
    %----------------------------------------------------------------------
    % A function to parse and implement the ECOMO distributed parameter
    % constraints.
    %
    % [ C, Ceq ] = ecomoBsplineConstraintHandler( E, D );
    %
    % Input Arguments:
    %
    % E     --> (ecomoInterface) object
    % D     --> (SobolSequence) object containing DoE data
    %
    % Output Arguments:
    %
    % C     --> Nonlinear inequality constraints
    % Ceq   --> Nonlinear equality constraints
    %----------------------------------------------------------------------
    arguments
        E (1,1)     ecomoInterface  { mustBeNonempty( E ) }
        D (1,1)     SobolSequence   { mustBeNonempty( D ) }
    end
    
end
