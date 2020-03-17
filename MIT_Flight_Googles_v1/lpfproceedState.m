function [filterState,filterStateDer] = lpfproceedState(input,filterState,filterStateDer,gainP,gainQ,dt_secs)

% -- Determinant:
det = gainP*dt_secs^2 + gainQ*dt_secs + 1;

% -- New filtered state derivative:
stateDer = (filterStateDer + gainP*dt_secs.*input)/det - dt_secs*gainP*filterState/det;

% -- New state:
filterState = dt_secs*(filterStateDer + gainP.*dt_secs.*input)/det + ...
        ((dt_secs*gainQ + 1).*filterState)/det;

% -- Attribute the new filtered state derivative
filterStateDer = stateDer;

end