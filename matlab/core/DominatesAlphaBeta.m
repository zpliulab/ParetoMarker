function b = DominatesAlphaBeta(x, y, alpha, beta)
%DOMINATESALPHABETA ParetoMarker dynamic dominance with explicit parameters.
%   This revision-only function is a parameterized copy of the archived
%   Dominates.m. The conventional Pareto branch and the dynamic branch are
%   unchanged. Alpha replaces only the archived literal 0.2, and beta
%   replaces only the archived literal 3.

    if isstruct(x)
        x = x.Cost;
    end
    if isstruct(y)
        y = y.Cost;
    end

    validateattributes(alpha, {'numeric'}, {'scalar','real','finite','nonnegative'});
    validateattributes(beta, {'numeric'}, {'scalar','real','finite','integer','>=',1,'<=',4});
    if numel(x) ~= 5 || numel(y) ~= 5
        error('ParetoMarker:ObjectiveCount', ...
            'Dynamic dominance requires one AUC plus four non-AUC objectives.');
    end

    condition1 = all(x >= y) && any(x > y);
    condition2 = x(1) > y(1);

    x = x(2:end);
    y = y(2:end);
    condition3 = sum(x >= y) >= 1 && any(x > y);
    condition4 = sum(abs(x - y) < alpha .* y) <= beta;

    b = condition1 || (condition2 && condition3 && condition4);
end
