function test_DominatesAlphaBeta_formula()
%TEST_DOMINATESALPHABETA_FORMULA Executable MATLAB equivalence checks.

    rng(20260826, 'twister');
    n = 100000;
    X = 2 .* rand(n, 5);
    Y = 2 .* rand(n, 5);

    for k = 1:n
        observed = DominatesAlphaBeta(X(k,:), Y(k,:), 0.4, 3);
        conventional = all(X(k,:) >= Y(k,:)) && any(X(k,:) > Y(k,:));
        xNonAUC = X(k,2:end);
        yNonAUC = Y(k,2:end);
        formalDynamic = X(k,1) > Y(k,1) && ...
            sum(xNonAUC >= yNonAUC) >= 1 && any(xNonAUC > yNonAUC) && ...
            sum(abs(xNonAUC - yNonAUC) < 0.4 .* yNonAUC) <= 3;
        expected = conventional || formalDynamic;
        assert(observed == expected, 'Formula equivalence failed at pair %d.', k);
    end

    fprintf('PASS: %d alpha=0.4, beta=3 dominance pairs matched the formal rule.\n', n);
end
