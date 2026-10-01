function [pop, F] = NonDominatedSortingAlphaBeta(pop, alpha, beta)
%NONDominatedSortingAlphaBeta Archived sorter with alpha/beta passed through.

    nPop = numel(pop);
    for i = 1:nPop
        pop(i).Rank = []; % reset stale ranks before every sort
        pop(i).DominationSet = [];
        pop(i).DominatedCount = 0;
    end

    F{1} = [];
    for i = 1:nPop
        for j = i+1:nPop
            p = pop(i);
            q = pop(j);

            if DominatesAlphaBeta(p, q, alpha, beta)
                p.DominationSet = [p.DominationSet j];
                q.DominatedCount = q.DominatedCount + 1;
            end
            if DominatesAlphaBeta(q.Cost, p.Cost, alpha, beta)
                q.DominationSet = [q.DominationSet i];
                p.DominatedCount = p.DominatedCount + 1;
            end

            pop(i) = p;
            pop(j) = q;
        end

        if pop(i).DominatedCount == 0
            F{1} = [F{1} i]; %#ok<AGROW>
            pop(i).Rank = 1;
        end
    end

    k = 1;
    while true
        Q = [];
        result = [];
        for mi = 1:length(pop)
            result = [result; pop(mi).Rank]; %#ok<AGROW>
        end
        if length(result) == length(pop)
            break;
        end

        for i = F{k}
            p = pop(i);
            for j = p.DominationSet
                q = pop(j);
                q.DominatedCount = q.DominatedCount - 1;
                if q.DominatedCount == 0
                    Q = [Q j]; %#ok<AGROW>
                    q.Rank = k + 1;
                end
                pop(j) = q;
            end
        end

        if isempty(Q)
            cost = [];
            for mj = 1:length(pop)
                cost = [cost; pop(mj).DominatedCount]; %#ok<AGROW>
            end
            if isempty(cost(cost ~= 0))
                break;
            else
                indices = find(cost == min(cost(cost ~= 0)));
                Q = [Q indices']; %#ok<AGROW>
                if isempty(cost == min(cost(cost ~= 0)))
                    for mi = 1:length(Q)
                        pop(Q(mi)).Rank = 1;
                    end
                else
                    for mi = 1:length(Q)
                        pop(Q(mi)).Rank = k + 1;
                    end
                end
            end
        end

        F{k+1} = Q; %#ok<AGROW>
        k = k + 1;
    end
end
