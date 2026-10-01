function Population = EnvironmentalSelectionNSGA2AlphaBeta(Population, N, alpha, beta)
%ENVIRONMENTALSELECTIONNSGA2ALPHABETA Archived selection with parameter pass-through.

    result = [];
    for i = 1:length(Population)
        result = [result; Population(i).dec]; %#ok<AGROW>
    end
    [~, uniqueIndices, ~] = unique(result, 'rows', 'stable');
    Population = Population(uniqueIndices);
    result = result(uniqueIndices,:);

    empty_individual.Position = [];
    empty_individual.Cost = [];
    empty_individual.Rank = [];
    empty_individual.DominationSet = [];
    empty_individual.DominatedCount = [];
    empty_individual.CrowdingDistance = [];
    pop = repmat(empty_individual, length(Population), 1);

    for i = 1:length(Population)
        pop(i).Position = Population(i).dec;
        pop(i).Cost = Population(i).obj;
    end

    [pop, F] = NonDominatedSortingAlphaBeta(pop, alpha, beta);
    pop = CalcCrowdingDistance(pop, F);
    [pop, F] = SortPopulation(pop);
    pop = pop(1:min(N,numel(pop)));
    [pop, F] = NonDominatedSortingAlphaBeta(pop, alpha, beta);
    pop = CalcCrowdingDistance(pop, F);
    [pop, F] = SortPopulation(pop);

    F1 = pop(F{1});
    disp([' Number of F1 Members = ' num2str(numel(F1))]);

    cost = [];
    for j = 1:length(pop)
        cost = [cost; pop(j).Position]; %#ok<AGROW>
    end
    matchingIndices = [];
    for i = 1:size(result, 1)
        if any(all(cost == result(i, :), 2))
            matchingIndices = [matchingIndices; i]; %#ok<AGROW>
        end
    end
    Population = Population(matchingIndices);
    if numel(Population)<N
        % Keep N members even if mutation produced too few distinct masks.
        Population = Population(mod(0:N-1,numel(Population))+1);
    end

    resul = [];
    for PI = 1:length(Population)
        resul = [resul; Population(PI).obj]; %#ok<AGROW>
    end
    AUC = resul(:,1);
    outputStr = '****** AUC is ';
    for i = 1:length(AUC)
        outputStr = [outputStr, num2str(AUC(i)), '***']; %#ok<AGROW>
    end
    disp(outputStr);
end
