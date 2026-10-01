function Offspring = DEgenerator(Problem,Population,mlp)

    %% Parameter setting
    [CR,F,proM,disM] = deal(1,0.8,0.2,20);

    Lower = Problem.lower;
	Upper = Problem.upper;

    %% 
    % For each solution in the current population
    for i = 1 : Problem.N
		
		% Choose two different random parents
		p = randperm(Problem.N, 2);
		while p(1)==i || p(2)==i
			p = randperm(Problem.N, 2); 
		end		
       
     	% Generate an child
		Parent1 = Population(i).decs;
		Parent2 = Population(p(1)).decs;
		Parent3 = Population(p(2)).decs;
		[N, D] = size(Parent1);

        %learn the gdv
        [GDV, ~] = mlp.forward(Population(i).decs);

		Site = rand(N,D) < CR; %CR = 1.0
        child = Parent1;

        child(Site) = child(Site) + 0.9*(GDV(Site)-Parent1(Site)) +0.1*(Parent2(Site)-Parent3(Site)); 
		
		%% Polynomial mutation
% 		Site  = rand(N,D) < proM/D;
        Site  = rand(N,D) < proM;
		mu    = rand(N,D);
		temp  = Site & mu<=0.5;
		child       = min(max(child,Lower),Upper);
		child(temp) = child(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
					  (1-(child(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
		temp = Site & mu>0.5; 
		child(temp) = child(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
					  (1-(Upper(temp)-child(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
		
		%Evaluation of the new child
		child = Problem.Evaluation(child);
		
		%add the new child to the offspring population
		Offspring(i) = child;
	end  
end
