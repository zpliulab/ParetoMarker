classdef ALMOEAAlphaBeta < ALGORITHM
    methods
        function main(Algorithm,Problem)
            [tec,alpha,beta]=Algorithm.ParameterSet(1,0.4,3);
            assert(tec==1); validateattributes(alpha,{'numeric'},{'scalar','nonnegative'});
            validateattributes(beta,{'numeric'},{'scalar','integer','>=',1,'<=',4});
            Population=Problem.Initialization(); population={}; generation=0;
            while Algorithm.NotTerminated(Population,population)
                mlp=ModelLearning(Population);
                Offspring=DEgenerator(Problem,Population,mlp);
                Population=EnvironmentalSelectionNSGA2AlphaBeta([Population,Offspring],Problem.N,alpha,beta);
                assert(numel(Population)==Problem.N);
                generation=generation+1; population{generation}=Population;
            end
        end
    end
end
