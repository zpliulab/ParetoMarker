classdef ParetoProblem < PROBLEM
% Revision implementation. Search fitness is NOT unbiased predictive validation.
    properties
        Config; Data; Labels; Genes; Samples; Network; Fold; Cache;
        CacheHits=0; ComputedMasks=0;
    end
    methods
        function Setting(obj)
            cfg=obj.ParameterSet(struct());
            obj.Config=cfg; obj.M=5; obj.N=cfg.population_size;
            raw=readcell(fullfile(cfg.input_dir,'expression.csv'));
            obj.Samples=string(raw(2:end,1));obj.Genes=string(raw(1,2:end))';
            obj.Data=cell2mat(raw(2:end,2:end));obj.D=size(obj.Data,2);
            C=readtable(fullfile(cfg.input_dir,'labels.csv'),'TextType','string');
            assert(numel(unique(C.sample_id))==height(C),'Duplicate label IDs.');
            [found,loc]=ismember(obj.Samples,C.sample_id);assert(all(found),'Unmatched sample IDs.');
            obj.Labels=double(C.label(loc));
            assert(all(ismember(obj.Labels,[0 1])) && numel(unique(obj.Labels))==2,'Labels must be 0/1.');
            assert(numel(unique(obj.Samples))==size(obj.Data,1) && numel(unique(obj.Genes))==obj.D,'Duplicate IDs.');
            assert(obj.D>=2 && all(isfinite(obj.Data(:))),'Finite numeric expression and >=2 genes required.');
            E=readtable(fullfile(cfg.input_dir,'network.csv'),'TextType','string');
            obj.Network=[cellstr(E.source),repmat({''},height(E),1),cellstr(E.target)];
            obj.Cache=containers.Map('KeyType','char','ValueType','any');
            obj.lower=zeros(1,obj.D);obj.upper=ones(1,obj.D);obj.encoding=4*ones(1,obj.D);
            obj.Fold=zeros(size(obj.Labels));stream=RandStream('mt19937ar','Seed',cfg.fold_seed);
            for label=0:1
                idx=find(obj.Labels==label);assert(numel(idx)>=cfg.folds,'Too few samples per class.');
                idx=idx(randperm(stream,numel(idx)));obj.Fold(idx)=mod(0:numel(idx)-1,cfg.folds)+1;
            end
        end
        function Population=Initialization(obj,N)
            % Initialize the evolutionary population from the source-centered
            % modules exported by python/flow.py.  No random fallback is used:
            % a mismatch is treated as an input/configuration error.
            if nargin < 2
                N=obj.N;
            end
            initFile=fullfile(obj.Config.input_dir,'initial_population.csv');
            assert(isfile(initFile), ...
                'Missing initial_population.csv. Run python/flow.py before MATLAB optimization.');
            raw=readcell(initFile);
            assert(size(raw,1)-1==N, ...
                'Initial module count (%d) must equal population_size (%d).',size(raw,1)-1,N);
            initGenes=string(raw(1,2:end))';
            assert(numel(unique(initGenes))==numel(initGenes),'Duplicate genes in initial population header.');
            [found,loc]=ismember(obj.Genes,initGenes);
            assert(all(found) && numel(initGenes)==numel(obj.Genes), ...
                'Initial-population genes must exactly match expression.csv candidate genes.');
            X=cell2mat(raw(2:end,2:end));
            X=X(:,loc);
            assert(all(ismember(X(:),[0 1])),'Initial population must contain binary masks.');
            assert(all(sum(X,2)>0),'Every initial module must contain at least one candidate gene.');
            Population=obj.Evaluation(double(X));
        end

        function F=CalObj(obj,X)
            F=zeros(size(X,1),5);
            for i=1:size(X,1)
                mask=logical(X(i,:)); key=char(48+mask);
                if isKey(obj.Cache,key), F(i,:)=obj.Cache(key); obj.CacheHits=obj.CacheHits+1; continue; end
                obj.ComputedMasks=obj.ComputedMasks+1;
                if ~any(mask)
                    F(i,:)=[-1,-1e12,-1e12,-1e12,-1e12]; continue;
                end
                % Isolate evaluation/calibration/clustering randomness from the search.
                prior=rng; guard=onCleanup(@() rng(prior));
                rng(obj.Config.fold_seed,'twister');
                [score,auc]=pm_cv(obj.Data(:,mask),obj.Labels,obj.Fold);
                n=numel(obj.Labels); k=sum(mask);
                % Negative residual-based AIC-like criterion evaluated on pooled OOF scores.
                aicSurrogate=-(2*k+n*log(max(sum((obj.Labels-score).^2)/n,eps)));
                data=obj.Data(:,mask)';
                H=0; MI=0;
                for j=1:size(data,1)
                    [~,~,bin]=unique(data(j,:)); counts=accumarray(bin(:),1);
                    p=counts/sum(counts); H=H-sum(p.*log2(p));
                    MI=MI+pm_mutual(data(j,:)',obj.Labels,n);
                end
                genes=obj.Genes(mask); adj=zeros(k);
                [ok1,a]=ismember(string(obj.Network(:,1)),genes);
                [ok2,b]=ismember(string(obj.Network(:,3)),genes); valid=ok1 & ok2;
                adj(sub2ind([k k],a(valid),b(valid)))=1;
                edges=nnz(adj); Q=0;
                if edges>0
                    clusters=min([3,k,size(unique(adj,'rows'),1)]);
                    if clusters==1, community=ones(k,1);
                    else, community=kmeans(adj,clusters,'Replicates',1,'MaxIter',100,'EmptyAction','singleton'); end
                    degrees=sum(adj,2)+sum(adj,1)';
                    for c=unique(community)'
                        members=community==c;
                        Q=Q+sum(sum(adj(members,members)))/edges-(sum(degrees(members))/(2*edges))^2;
                    end
                end
                F(i,:)=[auc,aicSurrogate,MI/k,H,Q];
                assert(all(isfinite(F(i,:))),'Nonfinite objective: stop rather than silently omit.');
                obj.Cache(key)=F(i,:);
                clear guard;
            end
        end
    end
end
function MI=pm_mutual(x,y,n)
    if numel(unique(x))<2 || numel(unique(y))<2, MI=0; return; end
    bins=zeros(numel(y),2); values=[x,y];
    for j=1:2
        e=linspace(min(values(:,j)),max(values(:,j)),n+1);
        [~,bins(:,j)]=histc(values(:,j),[-Inf,e(2:end-1),Inf]);
    end
    joint=accumarray(bins,1,[n n])/numel(y);
    px=sum(joint,2); py=sum(joint,1);
    MI=-sum(px.*log2(px+eps))-sum(py.*log2(py+eps))+sum(joint(:).*log2(joint(:)+eps));
end
