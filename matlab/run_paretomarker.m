function run_paretomarker(configFile)
% Run corrected Pareto search using user data and a JSON configuration.
root=fileparts(mfilename('fullpath'));oldPath=path;oldDir=pwd;oldRng=rng;
restore=onCleanup(@() restoreState(oldPath,oldDir,oldRng)); %#ok<NASGU>
if nargin<1,configFile=fullfile(root,'..','config','optimization.example.json');end
cfg=jsondecode(fileread(configFile));
% Relative paths are interpreted from the caller's working directory.
if ~isfolder(cfg.input_dir),error('Input directory does not exist.');end
[ok,a]=fileattrib(cfg.input_dir);assert(ok);cfg.input_dir=a.Name;
validateattributes(cfg.population_size,{'numeric'},{'scalar','integer','>=',4});
validateattributes(cfg.max_evaluations,{'numeric'},{'scalar','integer','>=',cfg.population_size});
assert(mod(cfg.max_evaluations,cfg.population_size)==0,'Budget must be a multiple of population size.');
validateattributes(cfg.folds,{'numeric'},{'scalar','integer','>=',2});
validateattributes(cfg.vote_fraction,{'numeric'},{'scalar','>',0,'<=',1});
assert(~isfolder(cfg.output_dir),'Choose a fresh output directory.');mkdir(cfg.output_dir);
[ok,a]=fileattrib(cfg.output_dir);assert(ok);cfg.output_dir=a.Name;
addpath(genpath(fullfile(root,'vendor')));addpath(fullfile(root,'core'),'-begin');
assert(~isempty(which('fitcsvm')),'Statistics and Machine Learning Toolbox is required.');
fid=fopen(fullfile(cfg.output_dir,'config_used.json'),'w');fprintf(fid,'%s',jsonencode(cfg));fclose(fid);
cd(cfg.output_dir);rng(cfg.seed,'twister');
problem=ParetoProblem('maxFE',cfg.max_evaluations,'parameter',{cfg});
alg=ALMOEAAlphaBeta('parameter',{1,cfg.alpha,cfg.beta},'save',100,'outputFcn',@(~,~) []);
t=tic;alg.Solve(problem);search_seconds=toc(t);P=alg.result{end,2};
assert(problem.FE==cfg.max_evaluations);decisions=P.decs;objectives=P.objs;
votes=sum(decisions,1)';selected=votes>=ceil(cfg.vote_fraction*problem.N);
writetable(table(problem.Genes,votes,selected,'VariableNames',{'gene','vote_count','selected'}),'gene_votes.csv');
writematrix(decisions,'final_decisions.csv');
objective_names={'AUC','negative_AIC_like','PAS','ExpressionEntropy','Modularity'};
writetable(array2table(objectives,'VariableNames',objective_names),'final_objectives.csv');
matlab_version=version;products=ver;save('run_info.mat','cfg','search_seconds','matlab_version','products');
fprintf('Finished. Search fitness is descriptive, not independent validation.\n');
end
function restoreState(p,d,r)
path(p);cd(d);rng(r);
end
