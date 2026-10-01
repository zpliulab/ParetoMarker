function [score,auc,pred]=pm_cv(X,y,fold)
% User-supplied folds, training-only feature standardization, explicit positive class.
score=nan(size(y)); pred=nan(size(y));
for f=unique(fold(:))'
    tr=fold~=f; te=fold==f;
    assert(numel(unique(y(tr)))==2 && numel(unique(y(te)))==2);
    mu=mean(X(tr,:),1); sd=std(X(tr,:),0,1); keep=sd>1e-10;
    if ~any(keep)
        score(te)=mean(y(tr)); pred(te)=double(score(te)>=.5); continue;
    end
    train=(X(tr,keep)-mu(keep))./sd(keep);
    test=(X(te,keep)-mu(keep))./sd(keep);
    model=fitcsvm(train,y(tr),'KernelFunction','linear','BoxConstraint',1,'KernelScale',1,'ClassNames',[0;1]);
    model=fitPosterior(model,train,y(tr));
    [~,p]=predict(model,test); positive=find(model.ClassNames==1);
    assert(isscalar(positive)); score(te)=p(:,positive); pred(te)=double(score(te)>=.5);
end
assert(all(isfinite(score)) && all(score>=0 & score<=1));
[~,~,~,auc]=perfcurve(y,score,1);
end
