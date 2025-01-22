%% manopt n = 1500, p = 60 注意需要在steepestdescent文件中更改停机准则
% localdefaults.tolgradnorm = 1e-4*300;
n = 1500;
p = 60;


% X = orth(randn(n,p));
load("LEP_X_1500_60.mat","X");

beta = 1;
L = gallery('tridiag', n, -1, 2, -1);
[Ll,Lu] = lu(L);


manifold = stiefelfactory(n,p);
problem.M = manifold;

problem.cost = @(X) 0.5*sum(sum(X.*(L*X))) + 1/4*(sum(X.^2, 2)'*(Lu\(Ll\sum(X.^2, 2))));
problem.egrad = @(X) L*X + bsxfun(@times,Lu\(Ll\sum(X.^2, 2)),X);

opts.verbosity = 0;
[x, xcost, info, options] = steepestdescent(problem,X, opts);
out_manopt.time = [info.time];
out_manopt.record = [info.cost];
out_manopt.feasi = [info.gradnorm];


%% manopt n = 2000, p = 80
% localdefaults.tolgradnorm = 1e-4*400;
n = 2000;
p = 80;


% X = orth(randn(n,p));
load("LEP_X_2000_80.mat","X");

beta = 1;
L = gallery('tridiag', n, -1, 2, -1);
[Ll,Lu] = lu(L);


manifold = stiefelfactory(n,p);
problem.M = manifold;

problem.cost = @(X) 0.5*sum(sum(X.*(L*X))) + 1/4*(sum(X.^2, 2)'*(Lu\(Ll\sum(X.^2, 2))));
problem.egrad = @(X) L*X + bsxfun(@times,Lu\(Ll\sum(X.^2, 2)),X);

opts.verbosity = 0;
[x, xcost, info, options] = steepestdescent(problem,X, opts);
out_manopt.time = [info.time];
out_manopt.record = [info.cost];
out_manopt.feasi = [info.gradnorm];

%% manopt n = 1000, p = 40
% localdefaults.tolgradnorm = 1e-4*200;
n = 1000;
p = 40;


% X = orth(randn(n,p));
load("LEP_X_1000_40.mat",'X');

beta = 1;
L = gallery('tridiag', n, -1, 2, -1);
[Ll,Lu] = lu(L);


manifold = stiefelfactory(n,p);
problem.M = manifold;

problem.cost = @(X) 0.5*sum(sum(X.*(L*X))) + 1/4*(sum(X.^2, 2)'*(Lu\(Ll\sum(X.^2, 2))));
problem.egrad = @(X) L*X + bsxfun(@times,Lu\(Ll\sum(X.^2, 2)),X);

opts.verbosity = 0;
[x, xcost, info, options] = steepestdescent(problem,X, opts);
out_manopt.time = [info.time];
out_manopt.record = [info.cost];
out_manopt.feasi = [info.gradnorm];



%%







%%
n = 500;
p = 20;
load('LEP_A_500_500.mat','A');
load('LEP_X_500_20.mat','X');

manifold = stiefelfactory(n,p);
problem.M = manifold;
problem.cost  = @(x) -.5*trace(x'*(A*x));
problem.egrad = @(x) -A*x;
opts.verbosity = 0;
[x, xcost, info, options] = steepestdescent(problem,X, opts);
out_manopt.time = [info.time];
out_manopt.record = [info.cost];
out_manopt.feasi = [info.gradnorm];


%%
n = 1000;
p = 40;
load('LEP_A_1000_1000.mat','A');
load('LEP_X_1000_40.mat','X');
% X = orth(randn(n,p));

manifold = stiefelfactory(n,p);
problem.M = manifold;
problem.cost  = @(x) -.5*trace(x'*(A*x));
problem.egrad = @(x) -A*x;
opts.verbosity = 0;
[x, xcost, info, options] = steepestdescent(problem,X, opts);
out_manopt.time = [info.time];
out_manopt.record = [info.cost];
out_manopt.feasi = [info.gradnorm];


%%
n = 1500;
p = 60;
load('LEP_A_1500_1500.mat','A');
load('LEP_X_1500_60.mat','X');
% X = orth(randn(n,p));

manifold = stiefelfactory(n,p);
problem.M = manifold;
problem.cost  = @(x) -.5*trace(x'*(A*x));
problem.egrad = @(x) -A*x;
opts.verbosity = 0;
[x, xcost, info, options] = steepestdescent(problem,X, opts);
out_manopt.time = [info.time];
out_manopt.record = [info.cost];
out_manopt.feasi = [info.gradnorm];

%%
n = 2000;
p = 80;
% A = randn(n);
% A = A*A';
load('LEP_A_2000_2000.mat','A');
load('LEP_X_2000_80.mat','X');
% X = orth(randn(n,p));

manifold = stiefelfactory(n,p);
problem.M = manifold;
problem.cost  = @(x) -.5*trace(x'*(A*x));
problem.egrad = @(x) -A*x;
opts.verbosity = 0;
[x, xcost, info, options] =  steepestdescent(problem,X, opts);
out_manopt.time = [info.time];
out_manopt.record = [info.cost];
out_manopt.feasi = [info.gradnorm];

%% NEP

%% manopt n = 1000, p = 40
% localdefaults.tolgradnorm = 1e-4*200;
n = 1000;
p = 40;


X = orth(randn(n,p));
% load("LEP_X_1000_40.mat",'X');

beta = 1;
L = gallery('tridiag', n, -1, 2, -1);
[Ll,Lu] = lu(L);


manifold = stiefelfactory(n,p);
problem.M = manifold;

problem.cost = @(X) 0.5*sum(sum(X.*(L*X))) + 1/4*(sum(X.^2, 2)'*(Lu\(Ll\sum(X.^2, 2))));
problem.egrad = @(X) L*X + bsxfun(@times,Lu\(Ll\sum(X.^2, 2)),X);

opts.verbosity = 0;
[x, xcost, info, options] = steepestdescent(problem,X, opts);
out_manopt.time = [info.time];
out_manopt.record = [info.cost];
out_manopt.feasi = [info.gradnorm];

%% manopt n = 1500, p = 60
% localdefaults.tolgradnorm = 1e-4*200;
n = 1500;
p = 60;


X = orth(randn(n,p));
% load("LEP_X_1000_40.mat",'X');

beta = 1;
L = gallery('tridiag', n, -1, 2, -1);
[Ll,Lu] = lu(L);


manifold = stiefelfactory(n,p);
problem.M = manifold;

problem.cost = @(X) 0.5*sum(sum(X.*(L*X))) + 1/4*(sum(X.^2, 2)'*(Lu\(Ll\sum(X.^2, 2))));
problem.egrad = @(X) L*X + bsxfun(@times,Lu\(Ll\sum(X.^2, 2)),X);

opts.verbosity = 0;
[x, xcost, info, options] = steepestdescent(problem,X, opts);
out_manopt.time = [info.time];
out_manopt.record = [info.cost];
out_manopt.feasi = [info.gradnorm];

%% manopt n = 2000, p = 80
% localdefaults.tolgradnorm = 1e-4*200;
n = 2000;
p = 80;


X = orth(randn(n,p));
% load("LEP_X_1000_40.mat",'X');

beta = 1;
L = gallery('tridiag', n, -1, 2, -1);
[Ll,Lu] = lu(L);


manifold = stiefelfactory(n,p);
problem.M = manifold;

problem.cost = @(X) 0.5*sum(sum(X.*(L*X))) + 1/4*(sum(X.^2, 2)'*(Lu\(Ll\sum(X.^2, 2))));
problem.egrad = @(X) L*X + bsxfun(@times,Lu\(Ll\sum(X.^2, 2)),X);

opts.verbosity = 0;
[x, xcost, info, options] = steepestdescent(problem,X, opts);
out_manopt.time = [info.time];
out_manopt.record = [info.cost];
out_manopt.feasi = [info.gradnorm];