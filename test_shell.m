%%
foldername = 'all_result_final';
if ~isfolder(foldername)
    mkdir(foldername)
end

tablename = [foldername,'/result_LEP','.xls'];

title = {'d','fval','normG','itr','time/s','feasi:X^TX-I'};

%% LEP : n=500,p=20
n = 500;
p = 20;
load('LEP_A_500_500.mat','A');
% load('LEP_X_500_20.mat','X');

row_sheet = 2;

writecell(title,tablename,'Sheet',sprintf('n=%d,p=%d',n,p),'Range','A1:Z1');

repnum = 10; % 重复试验次数

for d = [2,4,6,8,10,12,14,16,18,20]
    res = zeros(1,5);
    for repeat = 1:repnum
        X = orth(randn(n,p));
        [out] = test_LEP(X,n,p,d,1e-4,A,'rand_without_replacement');
        res_temp = [out.fval,out.normG,out.itr,out.time,norm(out.X'*out.X-eye(p))];
        res = res+res_temp;
    end
    res = res/repnum;
    writematrix([d,res],tablename,'Sheet',sprintf('n=%d,p=%d',n,p), ...
            'Range',sprintf('A%d:Z%d',row_sheet,row_sheet));
    row_sheet = row_sheet+1;
end

%% LEP : n=1000,p=40
n = 1000;
p = 40;
load('LEP_A_1000_1000.mat','A');
load('LEP_X_1000_40.mat','X');

row_sheet = 2;

writecell(title,tablename,'Sheet',sprintf('n=%d,p=%d',n,p),'Range','A1:Z1');

repnum = 10; % 重复试验次数

for d = [4,8,12,16,20,24,28,32,36,40]
    res = zeros(1,5);
    for repeat = 1:repnum
        % X = orth(randn(n,p));
        [out] = test_LEP(X,n,p,d,1e-4,A,'rand_without_replacement');
        res_temp = [out.fval,out.normG,out.itr,out.time,norm(out.X'*out.X-eye(p))];
        res = res+res_temp;
    end
    res = res/repnum;
    writematrix([d,res],tablename,'Sheet',sprintf('n=%d,p=%d',n,p), ...
            'Range',sprintf('A%d:Z%d',row_sheet,row_sheet));
    row_sheet = row_sheet+1;
end

%% LEP : n=1500,p=60
n = 1500;
p = 60;
load('LEP_A_1500_1500.mat','A');
load('LEP_X_1500_60.mat','X');

row_sheet = 2;

writecell(title,tablename,'Sheet',sprintf('n=%d,p=%d',n,p),'Range','A1:Z1');

repnum = 10; % 重复试验次数

for d = [6,12,18,24,30,36,42,48,54,60]
    res = zeros(1,5);
    for repeat = 1:repnum
        X = orth(randn(n,p));
        [out] = test_LEP(X,n,p,d,1e-4,A,'rand_without_replacement');
        res_temp = [out.fval,out.normG,out.itr,out.time,norm(out.X'*out.X-eye(p))];
        res = res+res_temp;
    end
    res = res/repnum;
    writematrix([d,res],tablename,'Sheet',sprintf('n=%d,p=%d',n,p), ...
            'Range',sprintf('A%d:Z%d',row_sheet,row_sheet));
    row_sheet = row_sheet+1;
end

%% LEP : n=2000,p=80
n = 2000;
p = 80;
load('LEP_A_2000_2000.mat','A');
load('LEP_X_2000_80.mat','X');

row_sheet = 2;

writecell(title,tablename,'Sheet',sprintf('n=%d,p=%d',n,p),'Range','A1:Z1');

repnum = 10; % 重复试验次数

for d = [8,16,24,32,40,48,56,64,72,80]
    res = zeros(1,5);
    for repeat = 1:repnum
        % X = orth(randn(n,p));
        [out] = test_LEP(X,n,p,d,1e-4,A,'rand_with_replacement');
        res_temp = [out.fval,out.normG,out.itr,out.time,norm(out.X'*out.X-eye(p))];
        res = res+res_temp;
    end
    res = res/repnum;
    writematrix([d,res],tablename,'Sheet',sprintf('n=%d,p=%d',n,p), ...
            'Range',sprintf('A%d:Z%d',row_sheet,row_sheet));
    row_sheet = row_sheet+1;
end



%% test_NEP

%%
foldername = 'all_result_final';
if ~isfolder(foldername)
    mkdir(foldername)
end

tablename = [foldername,'/result_NEP','.xls'];

title = {'d','fval','normG','itr','time/s','feasi:X^TX-I'};

%% NEP: n=1000,p=40
n = 1000;
p = 40;

% load('LEP_X_1000_40.mat','X');

row_sheet = 2;

writecell(title,tablename,'Sheet',sprintf('n=%d,p=%d',n,p),'Range','A1:Z1');

repnum = 10; % 重复试验次数

for d = [4,8,12,16,20,24,28,32,36,40]
    res = zeros(1,5);
    for repeat = 1:repnum
        X = orth(randn(n,p));
        [out] = test_NEP(X,n,p,d,1e-4,'rand_with_replacement');
        res_temp = [out.fval,out.feasi,out.itr,out.time,norm(out.X'*out.X-eye(p))];
        res = res+res_temp;
    end
    res = res/repnum;
    writematrix([d,res],tablename,'Sheet',sprintf('n=%d,p=%d',n,p), ...
            'Range',sprintf('A%d:Z%d',row_sheet,row_sheet));
    row_sheet = row_sheet+1;
end

%% NEP: n=1500,p=60
n = 1500;
p = 60;

% load('LEP_X_1000_40.mat','X');

row_sheet = 2;

writecell(title,tablename,'Sheet',sprintf('n=%d,p=%d',n,p),'Range','A1:Z1');

repnum = 10; % 重复试验次数

for d = [6,12,18,24,30,36,42,48,54,60]
    res = zeros(1,5);
    for repeat = 1:repnum
        X = orth(randn(n,p));
        [out] = test_NEP(X,n,p,d,1e-4,'rand_with_replacement');
        res_temp = [out.fval,out.feasi,out.itr,out.time,norm(out.X'*out.X-eye(p))];
        res = res+res_temp;
    end
    res = res/repnum;
    writematrix([d,res],tablename,'Sheet',sprintf('n=%d,p=%d',n,p), ...
            'Range',sprintf('A%d:Z%d',row_sheet,row_sheet));
    row_sheet = row_sheet+1;
end

%% NEP: n=2000,p=80
n = 2000;
p = 80;

% load('LEP_X_1000_40.mat','X');

row_sheet = 2;

writecell(title,tablename,'Sheet',sprintf('n=%d,p=%d',n,p),'Range','A1:Z1');

repnum = 10; % 重复试验次数

for d = [8,16,24,32,40,48,56,64,72,80]
    res = zeros(1,5);
    for repeat = 1:repnum
        X = orth(randn(n,p));
        [out] = test_NEP(X,n,p,d,1e-4,'rand_without_replacement');
        res_temp = [out.fval,out.feasi,out.itr,out.time,norm(out.X'*out.X-eye(p))];
        res = res+res_temp;
    end
    res = res/repnum;
    writematrix([d,res],tablename,'Sheet',sprintf('n=%d,p=%d',n,p), ...
            'Range',sprintf('A%d:Z%d',row_sheet,row_sheet));
    row_sheet = row_sheet+1;
end

%% test_OPP

%%
foldername = 'all_result_final';
if ~isfolder(foldername)
    mkdir(foldername)
end

tablename = [foldername,'/result_OPP','.xls'];

title = {'d','fval','normG','itr','time/s','feasi:X^TX-I'};

%% OPP: n=1000,p=40
n = 1000;
p = 40;

load('OPP_A_1000_1000.mat','A');
load('OPP_C_1000_40.mat','C');
% load('OPP_X_1000_40.mat','X');

row_sheet = 2;

writecell(title,tablename,'Sheet',sprintf('n=%d,p=%d',n,p),'Range','A1:Z1');

repnum = 10; % 重复试验次数

for d = [4,8,12,16,20,24,28,32,36,40]
    res = zeros(1,5);
    for repeat = 1:repnum
        X = orth(randn(n,p));
        [out] = test_OPP(X,n,p,d,1e-4,A,C,'rand_without_replacement');
        res_temp = [out.fval,out.feasi,out.itr,out.time,norm(out.X'*out.X-eye(p))];
        res = res+res_temp;
    end
    res = res/repnum;
    writematrix([d,res],tablename,'Sheet',sprintf('n=%d,p=%d',n,p), ...
            'Range',sprintf('A%d:Z%d',row_sheet,row_sheet));
    row_sheet = row_sheet+1;
end

%% OPP: n=1500,p=60
n = 1500;
p = 60;

load('OPP_A_1500_1500.mat','A');
load('OPP_C_1500_60.mat','C');
% load('OPP_X_1500_60.mat','X');

row_sheet = 2;

writecell(title,tablename,'Sheet',sprintf('n=%d,p=%d',n,p),'Range','A1:Z1');

repnum = 10; % 重复试验次数

for d = [6,12,18,24,30,36,42,48,54,60]
    res = zeros(1,5);
    for repeat = 1:repnum
        X = orth(randn(n,p));
        [out] = test_OPP(X,n,p,d,1e-4,A,C,'rand_without_replacement');
        res_temp = [out.fval,out.feasi,out.itr,out.time,norm(out.X'*out.X-eye(p))];
        res = res+res_temp;
    end
    res = res/repnum;
    writematrix([d,res],tablename,'Sheet',sprintf('n=%d,p=%d',n,p), ...
            'Range',sprintf('A%d:Z%d',row_sheet,row_sheet));
    row_sheet = row_sheet+1;
end

%% OPP: n=2000,p=80
n = 2000;
p = 80;

load('OPP_A_2000_2000.mat','A');
load('OPP_C_2000_80.mat','C');
% load('OPP_X_2000_80.mat','X');

row_sheet = 2;

writecell(title,tablename,'Sheet',sprintf('n=%d,p=%d',n,p),'Range','A1:Z1');

repnum = 10; % 重复试验次数

for d = [8,16,24,32,40,48,56,64,72,80]
    res = zeros(1,5);
    for repeat = 1:repnum
        X = orth(randn(n,p));
        [out] = test_OPP(X,n,p,d,1e-4,A,C,'rand_without_replacement');
        res_temp = [out.fval,out.feasi,out.itr,out.time,norm(out.X'*out.X-eye(p))];
        res = res+res_temp;
    end
    res = res/repnum;
    writematrix([d,res],tablename,'Sheet',sprintf('n=%d,p=%d',n,p), ...
            'Range',sprintf('A%d:Z%d',row_sheet,row_sheet));
    row_sheet = row_sheet+1;
end
