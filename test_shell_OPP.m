foldername = 'all_result_final';
if ~isfolder(foldername)
    mkdir(foldername)
end

tablename = [foldername,'/result_OPP','.xls'];

title = {'d','fval','normG','itr','time/s','feasi:X^TX-I'};

%% n = 1000, p = 500
n = 1000;
p = 500;
A = randn(n)/sqrt(n);
B = randn(n,p);

row_sheet = 2;

writecell(title,tablename,'Sheet',sprintf('n=%d,p=%d',n,p),'Range','A1:Z1');

repnum = 10; % 重复试验次数

for d = 50:50:500
    res = zeros(1,5);
    for repeat = 1:repnum
        X = orth(randn(n,p));
        [out] = test_OPP(X,n,p,d,1e-4,A,B,'rand_without_replacement');
        res_temp = [out.fval,out.normG,out.itr,out.time,norm(out.X'*out.X-eye(p))];
        res = res+res_temp;
    end
    res = res/repnum;
    writematrix([d,res],tablename,'Sheet',sprintf('n=%d,p=%d',n,p), ...
            'Range',sprintf('A%d:Z%d',row_sheet,row_sheet));
    row_sheet = row_sheet+1;
end

%% n = 10000, p = 10
n = 10000;
p = 10;
A = randn(n)/sqrt(n);
B = randn(n,p);

row_sheet = 2;

writecell(title,tablename,'Sheet',sprintf('n=%d,p=%d',n,p),'Range','A1:Z1');

repnum = 10; % 重复试验次数

for d = 1:10
    res = zeros(1,5);
    for repeat = 1:repnum
        X = orth(randn(n,p));
        [out] = test_OPP(X,n,p,d,1e-4,A,B,'rand_without_replacement');
        res_temp = [out.fval,out.normG,out.itr,out.time,norm(out.X'*out.X-eye(p))];
        res = res+res_temp;
    end
    res = res/repnum;
    writematrix([d,res],tablename,'Sheet',sprintf('n=%d,p=%d',n,p), ...
            'Range',sprintf('A%d:Z%d',row_sheet,row_sheet));
    row_sheet = row_sheet+1;
end

%%
n = 1000;
p = 500;
A = randn(n)/sqrt(n);
B = randn(n,p);

x = 300:50:500;

for d = x
        X = orth(randn(n,p));
        [out] = test_OPP(X,n,p,d,1e-4,A,B,'rand_without_replacement');
        time1(round(d/50)) = out.time;
end

for d = x
        X = orth(randn(n,p));
        [out] = test_OPP(X,n,p,d,1e-4,A,B,'rand_without_replacement');
        time2(round(d/50)) = out.time;
end

for d = x
        X = orth(randn(n,p));
        [out] = test_OPP(X,n,p,d,1e-4,A,B,'rand_without_replacement');
        time3(round(d/50)) = out.time;
end

for d = x
        X = orth(randn(n,p));
        [out] = test_OPP(X,n,p,d,1e-4,A,B,'rand_without_replacement');
        time4(round(d/50)) = out.time;
end

for d = x
        X = orth(randn(n,p));
        [out] = test_OPP(X,n,p,d,1e-4,A,B,'rand_without_replacement');
        time5(round(d/50)) = out.time;
end

%%
time_average = (time1+time2+time3+time4+time5)/5;

%%
figure

p1 = plot(x,time1(6:10),'c-','LineWidth',1);
hold on;
p2 = plot(x,time2(6:10),'c-','LineWidth',1);
hold on;
p3 = plot(x,time3(6:10),'c-','LineWidth',1);
hold on;
p4 = plot(x,time4(6:10),'c-','LineWidth',1);
hold on;
p5 = plot(x,time5(6:10),'c-','LineWidth',1);
hold on;
p6 = plot(x,time_average(6:10),'b-o','LineWidth',2);

title('d vs time');
legend([p1,p6],{'sample','average'});

xlabel('d');
ylabel('time');

%%
n = 10000;
p = 10;
A = randn(n)/sqrt(n);
B = randn(n,p);

x = 6:10;

for d = x
        X = orth(randn(n,p));
        [out] = test_OPP(X,n,p,d,1e-4,A,B,'rand_without_replacement');
        time1(round(d-5)) = out.time;
end

for d = x
        X = orth(randn(n,p));
        [out] = test_OPP(X,n,p,d,1e-4,A,B,'rand_without_replacement');
        time2(round(d-5)) = out.time;
end

for d = x
        X = orth(randn(n,p));
        [out] = test_OPP(X,n,p,d,1e-4,A,B,'rand_without_replacement');
        time3(round(d-5)) = out.time;
end

for d = x
        X = orth(randn(n,p));
        [out] = test_OPP(X,n,p,d,1e-4,A,B,'rand_without_replacement');
        time4(round(d-5)) = out.time;
end

for d = x
        X = orth(randn(n,p));
        [out] = test_OPP(X,n,p,d,1e-4,A,B,'rand_without_replacement');
        time5(round(d-5)) = out.time;
end

time_average = (time1+time2+time3+time4+time5)/5;

%%
figure

p1 = plot(x,time1(1:5),'c-','LineWidth',1);
hold on;
p2 = plot(x,time2(1:5),'c-','LineWidth',1);
hold on;
p3 = plot(x,time3(1:5),'c-','LineWidth',1);
hold on;
p4 = plot(x,time4(1:5),'c-','LineWidth',1);
hold on;
p5 = plot(x,time5(1:5),'c-','LineWidth',1);
hold on;
p6 = plot(x,time_average(1:5),'b-o','LineWidth',2);

title('d vs time');
legend([p1,p6],{'sample','average'});

xlabel('d');
ylabel('time');