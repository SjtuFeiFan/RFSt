function [out] = test_LEP(X,n,p,d,stop_tol,A,sample_method)
% varargin只存放矩阵A

if(~exist('sample_method','var'))
    % 样本选取方式，默认为无放回随机抽样
    sample_method = 'rand_with_replacement';
end

flag1 = 0;%记录矩阵重新正交化的次数
flag2 = 0;

% 终止条件
inner_maxit = 5;
maxit = 1000;

sqrt_nd = sqrt(n*d);
stop_tol = stop_tol*sqrt_nd;

tiny = 1e-12;

alpha = 0.25;
rho = 1e-4;
gamma = 2;

record = zeros(maxit,1);
feasible_value = record;
tau_record = record;


break_judge = 0;
time = 0;

idx = 1:p;  % idx表示可以选择的列的指标，初始值是1:p

tic

tau = 1; % 初始步长

for itr=1:maxit
    switch sample_method
        case 'rand_with_replacement'
            % 有放回随机抽样
            index_B = sort(randperm(p,d));
        case 'rand_without_replacement'
            % 无放回抽样，使所有列都得到更新
            if length(idx) > d
                % 如果idx中指标个数大于d
                index_B = sort(randsample(idx,d)); % 从idx集合中随机选取d列
                idx = setdiff(idx,index_B); % 从idx集合中去除所选的列
            else
                % 如果idx中指标个数小于或等于d
                % 从setdiff(1:p,idx)中再选d-length(idx)个指标，与idx中所有元素组成index_B
                index_B = sort([idx,randsample(setdiff(1:p,idx),d-length(idx))]);
                idx = 1:p; % 将idx集合重置为1:p
            end
    end
    index_N = setdiff(1:p,index_B);

    Xb = X(:,index_B);
    Xn = X(:,index_N);

    Gb = -(A*Xb); % Gb是子模型梯度
    Fb = 0.5*sum(dot(Gb,Xb,1));
    Gb = Gb-Xn*(Xn'*Gb); % Gb是投影梯度
    GXb = Gb'*Xb;
    eye2d = eye(2*d); % 在SMW公式中需要用到这个2d维度的单位矩阵

    U =  [Gb, Xb];    V = [Xb, -Gb];       VU = V'*U;
    VXb = V'*Xb;

    dtX = Gb - Xb*GXb; 
    nrmG  = norm(dtX, 'fro');


    % RL = Xb'*Gb;
    % XnGb = Xn'*Gb;
    % minusLL = -Gb'*Gb;
    % VU = [RL,eye(d);minusLL,-RL'];
    % VXb = [eye(d);-RL'];
    % nrmGsquare = abs(norm(Gb,'fro')^2-norm(XnGb,'fro')^2-sum(dot(RL,RL,1)));
    nrmGsquare = nrmG^2;

    % tau = 3e-4; % 初始步长

    for inner_itr = 1:inner_maxit % 内部迭代过程
       
        XbP = Xb;
        % FP = funsub(Xb);
        % dtXbP = dtXb;
        FPsub = Fb;

        nls = 1; % 用于终止回溯法的指标
        deriv = rho*nrmGsquare;

        while 1 
            % aa = (eye2d + (0.5*tau)*VU)\VXb;
            aa = linsolve(eye2d + (0.5*tau)*VU,VXb);
            Xb = XbP - [Gb,XbP]*(tau*aa);
            % if norm(Xb'*Xb - eye(d),'fro') > tiny; Xb = myQR(Xb,d); flag1 = flag1+1; end
            Gbtest = -(A*Xb);
            F = 0.5*sum(dot(Gbtest,Xb,1));
            if F <= FPsub - tau*deriv || nls >= 5
                break
            end
            tau = alpha*tau;
            nls = nls+1;
        end
        tau_record(itr) = tau;

        X(:,index_B) = Xb;

        Gb = Gbtest;
        Fb = F;
        Gb = Gb-Xn*(Xn'*Gb);
        GXb = Gb'*Xb;

        U =  [Gb, Xb];    V = [Xb, -Gb];       VU = V'*U;
        VXb = V'*Xb;

        dtX = Gb - Xb*GXb; 
        nrmG  = norm(dtX, 'fro');

        % RL = Xb'*Gb;
        % XnGb = Xn'*Gb;
        % minusLL = -Gb'*Gb;
        % VU = [RL,eye(d);minusLL,-RL'];
        % VXb = [eye(d);-RL'];
        % nrmGsquare = abs(norm(Gb,'fro')^2-norm(XnGb,'fro')^2-sum(dot(RL,RL,1)));
        
        nrmGsquare = nrmG^2;
        
        if nrmGsquare<stop_tol^2
            out.msg = 'converge';
            break_judge = 1;
            break
        end

        
        
        tau = tau*gamma;

    end
    feasible_value(itr) = nrmGsquare;

    add = toc;
    time = time+add;

    if norm(X'*X - eye(p),'fro') > tiny
        X = myQR(X,p);
        flag2 = flag2+1;
        break_judge = 0;
    end

    if break_judge
        break
    end

    tic


end

if itr == maxit
    out.msg = 'exceed max iteration';
end

out.record = record(1:itr);
out.feasible_value = feasible_value(1:itr);
out.time = time;
out.fval = -0.5*sum(dot(A*X,X,1));
out.normG = (nrmGsquare^0.5)/sqrt_nd;
out.itr = itr;
out.X = X;
out.flag1 = flag1;
out.flag2 = flag2;
out.tau_record = tau_record(1:itr);

    function [f,g] = fun(X,A)
        g = -(A*X);
        f = 0.5*sum(dot(g,X,1));
    end

    function [Q, RR] = myQR(XX,k)
        [Q, RR] = qr(XX, 0);
        diagRR = sign(diag(RR)); ndr = diagRR < 0;
        if nnz(ndr) > 0
            Q = Q*spdiags(diagRR,0,k,k);
        end
    end


end