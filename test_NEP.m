function [out] = test_NEP(X,n,p,d,stop_tol,sample_method)

beta = 1;

if(~exist('sample_method','var'))
    % 样本选取方式，默认为无放回随机抽样
    sample_method = 'rand_without_replacement';
end

flag1 = 0;%记录矩阵重新正交化的次数
flag2 = 0;

L = gallery('tridiag', n, -1, 2, -1);
[Ll,Lu] = lu(L);

% X = orth(randn(n,p));

% 终止条件
inner_maxit = 5;
maxit = 30000;

% sqrt_np = sqrt(n*p);
sqrt_nd = sqrt(n*d);
% tol = 1e-6; % 外层迭代终止条件（feasi)
% tol = tol*sqrt_np;
% stop_tol = 1e-6;  % 子问题的终止条件
stop_tol = stop_tol*sqrt_nd;

tiny = 1e-12;

alpha = 0.2;
gamma = 2;
rho = 1e-4;

% record = [];
% feasible_value = [];
record = zeros(maxit,1);
feasible_value = zeros(maxit,1);
step_record = zeros(maxit,1);


break_judge = 0;
time = 0;

idx = 1:p;  % idx表示可以选择的列的指标，初始值是1:p

tau = 1;

tic

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

    [F,G] = fun(X); % 先求出原始问题函数值及其梯度
    Gb = G(:,index_B); % Gb是子模型梯度
    Gb = Gb-Xn*(Xn'*Gb); % Gb是投影梯度
    % dtXb = Gb - Xb*(Gb'*Xb); % dtXb是流形梯度
    % nrmG = norm(dtXb, 'fro'); % nrmG是流形梯度的模，后面会有先搜索判断用到

    eye2d = eye(2*d); % 在SMW公式中需要用到这个2d维度的单位矩阵
    % U =  [Gb, Xb];
    % V = [Xb, -Gb];
    % VU = V'*U;
    % VXb = V'*Xb;

    % RL = Xb'*Gb;
    % XnGb = Xn'*Gb;
    % minusLL = -Gb'*Gb;
    % VU = [RL,eye(d);minusLL,-RL'];
    % VXb = [eye(d);-RL'];
    % nrmGsquare = abs(norm(Gb,'fro')^2-norm(XnGb,'fro')^2-sum(dot(RL,RL,1)));

    GXb = Gb'*Xb;

    U =  [Gb, Xb];    V = [Xb, -Gb];       VU = V'*U;
    VXb = V'*Xb;

    dtX = Gb - Xb*GXb; 
    nrmG  = norm(dtX, 'fro');

    nrmGsquare = nrmG^2;



    % Q = 1; % 非单调线搜索方法
    % Cval = F;
    % tau = 1e-1; % 初始步长

    rhoXn = sum(Xn.^2,2);
    % FP = funsub(Xb);

    for inner_itr = 1:inner_maxit % 内部迭代过程
       
        XbP = Xb;
        FP = funsub(Xb);
        % dtXbP = dtXb;

        nls = 1; % 用于终止回溯法的指标
        % deriv = rho*nrmG^2; % 线搜索后面的值
        deriv = rho*nrmGsquare;

        while 1 %开始迭代
            % [aa, ~] = linsolve(eye2d + (0.5*tau)*VU, VXb);
            % d=1时报错"linsolve 不支持以稀疏矩阵作为输入项"
            % 反斜杠 \ 对应mldivide函数
            % mldivide 会对检查输入矩阵是否具有特殊属性，比linsolve慢一点

            % aa = (eye2d + (0.5*tau)*VU)\VXb;
            aa = linsolve(eye2d + (0.5*tau)*VU,VXb);

            Xb = XbP - [Gb,XbP]*(tau*aa);
            % if norm(Xb'*Xb - eye(d),'fro') > tiny; Xb = myQR(Xb,d); flag1 = flag1+1; end
            % X(:,index_B) = Xb;

            % [F,G] = fun(X);
            F = funsub(Xb);
            if F <= FP - tau*deriv || nls >= 10
                break
            end
            tau = alpha*tau;
            nls = nls+1;
        end

        step_record(itr) = tau;

        X(:,index_B) = Xb;
        [F,G] = fun(X);
        Gb = G(:,index_B);
        Gb = Gb-Xn*(Xn'*Gb);
        % dtXb = Gb - Xb*(Gb'*Xb);
        % nrmG  = norm(dtXb, 'fro');
        % RL = Xb'*Gb;
        % XnGb = Xn'*Gb;
        % minusLL = -Gb'*Gb;
        % VU = [RL,eye(d);minusLL,-RL'];
        % VXb = [eye(d);-RL'];
        % nrmGsquare = abs(norm(Gb,'fro')^2-norm(XnGb,'fro')^2-sum(dot(RL,RL,1)));

        GXb = Gb'*Xb;

        U =  [Gb, Xb];    V = [Xb, -Gb];       VU = V'*U;
        VXb = V'*Xb;

        dtX = Gb - Xb*GXb; 
        nrmG  = norm(dtX, 'fro');

        nrmGsquare = nrmG^2;


        % S = Xb - XbP;
        % nrmdX = norm(S,'fro');
        % diffF = abs(F-FP)/(abs(F)+1e-2);

        %          if nrmG<stop_tol ||(nrmdX < xtol && diffF < ftol)
        %             break_judge = 1;
        %             break
        %          end

        if nrmGsquare<stop_tol^2
            % 子问题 nrmG 足够小则停止整个算法
            out.msg = 'converge';
            break_judge = 1;
            break
        end

        % if nrmdX < xtol && diffF < ftol
        %     out.msg = 'nrmdX < xtol or diffF < ftol';
        %     break_judge = 1;
        %     break
        % end

        % U =  [Gb, Xb];
        % V = [Xb, -Gb];
        % VU = V'*U;
        % VXb = V'*Xb;

        % Qp = Q;
        % Q = gamma*Qp + 1;
        % Cval = (gamma*Qp*Cval + F)/Q;

        % Y = dtXb - dtXbP;
        % SY = abs(iprod(S,Y));
        % if mod(inner_itr,2)==0
        %     tau = (norm(S,'fro')^2)/SY;
        % else
        %     tau  = SY/(norm(Y,'fro')^2);
        % end
        % % tau = max(min(tau, 1e20), 1e-20);
        % tau = max(min(tau, 1e5), 1e-5);       %上述代码为使用BB步长
        
        tau = tau*gamma;
        % tau = 5e-2;

    end

    % record = [record,F];
    % feasible_value = [feasible_value,nrmG];
    % 重新调整数组大小需要MATLAB花费额外的时间来寻找更大的连续内存块，然后将数组移入新地址。
    % 可通过预分配数组所需的最大空间量来缩短代码的执行时间。
    record(itr) = F;
    % feasible_value(itr) = nrmG;
    feasible_value(itr) = nrmGsquare;

    add = toc;
    time = time+add;

    if norm(X'*X - eye(p),'fro') > tiny
        X = myQR(X,p);
        flag2 = flag2+1;
%         [F,G] = fun(X);
        break_judge = 0;
    end

    if break_judge
        break
    end

    tic

    % grad = G-X*G'*X;

    % diff_x = norm(X-XP,'fro')/sqrt(n*p);
    % diff_f = abs(F-FP)/(abs(F)+1e-2);
    %
    % if diff_x < xtol && diff_f < ftol
    %     out.msg = 'diff_x < xtol and diff_f < ftol';
    %     break
    % end
    %
    % if feasible_value(itr) < tol
    %     out.msg = 'converge';
    %     break
    % end
end

if itr == maxit
    out.msg = 'exceed max iteration';
end

out.record = record(1:itr);
out.feasible_value = feasible_value(1:itr);
out.time = time;
out.fval = F;
% out.feasi = nrmG;
out.feasi = (nrmGsquare^0.5)/sqrt_nd;
out.itr = itr;
out.X = X;
out.flag1 = flag1;
out.flag2 = flag2;
out.step_record = step_record(1:itr);

    function [f,g] = fun(X)
        LX = L*X;
        rhoX = sum(X.^2, 2);
        tempa = Lu\(Ll\rhoX);
        tempa = beta*tempa;
        f = 0.5*sum(sum(X.*(LX))) + 1/4*(rhoX'*tempa);
        g = LX + bsxfun(@times,tempa,X);
    end

    function [Q, RR] = myQR(XX,k)
        [Q, RR] = qr(XX, 0);
        diagRR = sign(diag(RR)); ndr = diagRR < 0;
        if nnz(ndr) > 0
            Q = Q*spdiags(diagRR,0,k,k);
        end
    end

    function a = iprod(x,y)
        a = real(sum(sum(conj(x).*y)));
    end

    function [f] = funsub(Xb)
        LXb = L*Xb;
        rhoXb = sum(Xb.^2,2);
        % rhoXn = sum(Xn.^2,2); 在全局计算
        tempa = Lu\(Ll\rhoXb);
        tempa = beta*tempa;
        f = 0.5*sum(sum(Xb.*(LXb))) + 1/4*((rhoXb+2*rhoXn)'*tempa);
    end

end
