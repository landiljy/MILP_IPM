function [ ] = MILP_IPM(  )
%MILP_IPM实现不计损耗的DED-VPE
%   此处显示详细说明
%5机组数据
%Dt:24小时负荷需求
%ai:机组系数alpha
%bi:机组系数beta
%ci:机组系数gama
%ei:机组系数
%fi:机组系数
%Pimin:最小实际输出功率
%Pimax:最大实际输出功率
%UR:上坡功率限制
%DR:下坡功率限制
%N:机组总数
%T:时段总数
%B:损耗公式系数
%PL;传输损耗
%Li:（23）式中的Li
%kli:（23）式中的kl,i
%bli:（23）式中的bl,i
%c:（23）式中的c（）
%ali:（23）式中的al,i
%zlit:(22)式中的zl,i,t
%plit:
%M:等分段数
%y:
%z:
%SRit:机组i在t时间的热备用
%P_it:

Dt = [410 435 475 530 558 608 626 654 690 704 720 740 704 690 654 580 558 608 654 704 680 605 527 463]';
ai = [25 60 100 120 40]';
bi = [2.0 1.8 2.1 2.0 1.8]';
ci = [0.0080 0.0030 0.0012 0.0010 0.0015]';
ei = [100 140 160 180 200]';
fi = [0.042 0.040 0.038 0.037 0.035]';
Pimin = [10 20 30 40 50]';
Pimax = [75 125 175 250 300]';
UR = [30 30 40 50 50]';
DR = [30 30 40 50 50]';
B = [0.000049 0.000014 0.000015 0.000015 0.000020;
    0.000014 0.000045 0.000016 0.000020 0.000018;
    0.000015 0.000016 0.000039 0.000010 0.000012;
    0.000015 0.000020 0.000010 0.000040 0.000014;
    0.000020 0.000018 0.000012 0.000014 0.000035]; 

N= 5;%机组数量
T = 24;%时间数量 
Pit = sdpvar(N,T);%机组i在t时间的输出功率
PL = sdpvar(T,1,'full');%传输损耗
s = sdpvar(N,T);%辅助变量
u = sdpvar(N,T);%松弛变量
v = sdpvar(N,T);%松弛变量
SRit = sdpvar(N,T);%机组i在t时间的热备用
y = 1;%(7)式中的热备用系数
z = 0.05;%(8)式中的热备用系数
M = 4;%分隔段数
kli = [];%(23)式中的kl,i
Li = [];%(23)式中的Li
ali = [];
bli = [];

%计算出每个机组的分割点数
for i = 1:N
    Li(i) = ceil(M*(fi(i)*(Pimax(i)-Pimin(i))/pi))
end

%计算出每个机组的分割点的x坐标
start = 1;
for i = 1:N
    add = (Pimax(i) - Pimin(i))/Li(i);
    ali(1,start) = Pimin(i);
    for l = start:Li(i)+start
        ali(1,l+1) = ali(1,l) + add;
        if i == size(Li,2) & l == Li(i) + start - 1
            break;
        end
    end
    start = start + Li(i) + 1;
end
ali
%计算出每个机组每个分割段中两两相邻线性点连线的斜率和(23)式中的bl,i
start = 1;
ed = Li(1);
for i = 1:N
    add = (Pimax(i) - Pimin(i))/Li(i)/M;
    for l = start:ed
        for m = 1: M
            kli(1,(l - 1)*M + m) = ((ai(i) + bi(i) * (ali(1,l) + m*add) + ...
            ci(i) * (ali(1,l) + m*add) * (ali(1,l) + m*add) ...
            + ei(i) * abs(sin(fi(i) * (ali(1,l) + m*add - Pimin(i)))))...
            - (ai(i) + bi(i) * (ali(1,l) + (m - 1)*add) + ...
            ci(i) * (ali(1,l) + (m - 1)*add) * (ali(1,l) + (m - 1) *add)...
            + ei(i) * abs(sin(fi(i) * (ali(1,l) + (m - 1)*add - Pimin(i)))))) ...
            / (ali(1,l) + m*add - ali(1,l) - (m - 1)*add);
        
            bli(1,(l - 1)*M + m) = (ai(i) + bi(i) * (ali(1,l) + (m - 1)*add) + ...
            ci(i) * (ali(1,l) + (m - 1)*add) * (ali(1,l) + (m - 1) *add)...
            + ei(i) * abs(sin(fi(i) * (ali(1,l) + (m - 1)*add - Pimin(i)))))...
            - kli(1,(l - 1)*M + m)*(ali(1,l) + (m - 1)*add);
        end
    end
    start = start + Li(i) + 1;
    if i ~= N
        ed = start + Li(i + 1) - 1;
    end
end
kli
bli
%计算（22）式中的约束条件
constraint = [];%约束条件
zlit = binvar(T,(sum(Li)+size(Li,2)-1)*M);%行为时间，列为每个机组的某个分割段的
plit = sdpvar(T,(sum(Li)+size(Li,2)-1)*M);
start = 1;
ed = Li(1);
for t = 1:T
    for i = 1:N
        constraint = [constraint,sum(zlit(t,(start - 1)*M + 1 : ed*M)) == 1];
        constraint = [constraint,sum(plit(t,(start - 1)*M + 1 : ed*M)) == Pit(i,t)];   
        if i == N
            start = 1;
            ed = Li(1);
            break;
        end
        start = start + Li(i) + 1;
        if i ~=N
            ed = start + Li(i+1) - 1;
        end
        
    end
end
start = 1;
ed = Li(1);
for t = 1:T
    for i = 1:N
        add = (Pimax(i) - Pimin(i))/Li(i)/M;
        for l = start:ed
            for m = 1:M
                constraint = [constraint,(ali(1,l) + add*(m-1))*zlit(t,(l - 1)*M + m) <= plit(t,(l - 1)*M + m) <= (ali(1,l) + add*m)*zlit(t,(l - 1)*M + m)];
            end
        end
        if i == N
            start = 1;
            ed = Li(1);
            break;
        end
        start = start + Li(i) + 1;
        if i ~= N
            ed = start + Li(i+1) - 1;
        end
        
    end
end

constraint = [constraint,sum(Pit)' - PL == Dt];%(3)式
for t = 2:T
    constraint = [constraint,-DR <= Pit(:,t)-Pit(:,t-1) <= UR];%爬坡约束(6)式
end

for t = 1:T
    constraint = [constraint,Pimin <= Pit(:,t) <= Pimax];%发电限制(5)式
    constraint = [constraint,SRit(:,t) <= Pimax - Pit(:,t)];%热备用约束（7）式
    constraint = [constraint,SRit(:,t) <= y * UR];%热备用约束(7)式
    constraint = [constraint,sum(SRit(:,t)) >= z*Dt(t)];%热备用约束（8）式
end
   
constraint = [constraint,Pit >= 0,s >= 0,PL >= 0,u >= 0,v >= 0,SRit >= 0,plit >= 0];

%计算（22）式中的目标函数
fun = 0;
start = 1;
ed = Li(1);
for t = 1:T
    for i = 1:N
        fun = fun + kli(1,(start - 1)*M + 1 : ed*M)*plit(t,(start - 1)*M + 1 : ed*M)' + bli(1,(start - 1)*M + 1 : ed*M)*zlit(t,(start - 1)*M + 1 : ed*M)';    
        if i == N
            start = 1;
            ed = Li(1);
            break;
        end
        start = start + Li(i) + 1;
        ed = start + Li(i+1) - 1;
    end
end
options = sdpsettings('verbose',1 ,'solver' ,'cplex');
options.cplex.exportmodel='model.lp';%输出数学模型
options.cplex.timelimit=50;%迭代时间限制，单位为秒
sol = optimize(constraint,fun,options);
disp(sol);
disp('fun=');
disp(value(fun));
disp(value(Pit));
disp(value(PL'));


P_it = sdpvar(N,T);%机组i在t时间的输出功率
PL = sdpvar(T,1,'full');%传输损耗
s = sdpvar(N,T);%辅助变量
u = sdpvar(N,T);%松弛变量
v = sdpvar(N,T);%松弛变量
SRit = sdpvar(N,T);%机组i在t时间的热备用

F = [];%约束条件
F = [F,sum(P_it)' - PL == Dt];%功率平衡方程
for t = 2:T
    F = [F,-DR <= P_it(:,t)-P_it(:,t-1) <= UR];%爬坡约束
end

for t = 1:T
    F = [F,Pimin <= P_it(:,t) <= Pimax];%发电限制
    F = [F,P_it(:,t)'*B*P_it(:,t) - PL(t) == 0];%损耗约束
    F = [F,sin(fi.*(P_it(:,t)-Pimin)) + u(:,t) - v(:,t) == 0];%等式(19)
    F = [F,s(:,t) == u(:,t) + v(:,t)];%等式(18) 
%     constraint = [constraint,SRit(:,t) <= Pimax - Pit(:,t)];%热备用约束（7）式
%     constraint = [constraint,SRit(:,t) <= y * UR];%热备用约束(7)式
%     constraint = [constraint,sum(SRit(:,t)) >= z*Dt(t)];%热备用约束（8）式
end
F = [F,P_it >= 0,s >= 0,PL >= 0,u >= 0,v >= 0];
f = 0;

for t = 1:T
    f = f + ai'*ones(N,1) + bi'* P_it(:,t) +  ci' * (P_it(:,t).*P_it(:,t)) + ei' * s(:,t);%目标函数
end
assign(P_it,double(Pit));%将cplex求出的Pit赋值给P_it；
options = sdpsettings('verbose',1 ,'solver' ,'ipopt','usex0',1);%usex0设为1，即为设置初值。
sol = optimize(F,f,options);
disp(sol);
disp('f=');
disp(value(f));
disp(value(P_it));
disp(value(PL'));
end

