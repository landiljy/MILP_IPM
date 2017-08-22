function [ ] = MILP_IPM(  )
%MILP_IPMʵ�ֲ�����ĵ�DED-VPE
%   �˴���ʾ��ϸ˵��
%5��������
%Dt:24Сʱ��������
%ai:����ϵ��alpha
%bi:����ϵ��beta
%ci:����ϵ��gama
%ei:����ϵ��
%fi:����ϵ��
%Pimin:��Сʵ���������
%Pimax:���ʵ���������
%UR:���¹�������
%DR:���¹�������
%N:��������
%T:ʱ������
%B:��Ĺ�ʽϵ��
%PL;�������
%Li:��23��ʽ�е�Li
%kli:��23��ʽ�е�kl,i
%bli:��23��ʽ�е�bl,i
%c:��23��ʽ�е�c����
%ali:��23��ʽ�е�al,i
%zlit:(22)ʽ�е�zl,i,t
%plit:
%M:�ȷֶ���
%y:
%z:
%SRit:����i��tʱ����ȱ���
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

N= 5;%��������
T = 24;%ʱ������ 
Pit = sdpvar(N,T);%����i��tʱ����������
PL = sdpvar(T,1,'full');%�������
s = sdpvar(N,T);%��������
u = sdpvar(N,T);%�ɳڱ���
v = sdpvar(N,T);%�ɳڱ���
SRit = sdpvar(N,T);%����i��tʱ����ȱ���
y = 1;%(7)ʽ�е��ȱ���ϵ��
z = 0.05;%(8)ʽ�е��ȱ���ϵ��
M = 4;%�ָ�����
kli = [];%(23)ʽ�е�kl,i
Li = [];%(23)ʽ�е�Li
ali = [];
bli = [];

%�����ÿ������ķָ����
for i = 1:N
    Li(i) = ceil(M*(fi(i)*(Pimax(i)-Pimin(i))/pi))
end

%�����ÿ������ķָ���x����
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
%�����ÿ������ÿ���ָ���������������Ե����ߵ�б�ʺ�(23)ʽ�е�bl,i
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
%���㣨22��ʽ�е�Լ������
constraint = [];%Լ������
zlit = binvar(T,(sum(Li)+size(Li,2)-1)*M);%��Ϊʱ�䣬��Ϊÿ�������ĳ���ָ�ε�
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

constraint = [constraint,sum(Pit)' - PL == Dt];%(3)ʽ
for t = 2:T
    constraint = [constraint,-DR <= Pit(:,t)-Pit(:,t-1) <= UR];%����Լ��(6)ʽ
end

for t = 1:T
    constraint = [constraint,Pimin <= Pit(:,t) <= Pimax];%��������(5)ʽ
    constraint = [constraint,SRit(:,t) <= Pimax - Pit(:,t)];%�ȱ���Լ����7��ʽ
    constraint = [constraint,SRit(:,t) <= y * UR];%�ȱ���Լ��(7)ʽ
    constraint = [constraint,sum(SRit(:,t)) >= z*Dt(t)];%�ȱ���Լ����8��ʽ
end
   
constraint = [constraint,Pit >= 0,s >= 0,PL >= 0,u >= 0,v >= 0,SRit >= 0,plit >= 0];

%���㣨22��ʽ�е�Ŀ�꺯��
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
options.cplex.exportmodel='model.lp';%�����ѧģ��
options.cplex.timelimit=50;%����ʱ�����ƣ���λΪ��
sol = optimize(constraint,fun,options);
disp(sol);
disp('fun=');
disp(value(fun));
disp(value(Pit));
disp(value(PL'));


P_it = sdpvar(N,T);%����i��tʱ����������
PL = sdpvar(T,1,'full');%�������
s = sdpvar(N,T);%��������
u = sdpvar(N,T);%�ɳڱ���
v = sdpvar(N,T);%�ɳڱ���
SRit = sdpvar(N,T);%����i��tʱ����ȱ���

F = [];%Լ������
F = [F,sum(P_it)' - PL == Dt];%����ƽ�ⷽ��
for t = 2:T
    F = [F,-DR <= P_it(:,t)-P_it(:,t-1) <= UR];%����Լ��
end

for t = 1:T
    F = [F,Pimin <= P_it(:,t) <= Pimax];%��������
    F = [F,P_it(:,t)'*B*P_it(:,t) - PL(t) == 0];%���Լ��
    F = [F,sin(fi.*(P_it(:,t)-Pimin)) + u(:,t) - v(:,t) == 0];%��ʽ(19)
    F = [F,s(:,t) == u(:,t) + v(:,t)];%��ʽ(18) 
%     constraint = [constraint,SRit(:,t) <= Pimax - Pit(:,t)];%�ȱ���Լ����7��ʽ
%     constraint = [constraint,SRit(:,t) <= y * UR];%�ȱ���Լ��(7)ʽ
%     constraint = [constraint,sum(SRit(:,t)) >= z*Dt(t)];%�ȱ���Լ����8��ʽ
end
F = [F,P_it >= 0,s >= 0,PL >= 0,u >= 0,v >= 0];
f = 0;

for t = 1:T
    f = f + ai'*ones(N,1) + bi'* P_it(:,t) +  ci' * (P_it(:,t).*P_it(:,t)) + ei' * s(:,t);%Ŀ�꺯��
end
assign(P_it,double(Pit));%��cplex�����Pit��ֵ��P_it��
options = sdpsettings('verbose',1 ,'solver' ,'ipopt','usex0',1);%usex0��Ϊ1����Ϊ���ó�ֵ��
sol = optimize(F,f,options);
disp(sol);
disp('f=');
disp(value(f));
disp(value(P_it));
disp(value(PL'));
end

