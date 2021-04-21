clear all
e=10;
k1 = 1; k2 = 1; m1 = 1; m2 = 0.5; c0 = 2;
An{1} = [0 0 1/m1 0; 0 0 0 1/m2; -(k1+k2)  k2 -c0/m1 0; k2 -k2 0 -c0/m2];
An{2} = [0 0 1/m1 0; 0 0 0 1/m2; -(k1+k2)  k2 -c0/m1 0; k2  k2 0 -c0/m2];
An{3} = [0 0 1/m1 0; 0 0 0 1/m2; -(k1+k2) -k2 -c0/m1 0; k2  k2 0 -c0/m2];
An{4} = [0 0 1/m1 0; 0 0 0 1/m2; -(k1+k2) -k2 -c0/m1 0; k2 -k2 0 -c0/m2];
B{1} = [0; 0; 1; 0];
C{1} = [0 1 0 0];
Bw{1} = [0; 0; 1; 0];
H = [0 0 0 2; 0 0 4 0; 0 0 0 -2; 0 0 -4 0];
E = [0 0 0 0; 0 0 0 0; 0 0 0 1; 0 0 1 0];
mumelhor=200;
N = size(An,2);
nx = size(An{1},1);
q = size(B{1},2);
Xmelhor = 15206;
alfa=0;
for i=2:N
    B{i}=B{1};
    C{i}=C{1};
    Bw{i}=Bw{1};
end
%% Declara��o das vari�veis da LMI

for n1 = 1:length(alfa)
    alpha=alfa(n1);
    a=1;
    
    while a>0
        
        P = sdpvar(nx,nx,'symmetric');
        X = sdpvar(1,1,'full');
        F1 = sdpvar(1,nx,'full');
        F2 = sdpvar(nx,nx,'full');
        F3 = sdpvar(nx,nx,'full');
        W = sdpvar(q,nx,'full');
        gama=sdpvar(1,1);
        %F4=alpha*F3';
        
        if a==1
            obj = trace(X);
        else
            X = Xf;
            obj = [];
        end
        
        
        %% LMIs programadas
        LMIs = [];
        for i = 1:N
        LMI1 = [C{i}*F1'+F1*C{i}'-X       C{i}*F2-F1;
                F2'*C{i}'-F1'             P-F2-F2'];
        
        J = [An{i}*F3'+B{i}*W+F3*An{i}'+W'*B{i}'          P-F3+An{i}*alpha*F3'+alpha*B{i}*W                   Bw{i};
            P-F3'+alpha*F3*An{i}'+alpha*W'*B{i}'               -alpha*F3'-alpha*F3                                      zeros(nx,q);
            Bw{i}'                                        zeros(q,nx)                                  -diag(ones(q,1))];
        M=[E*F3' E*alpha*F3' zeros(nx,q)];
        N1=[H' zeros(nx,nx) zeros(nx,q)];
        LMI2 = [J M' N1';
            M -gama*diag(ones(nx,1)) zeros(nx,nx);
            N1 zeros(nx,nx) -gama*diag(ones(nx,1))];
        LMIs = LMIs+(LMI1<=0)+(LMI2<=0)+(P>=0);
        end
        % solver
        options = sdpsettings('savesolveroutput',1,'verbose',1,'warning',1,'solver','sedumi','showprogress',0);
        %solucao=
        % Tentativa de solu��o
        tempo = cputime;
        solvesdp(LMIs,obj,options);
        tempo = cputime-tempo;
        
        %% Verifica a necessidade de outro loop
        teste = min(checkset(LMIs));
        if teste <0 && a==1
            Xf = 1.001*double(trace(X));
            a=a+1;
        else
            break
        end
    end
    
    if teste > -10^(-7)
        Resp = 1; % flag de factibilidade
        if(double(trace(X))<Xmelhor)
            Xmelhor=double(trace(X));
            Pmelhor=double(P);
            Wmelhor=double(W);
            F1melhor=double(F1);
            F2melhor=double(F2);
            F3melhor=double(F3);
            gamamelhor=double(gama);
            alphamelhor=alpha;
            Kmelhor =  Wmelhor/F3melhor';
            h2melhor = sqrt(double(trace(X)));
        end
        
    else
        Resp = -1;
    end 
end