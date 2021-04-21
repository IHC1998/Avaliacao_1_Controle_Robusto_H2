clear all
An{1} = [-0.02 0.16-0.1 -0.6-0.1;-0.45+0.2 -0.05-0.1 0.59+0.3; -0.19-0.3 0.07-0.4 -0.57];
An{2} = [-0.02 0.16-0.1 -0.6+0.1;-0.45+0.2 -0.05+0.1 0.59+0.3; -0.19+0.3 0.07-0.4 -0.57];
An{3} = [-0.02 0.16+0.1 -0.6-0.1;-0.45-0.2 -0.05-0.1 0.59-0.3; -0.19-0.3 0.07+0.4 -0.57];
An{4} = [-0.02 0.16+0.1 -0.6+0.1;-0.45-0.2 -0.05+0.1 0.59-0.3; -0.19+0.3 0.07+0.4 -0.57];
B{1} = [1; 0; 0]; B{2}=B{1}; B{3}=B{1}; B{4}=B{1};
Xmelhor = 12598;
C{1} = [0 0 1]; C{2}=C{1}; C{3}=C{1}; C{4}=C{1};
Bw = B;
H = [0 1 0; -2 3 -2; -3 0 4];
E = [1 0 0; 0 1 1; 0 1 -1];
N = size(An,2);
nx = size(An{1},1);
q = size(B{1},2);
%% Declara��o das vari�veis da LMI

for alpha = 0
    for i=1:N
        P{i} = sdpvar(nx,nx,'symmetric');
    end
    X = sdpvar(1,1,'full');
    F1 = sdpvar(nx,nx,'full');
    F2 = sdpvar(1,nx,'full');
    F3 = sdpvar(nx,nx,'full');
    F4 = sdpvar(nx,nx,'full');
    W = sdpvar(q,nx,'full');
    gama=sdpvar(1,1);
    F5=alpha*F4';
    %% LMIs programadas
    LMIs = [];
    for i= 1:N
        for l =1:N
            LMI1 = [P{i}                                      F1*C{i}'                      P{i}-F1;
                C{i}*F1'                             X+F2*C{i}'+C{i}*F2'                -F2+C{i}*F3';
                P{i}-F1'                              -F2'+F3*C{i}'                   -F3-F3'];
            
            J = [P{i}-An{i}*F4'-B{i}*W-F4*An{i}'-W'*B{i}'          F4-An{i}*F5-alpha*B{i}*W                   Bw{i};
                F4'-F5'*An{i}'-alpha*W'*B{i}'                 -P{l}+F5+F5'                                    zeros(nx,q);
                Bw{i}'                                          zeros(q,nx)                                diag(ones(q,1))];
            
            M=[E*F4' E*F5 zeros(nx,q)];
            N1=[H' zeros(nx,nx) zeros(nx,q)];
            LMI2 = [-J M' N1';
                M -gama*diag(ones(nx,1)) zeros(nx,nx);
                N1 zeros(nx,nx) -gama*diag(ones(nx,1))];
            LMIs = LMIs+(LMI1>=0)+(LMI2<=0)+(P{i}>=0);
        end
    end
    % solver
    options = sdpsettings('savesolveroutput',1,'verbose',1,'warning',1,'solver','sedumi','showprogress',0);
    %solucao=
    % Tentativa de solu��o
    tempo = cputime;
    solvesdp(LMIs,trace(X),options);
    tempo = cputime-tempo;
    
    %% Verifica a necessidade de outro loop
    teste = min(checkset(LMIs));
    
    if teste > -10^-6
        Resp = 1; % flag de factibilidade
        if(double(trace(X))<Xmelhor)
            Xmelhor=double(trace(X));
            for i=1:N
                Pmelhor{i}=double(P{i});
            end
            Wmelhor=double(W);
            F1melhor=double(F1);
            F2melhor=double(F2);
            F3melhor=double(F3);
            F4melhor=double(F4);
            F5melhor=double(F5);
            gamamelhor=double(gama);
            alphamelhor=alpha;
            Kmelhor =  Wmelhor/F4melhor';
            H2=sqrt(double(trace(X)));
        end
        
    else
        Resp = -1;
    end
    
    
end