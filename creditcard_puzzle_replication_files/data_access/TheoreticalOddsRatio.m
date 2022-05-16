clc
clear

%% Transition Matrices

% probabilities
gain = 0.0607;
lose = 0.0263;
uast = 0.07;
chi_u = 4;

lose_u = chi_u*lose;
lose_w = (lose - lose_u * uast) / (1-uast);

% transition matrices
A = zeros(2);
A(1,1) = (1-lose);
A(2,1) = lose;
A(2,2) = (1-gain);
A(1,2) = gain;

A_w = A;
A_w(1,1) = (1-lose_w);
A_w(2,1) = lose_w;

A_u = A;
A_u(1,1) = (1-lose_u);
A_u(2,1) = lose_u;

%% Unemployment Paths

umat = zeros(4,16);

umat(:,1) = [0, 0, 0, 0];

umat(:,2) = [1, 0, 0, 0];
umat(:,3) = [0, 1, 0, 0];
umat(:,4) = [0, 0, 1, 0];
umat(:,5) = [0, 0, 0, 1];

umat(:,6) = [1, 1, 0, 0];
umat(:,7) = [1, 0, 1, 0];
umat(:,8) = [1, 0, 0, 1];
umat(:,9) = [0, 1, 1, 0];
umat(:,10) = [0, 1, 0, 1];
umat(:,11) = [0, 0, 1, 1];

umat(:,12) = [0, 1, 1, 1];
umat(:,13) = [1, 0, 1, 1];
umat(:,14) = [1, 1, 0, 1];
umat(:,15) = [1, 1, 1, 0];

umat(:,16) = [1, 1, 1, 1];

% associated weights
ws = zeros(4,1);
w = zeros(16,1);

for i = 1:16
    
    I = umat(:,i) == 1;
    ws(I) = uast;
    
    I = umat(:,i) == 0;    
    ws(I) = 1-uast;
    
    w(i) = prod(ws);
end

if sum(w) ~= 1
    error('probabilities does not sum to 1');
end

%% "Simulation"

x = zeros(2,16,8);

for t = 1:8
for i = 1:16 
    
    if t == 1
    
        x(:,i,t) = [1, 0];
    
    elseif t <= 4
        
        x(:,i,t) = A*x(:,i,t-1);
    
    else
        
        if umat(t-4, i) == 1
            x(:,i,t) = A_u*x(:,i,t-1);
        else
            x(:,i,t) = A_w*x(:,i,t-1);        
        end
        
    end
       
end;

    if abs(sum(sum(x(:,:,t)))-16) > 10^(-8)
        sum(x(:,:,t))
        error('probabilities does not sum to 1');
    end
    
end;


%% Results

x_final = x(:,:,end);
kept = x_final(1,:);
lost = x_final(2,:);

lost_tot = lost*w;

kept_nottreated = kept(1);
lost_nottreated = lost(1);

kept_treated = kept(2:end)*w(2:end) / sum(w(2:end));
lost_treated = lost(2:end)*w(2:end) / sum(w(2:end));

OR = (lost_treated/kept_treated) / (lost_nottreated/kept_nottreated);

[lost_tot, lost_nottreated, lost_treated, OR]