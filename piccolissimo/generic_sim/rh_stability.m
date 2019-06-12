syms a b c d e real
syms g positive
syms s

assumeAlso(a < 0);
assumeAlso(d < 0);

assumeAlso(b*e < 0);

%% Define u and v
u = g/s/(a-s); %without f
v = g/s/(s-a); %without f
% u = (f+g/s)/(a-s); %with f
% v = (f+g/s)/(s-a); %with f

% A = [c*v + d , b*u - H+e; b*v + H-e, -c*u + d];
% A = [c*v+d,b*u+e;b*v-e,-c*u+d]; % no H
% A = [d,b*u+e;b*v-e, d]; % no H, no c
% A = [c*v+d,0;0,-c*u+d]; % no H, no b, no e
% A = [c*v + d , b*u - H; b*v + H, -c*u + d]; % no e
% A = [d , b*u - H; b*v + H, d]; %no c and no e
% A = [a,0,0,0,0,-g;0,a,0,0,g,0;b,c,d,e,0,0;-c,b,-e,d,0,0;0,0,1,0,0,0;0,0,0,1,0,0];
% A = [a,0,0,0,0,-g;0,a,0,0,g,0;b,0,d,e,0,0;0,b,-e,d,0,0;0,0,1,0,0,0;0,0,0,1,0,0]; % no c
% A = [a,0,0,0,0,-g;0,a,0,0,g,0;b,0,0,e,0,0;0,b,-e,0,0,0;0,0,1,0,0,0;0,0,0,1,0,0]; % no c and no d
A = [0,0,0,0,0,-g;0,0,0,0,g,0;b,0,0,e,0,0;0,b,-e,0,0,0;0,0,1,0,0,0;0,0,0,1,0,0]; % no a, c, d

% A = [a, 0, -g; ...
%      c, d, 0; ...
%      0, 1, 0];

lambda = det(A - s*eye(length(A)));
lambda2 = eig(A);

%you'll probably have to multiply by some denominator, here's help to
%figure out what it is
lambda_simple = factor(lambda);
% pretty(lambda_simple)
disp('Lambda2 pretty');
pretty(lambda2)
disp('Lambda2 simplify pretty');
pretty(simplify(lambda2,1000))


% lambda_simple = factor(lambda*(s^2*(a-s)^2)); %do the multiply here

[alphas, sierras] = coeffs(lambda_simple,s);
alphas = simplify(alphas, 10000);
rh = routh_hurwitz(alphas);
% rh_simp = simplify(rh(:,1),10000);
% pretty(-alphas)
% pretty(-rh_simp)