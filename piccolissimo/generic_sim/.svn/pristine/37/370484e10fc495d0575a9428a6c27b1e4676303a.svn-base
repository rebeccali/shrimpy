%newton/euler
syms nu u v w p q r omg_r omg_s phi theta psi psi_r psi_s xe ye ze i F_p_x F_p_y F_p_z F_d_x F_d_y F_d_z Is_sx Is_sy Is_sz Ir_rx Ir_ry Ir_rz h_s h_r m_s m_r M_p_x M_p_y M_p_z M_d_x M_d_y M_d_z g omg_sDot omg_rDot real
Is_s = [Is_sx 0 0; 0 Is_sy 0; 0 0 Is_sz];
Ir_r = [Ir_rx 0 0; 0 Ir_ry 0; 0 0 Ir_rz];
m = m_s + m_r;
S_s_f = [0 0 h_s];
S_r_f = [0 0 h_r];
X = [nu; u; v; w; p; q; r; omg_r; omg_s; phi; theta; psi; psi_r; psi_s; xe; ye; ze; i];
% V = [u;v;w];
V = [u;v;0];
F_p = [F_p_x; F_p_y; F_p_z];
F_d = [F_d_x; F_d_y; F_d_z];
M_p = [M_p_x; M_p_y; M_p_z];
M_d = [M_d_x; M_d_y; M_d_z];

V_dot = (F_p + F_d + m*[0; 0; g] - m*cross([p; q; r],V))/m

omg_dot = (Is_s + Ir_r + m_s*(S_s_f')*S_s_f + m_r*(S_r_f')*S_r_f)\...
    (M_d + M_p - cross([p;q;r],Is_s*[p;q;r+omg_s] )- cross([p;q;r],Ir_r*[p;q;r+omg_r]) ...
    - Is_s*[0;0;omg_sDot] -Ir_r*[0;0;omg_rDot])

%blade
syms A h_p d rho c CL CD real
% V_A = V - [0;0;nu] + cross([p;q;r+omg_r],[0; d; h_p])
V_A = [u;v;0] - [0;0;nu] + cross([p;q;omg_r],[0; d; h_p])
%u, v, p, q negative on other side
%rot 90 deg u>-v, v>u, p>-q, q>p
V_A = [V_A(1); 0; V_A(3)];
alpha = V_A(3)/V_A(1); %small angle approx
lift = rho/2*dot(V_A,V_A)*c*CL*alpha;
drag = rho/2*dot(V_A,V_A)*c*CD*alpha;
% N = lift*cos(atan(V_A(3)/V_A(1))) + drag*sin(atan(V_A(3)/V_A(1)));
% T = -lift*sin(atan(V_A(3)/V_A(1)))+ drag*cos(atan(V_A(3)/V_A(1)));
% N = lift*cos(atan(alpha)) + drag*sin(atan(alpha));
% T = -lift*sin(atan(alpha))+ drag*cos(atan(alpha));
N = lift; %drag is small, so we don't care
T = drag + lift*alpha;
f1 = [-T; 0; -N];
m1 = cross([0; d; h_p],f);

% V_A = V - [0;0;nu] + cross([p;q;r+omg_r],[0; d; h_p])
V_A = [-u;-v;0] - [0;0;nu] + cross([-p;-q;omg_r],[0; d; h_p])
%u, v, p, q negative on other side
%rot 90 deg u>-v, v>u, p>-q, q>p
V_A = [V_A(1); 0; V_A(3)];
alpha = V_A(3)/V_A(1); %small angle approx
lift = rho/2*dot(V_A,V_A)*c*CL*alpha;
drag = rho/2*dot(V_A,V_A)*c*CD*alpha;
% N = lift*cos(atan(V_A(3)/V_A(1))) + drag*sin(atan(V_A(3)/V_A(1)));
% T = -lift*sin(atan(V_A(3)/V_A(1)))+ drag*cos(atan(V_A(3)/V_A(1)));
% N = lift*cos(atan(alpha)) + drag*sin(atan(alpha));
% T = -lift*sin(atan(alpha))+ drag*cos(atan(alpha));
N = lift; %drag is small, so we don't care
T = drag + lift*alpha;
f2 = [T; 0; -N];
m2 = cross([0; d; h_p],f);

% V_A = V - [0;0;nu] + cross([p;q;r+omg_r],[0; d; h_p])
V_A = [-v;u;0] - [0;0;nu] + cross([-q;p;omg_r],[0; d; h_p])
%u, v, p, q negative on other side
%rot 90 deg u>-v, v>u, p>-q, q>p
V_A = [V_A(1); 0; V_A(3)];
alpha = V_A(3)/V_A(1); %small angle approx
lift = rho/2*dot(V_A,V_A)*c*CL*alpha;
drag = rho/2*dot(V_A,V_A)*c*CD*alpha;
% N = lift*cos(atan(V_A(3)/V_A(1))) + drag*sin(atan(V_A(3)/V_A(1)));
% T = -lift*sin(atan(V_A(3)/V_A(1)))+ drag*cos(atan(V_A(3)/V_A(1)));
% N = lift*cos(atan(alpha)) + drag*sin(atan(alpha));
% T = -lift*sin(atan(alpha))+ drag*cos(atan(alpha));
N = lift; %drag is small, so we don't care
T = drag + lift*alpha;
f3 = [0; T; -N];
m3 = cross([d; 0; h_p],f);

% V_A = V - [0;0;nu] + cross([p;q;r+omg_r],[0; d; h_p])
V_A = [v;-u;0] - [0;0;nu] + cross([q;-p;omg_r],[0; d; h_p])
%u, v, p, q negative on other side
%rot 90 deg u>-v, v>u, p>-q, q>p
V_A = [V_A(1); 0; V_A(3)];
alpha = V_A(3)/V_A(1); %small angle approx
lift = rho/2*dot(V_A,V_A)*c*CL*alpha;
drag = rho/2*dot(V_A,V_A)*c*CD*alpha;
% N = lift*cos(atan(V_A(3)/V_A(1))) + drag*sin(atan(V_A(3)/V_A(1)));
% T = -lift*sin(atan(V_A(3)/V_A(1)))+ drag*cos(atan(V_A(3)/V_A(1)));
% N = lift*cos(atan(alpha)) + drag*sin(atan(alpha));
% T = -lift*sin(atan(alpha))+ drag*cos(atan(alpha));
N = lift; %drag is small, so we don't care
T = drag + lift*alpha;
f4 = [0; -T; -N];
m4 = cross([d; 0; h_p],f);

f = f1+f2+f3+f4
m = m1+m2+m3+m4