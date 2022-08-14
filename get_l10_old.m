function l10=get_l10_old(gamma,alpha)
% Calculates Lyapunov coefficient using method from Kusznetsov's textbook
num_var=2; % CHANGE
% currently implemented for 2 variable model, you need to change this
% manually for higher dimensional systems
x = sym('x',[1 num_var],'real');
y = sym('y',[1 num_var],'real');
z = sym('z',[1 num_var],'real');
X=sym2cell(x);
Y=sym2cell(y);
Z=sym2cell(z);
% specify RHS of model and equilibrium of interest (could make these
% arguments of function)
F=[1-y(1)*y(2)^gamma; alpha*y(2)*(y(1)*y(2)^(gamma-1)-1)]; % CHANGE
y_eq=[1 1]; % CHANGE
A_func=jacobian(F,y);

id=eye(num_var);

% Get B(x,y)
B_vec=sym(zeros(num_var,1));
for i=1:num_var
    Fi=sum(F.*id(:,i));
    for j=1:num_var
        for k=1:num_var
            B_vec(i) = B_vec(i) + subs(diff(diff(Fi,y(j)),y(k)),y,y_eq)*x(j)*y(k);
        end
    end
end
B_func(X{:},Y{:})=B_vec;

% Get C(x,y,z)
C_vec=sym(zeros(num_var,1));
for i=1:num_var
    Fi=sum(F.*id(:,i));
    for j=1:num_var
        for k=1:num_var
            for l=1:num_var
                C_vec(i) = C_vec(i) + subs(diff(diff(diff(Fi,y(j)),y(k)),y(l)),y,y_eq)*x(j)*y(k)*z(l);
            end
        end
    end
end
C_func(X{:},Y{:},Z{:})=C_vec;

% make sure to set this to your equilibrium values
A_numerical=subs(A_func,y,y_eq)

% get q eigenvector
[ev,ea]=eig(A_numerical);
ea=diag(ea);
[~,ind]=max(imag(ea));
omega0=abs(imag(ea(ind)));
q=ev(:,ind);

% get p eigenvector
[ev,ea]=eig(A_numerical');
ea=diag(ea);
[~,ind]=min(imag(ea));
p=ev(:,ind);
p=p/conj(p'*q);

Q=num2cell(q);
Qp=num2cell(q');

% quantities needed for Lyapunov calculation
% make sure to modify these so that they include q1,q2,...,qnumvar
ABq=A_numerical\B_func(Q{:},Qp{:});
ABQ=num2cell(ABq);
IABq=(2*1i*omega0*eye(num_var)-A_numerical)\B_func(Q{:},Q{:});
IABQ=num2cell(IABq);

% make sure to evaluate at equilibrium for your model
l10=real(subs(p'*C_func(Q{:},Q{:},Qp{:})...
    - 2*p'*B_func(Q{:},ABQ{:})...
    +p'*B_func(Qp{:},IABQ{:}),y,y_eq))/2/omega0;
end