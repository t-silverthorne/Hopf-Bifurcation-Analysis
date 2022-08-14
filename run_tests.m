clear
%% test case from Kuznetsov
% sign should agree, magnitude should not (unless you manually pick same
% eigenvectors as in the textbook)
x0=dlarray([0;0;0]);
beta=.1;
alpha=1/beta+10;
Frhs1=@(x) x(2);
Frhs2=@(x) x(3);
Frhs3=@(x) -alpha*x(3) - beta*x(2) -x(1) +x(1).^2;
Frhs={Frhs1,Frhs2,Frhs3};
get_l10_autodiff_complex(Frhs,x0)
-(1+8*beta^3)*beta*sqrt(beta)/(1+4*beta^3)/(1+beta^3)
