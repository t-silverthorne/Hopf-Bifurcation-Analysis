function l10 = get_l10_autodiff_complex(Frhs,x0)
% construct B(x,y,ind) and C(x,y,z,ind)
for ind=1:numel(Frhs)
    [~,dydx(ind,:),dy2dx(:,:,ind),dy3dx(:,:,:,ind)]=dlfeval(@(x) get_diffmat_realpart(x,@(w) Frhs{ind}(w)),x0);
end
dy2dx=extractdata(dy2dx);
dy3dx=extractdata(dy3dx);
B=@(x,y,ind) x'*dy2dx(:,:,ind)*y;
C=@(x,y,z,ind) sum(arrayfun( @(kk) x'*dy3dx(:,:,kk,ind)*y*z(kk),1:numel(z)));
B_func=@(x,y)   arrayfun(@(ind) B(x,y,ind),transpose(1:numel(x0)));
C_func=@(x,y,z)   arrayfun(@(ind) C(x,y,z,ind),transpose(1:numel(x0)));

% convert Jacobian to matrix
Amat=extractdata(dydx);

% get q eigenvector
[ev,ea]=eig(Amat);
ea=diag(ea);
[~,ind]=max(imag(ea));
omega0=abs(imag(ea(ind)));
q=ev(:,ind);

% get p eigenvector
[ev,ea]=eig(Amat');
ea=diag(ea);
[~,ind]=min(imag(ea));
p=ev(:,ind);
p=p/conj(p'*q);

ABQ=Amat\B_func(q,q);
IABQ=(2*1i*omega0*eye(numel(x0))-Amat)\B_func(q,q);
l10=real(p'*C_func(q,q,q) - 2*p'*B_func(q,ABQ) + p'*B_func(q,IABQ))/2/omega0;
%%

function [y,dydx,dy2dx,dy3dx] = get_diffmat(x,Frhs)
% TODO: figure out how to pre-allocate a dlarray
n=numel(x);
y = Frhs(x);
dydx = dlgradient(y,x,'EnableHigherDerivatives',true);
for ii=1:n
    dy2dx(ii,1:n)= dlgradient(dydx(ii),x,'EnableHigherDerivatives',true);
end
for ii=1:n
    for jj=1:n
        dy3dx(ii,jj,:)=dlgradient(dy2dx(ii,jj),x);
    end
end
end


function [y,dydx,dy2dx,dy3dx] = get_diffmat_realpart(x,Frhs)
% TODO: figure out how to pre-allocate a dlarray
n=numel(x);
y = Frhs(x);
dydx = real(dlgradient(y,x,'EnableHigherDerivatives',true));
for ii=1:n
    dy2dx(ii,1:n)= real(dlgradient(dydx(ii),x,'EnableHigherDerivatives',true));
end
for ii=1:n
    for jj=1:n
        dy3dx(ii,jj,:)=real(dlgradient(dy2dx(ii,jj),x));
    end
end
end
end