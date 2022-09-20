% Detect Hopf bifurcations using predictor-corrector method on the IT-MFL model 
% Classify them using the auto diff implementation of the Lyapunov coefficient
clear
rho=0
K2=0
K2min = 0; K2max = 10;
Atot=1.01*1e-3;
gtot=1e-4;
Kd=10^(-6.73);
K1=1/Atot;
[K2,Btot,Aqss] = K2search(Atot,gtot,Kd,K1,rho,K2min,K2max);
computeLyapunovCoefficient(Btot,Atot,gtot,Kd,K1,K2,rho)
%%
close all

warning('on')
findHopf_continuation()
%%
fnamestr='hopf_fig_updated';
% generic template for nicely formatted matlab plot
fh=gcf
plot_filename=fnamestr;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
% save
ht=3.4; % height
wd=4; % width
set(gcf,'PaperUnits','inches')
set(gcf,'PaperPositionMode','manual','PaperSize',[wd,ht],'PaperPosition',[0 0 wd ht])
print(gcf,plot_filename,'-dpng','-r600') % -r sets the resolution
savefig(gcf,strcat(plot_filename,'.fig'))% save matlab .fig too
function [K1_values, rho_HB, Btoteq_HB, K2_HB, l10] = findHopf_continuation()
    display('called')
    % parameters
	Atot=1.31e-1;
    gtot=1.31e-2;
    Kd=1e-5;
%    Atot=1e-3;
%    gtot=1e-4;
%    Kd=10^(-6.73);
    
    numK1=1e2;
    K1min = 1e-1; K1max = 1e3;
    K1_values = logspace(log10(K1min),log10(K1max),numK1);
    
    K2min = 0; K2max = 2;
    rhomin = 0; rhomax = 1;
    
    Aqss_HB = cell(numK1,1);
    Btoteq_HB = cell(numK1,1);
    rho_HB = cell(numK1,1);
    K2_HB = cell(numK1,1);
    l10_HB = cell(numK1,1);
    
    h=1e-2; % parameterization step length
    
    for K1i=1:numK1
        K1 = K1_values(K1i);
        
        % initial point
        rho = rhomax;
        [K2,Btot,Aqss] = K2search(Atot,gtot,Kd,K1,rho,K2min,K2max);
        if isnan(K2)
            continue;
        end
        l10=computeLyapunovCoefficientOLD(Btot,Atot,gtot,Kd,K1,K2,rho);
        
        step = 1;
        Aqss_HB{K1i}(step) = Aqss;
        Btoteq_HB{K1i}(step) = Btot;
        rho_HB{K1i}(step) = rho;
        K2_HB{K1i}(step) = K2;
        l10_HB{K1i}(step) = l10;
        
        % continue the curve
        done = false;
        while ~done
            display(step)
            step = step + 1;
            
            % predict
            if step<=2
                rho = rho - h;
            else
                tau  = [ diff( K2_HB{K1i}([step-2,step-1]) ), diff( rho_HB{K1i}([step-2,step-1]) ) ];
                tau = tau / norm(tau,2);
                K2   = K2  + h*tau(1);
                rho  = rho + h*tau(2);
            end
            
            done = rho<rhomin || rho>rhomax || isnan(K2); % || K2<K2min || K2>K2max
            if done
                break;
            end
            
            % correct
            if step<=2
                [K2,Btot,Aqss] = K2search(Atot,gtot,Kd,K1,rho,K2min,K2max);
            else
                % use fsolve to find a solution
                Btot = fzero(@(Btot) T_fun(Btot,Atot,gtot,Kd,K1,K2,rho) - Btot, Atot); % update Btot to be an equilibrium
                options = optimset('Display','off');
                [x,~,EXITFLAG]=fsolve( @(x) [Phifun(x(1),x(2),Atot,gtot,Kd,K1,x(3));...
                                             T_fun(x(2),Atot,gtot,Kd,K1,x(3),x(4)) - x(2);...
                                            -(-nth_output( 2, @T_fun, x(2),Atot,gtot,Kd,K1,x(3),x(4))).^(1/3)+2],...
                                      [Aqss,Btot,K2,rho],options);
                Aqss=x(1); Btot=x(2); K2=x(3); rho=x(4);
            end
            
            % evaluate l1(0)
            l10=computeLyapunovCoefficientOLD(Btot,Atot,gtot,Kd,K1,K2,rho);
            
            % save the point
            Aqss_HB{K1i}(step) = Aqss;
            Btoteq_HB{K1i}(step) = Btot;
            rho_HB{K1i}(step) = rho;
            K2_HB{K1i}(step) = K2;
            l10_HB{K1i}(step) = l10;
            
            %figure(42); hold on;
            %plot(K2,rho,'.k')
            %drawnow

            done = rho<rhomin || rho>rhomax || isnan(K2); % || K2<K2min || K2>K2max
        end
        figure(42); 
        set(gca, 'XScale', 'log', 'YScale', 'log');
        loglog(K2_HB{K1i},rho_HB{K1i},'-k')
        
        hold on;
        grid on
        drawnow
        axis([K2min K2max rhomin rhomax])
    end
    
    figure(); hold on;
    lw = 2;
    clrs = gray(numK1+1); clrs = clrs(1:end-1,:);
    for K1i=numK1:-1:1
        loglog(K2_HB{K1i},rho_HB{K1i},'-', 'Color', clrs(K1i,:), 'LineWidth',lw)
    end
    xlabel('$$K_2$$','interpreter','latex')
    ylabel('$$\rho$$','interpreter','latex')
    box on;
    set(gca,'TickLabelInterpreter','latex')
    axis([K2min K2max rhomin rhomax])
    cbh=colorbar;
    ylabel(cbh,'$$\tilde K_1$$','interpreter','latex')
    set(cbh,'TickLabelInterpreter','latex')
    set(gca,'ColorScale','log')
    caxis([K1min,K1max])
    colormap(clrs)
    
end

function [K2,Btot,Aqss] = K2search(Atot,gtot,Kd,K1,rho,K2min,K2max)
    
    Btoteq_fun = @(K2) fzero(@(Btot) T_fun(Btot,Atot,gtot,Kd,K1,K2,rho) - Btot, Atot);
    Tp_fun = @(K2) nth_output( 2, @T_fun, Btoteq_fun(K2),Atot,gtot,Kd,K1,K2,rho);
    realeigzero_fun = @(K2) -(-Tp_fun(K2)).^(1/3)+2;
    if realeigzero_fun(K2min)*realeigzero_fun(K2max)<0 % likely no Hopf bifurcations for this value of K1
        K2 = fzero( realeigzero_fun, 0); % find the value of K2 corresponding to the Hopf bifurcation
        Btot = Btoteq_fun(K2);
        Aqss = Aqssfun(Btot,Atot,gtot,Kd,K1,K2);
    else
        K2=NaN; Btot=NaN; Aqss=NaN;
    end
end

% function l10 = computeLyapunovCoefficient(Btoteq,Atot,gtot,Kd,K1,K2,rho)
% x0=dlarray([Btoteq;Btoteq;Btoteq]);
% Frhs1=@(x) real(K1/((1+K2)/( ((K2 + K1*gtot - Atot*K1 + x(3)*K1 + Kd*K1 + 1)^2/(9*K1^2) - (x(3) - Atot + Kd - Atot*K2 + x(3)*K2 + K2*Kd - Atot*Kd*K1 + Kd*K1*gtot)/(3*K1))/((Atot*Kd + Atot*K2*Kd)/(2*K1) + (((Atot*Kd + Atot*K2*Kd)/(2*K1) - (K2 + K1*gtot - Atot*K1 + x(3)*K1 + Kd*K1 + 1)^3/(27*K1^3) + ((K2 + K1*gtot - Atot*K1 + x(3)*K1 + Kd*K1 + 1)*(x(3) - Atot + Kd - Atot*K2 + x(3)*K2 + K2*Kd - Atot*Kd*K1 + Kd*K1*gtot))/(6*K1^2))^2 - ((K2 + K1*gtot - Atot*K1 + x(3)*K1 + Kd*K1 + 1)^2/(9*K1^2) - (x(3) - Atot + Kd - Atot*K2 + x(3)*K2 + K2*Kd - Atot*Kd*K1 + Kd*K1*gtot)/(3*K1))^3)^(1/2) - (K2 + K1*gtot - Atot*K1 + x(3)*K1 + Kd*K1 + 1)^3/(27*K1^3) + ((K2 + K1*gtot - Atot*K1 + x(3)*K1 + Kd*K1 + 1)*(x(3) - Atot + Kd - Atot*K2 + x(3)*K2 + K2*Kd - Atot*Kd*K1 + Kd*K1*gtot))/(6*K1^2))^(1/3) + ((Atot*Kd + Atot*K2*Kd)/(2*K1) + (((Atot*Kd + Atot*K2*Kd)/(2*K1) - (K2 + K1*gtot - Atot*K1 + x(3)*K1 + Kd*K1 + 1)^3/(27*K1^3) + ((K2 + K1*gtot - Atot*K1 + x(3)*K1 + Kd*K1 + 1)*(x(3) - Atot + Kd - Atot*K2 + x(3)*K2 + K2*Kd - Atot*Kd*K1 + Kd*K1*gtot))/(6*K1^2))^2 - ((K2 + K1*gtot - Atot*K1 + x(3)*K1 + Kd*K1 + 1)^2/(9*K1^2) - (x(3) - Atot + Kd - Atot*K2 + x(3)*K2 + K2*Kd - Atot*Kd*K1 + Kd*K1*gtot)/(3*K1))^3)^(1/2) - (K2 + K1*gtot - Atot*K1 + x(3)*K1 + Kd*K1 + 1)^3/(27*K1^3) + ((K2 + K1*gtot - Atot*K1 + x(3)*K1 + Kd*K1 + 1)*(x(3) - Atot + Kd - Atot*K2 + x(3)*K2 + K2*Kd - Atot*Kd*K1 + Kd*K1*gtot))/(6*K1^2))^(1/3) - (K2 + K1*gtot - Atot*K1 + x(3)*K1 + Kd*K1 + 1)/(3*K1)  )+K1) * (1-rho*K2/(1+K2)) + rho*K2/(1+K2)-x(1));
% Frhs2=@(x) x(1)-x(2);
% Frhs3=@(x) x(2)-x(3);
% Frhs={Frhs1,Frhs2,Frhs3};
% l10=get_l10_autodiff_complex(Frhs,x0);
% end

function l10 = computeLyapunovCoefficient(Btoteq,Atot,gtot,Kd,K1,K2,rho)
x0=dlarray([Btoteq;Btoteq;Btoteq]);
Aqss=@(x) real(((K2 + K1*gtot - Atot*K1 + x(3)*K1 + Kd*K1 + 1)^2/(9*K1^2) - (x(3) - Atot + Kd - Atot*K2 + x(3)*K2 + K2*Kd - Atot*Kd*K1 + Kd*K1*gtot)/(3*K1))/((Atot*Kd + Atot*K2*Kd)/(2*K1) + (((Atot*Kd + Atot*K2*Kd)/(2*K1) - (K2 + K1*gtot - Atot*K1 + x(3)*K1 + Kd*K1 + 1)^3/(27*K1^3) + ((K2 + K1*gtot - Atot*K1 + x(3)*K1 + Kd*K1 + 1)*(x(3) - Atot + Kd - Atot*K2 + x(3)*K2 + K2*Kd - Atot*Kd*K1 + Kd*K1*gtot))/(6*K1^2))^2 - ((K2 + K1*gtot - Atot*K1 + x(3)*K1 + Kd*K1 + 1)^2/(9*K1^2) - (x(3) - Atot + Kd - Atot*K2 + x(3)*K2 + K2*Kd - Atot*Kd*K1 + Kd*K1*gtot)/(3*K1))^3)^(1/2) - (K2 + K1*gtot - Atot*K1 + x(3)*K1 + Kd*K1 + 1)^3/(27*K1^3) + ((K2 + K1*gtot - Atot*K1 + x(3)*K1 + Kd*K1 + 1)*(x(3) - Atot + Kd - Atot*K2 + x(3)*K2 + K2*Kd - Atot*Kd*K1 + Kd*K1*gtot))/(6*K1^2))^(1/3) + ((Atot*Kd + Atot*K2*Kd)/(2*K1) + (((Atot*Kd + Atot*K2*Kd)/(2*K1) - (K2 + K1*gtot - Atot*K1 + x(3)*K1 + Kd*K1 + 1)^3/(27*K1^3) + ((K2 + K1*gtot - Atot*K1 + x(3)*K1 + Kd*K1 + 1)*(x(3) - Atot + Kd - Atot*K2 + x(3)*K2 + K2*Kd - Atot*Kd*K1 + Kd*K1*gtot))/(6*K1^2))^2 - ((K2 + K1*gtot - Atot*K1 + x(3)*K1 + Kd*K1 + 1)^2/(9*K1^2) - (x(3) - Atot + Kd - Atot*K2 + x(3)*K2 + K2*Kd - Atot*Kd*K1 + Kd*K1*gtot)/(3*K1))^3)^(1/2) - (K2 + K1*gtot - Atot*K1 + x(3)*K1 + Kd*K1 + 1)^3/(27*K1^3) + ((K2 + K1*gtot - Atot*K1 + x(3)*K1 + Kd*K1 + 1)*(x(3) - Atot + Kd - Atot*K2 + x(3)*K2 + K2*Kd - Atot*Kd*K1 + Kd*K1*gtot))/(6*K1^2))^(1/3) - (K2 + K1*gtot - Atot*K1 + x(3)*K1 + Kd*K1 + 1)/(3*K1));
T=@(x) real(K1/((1+K2)/Aqss(x)+K1) * (1-rho*K2/(1+K2)) + rho*K2/(1+K2));
Frhs1=@(x) T(x)-x(1);
Frhs2=@(x) x(1)-x(2);
Frhs3=@(x) x(2)-x(3);
Frhs={Frhs1,Frhs2,Frhs3};
l10=get_l10_autodiff_complex(Frhs,x0);
end


function l10 = computeLyapunovCoefficientOLD(Btoteq,Atot,gtot,Kd,K1,K2,rho)
    % compute the first Lyapunov coefficient
        [T,Tp,Tpp,Tppp] = T_fun(Btoteq,Atot,gtot,Kd,K1,K2,rho);
        L = [-1,0,Tp;1,-1,0;0,1,-1]; % linearization
        [ev,ea]=eig(L);
        ea = diag(ea);
        [~,ind] = max( imag(ea) );
        %omega0 = imag(ea(ind));
        omega0 = (-Tp).^(1/3)*sqrt(3)/2;
        q = ev(:,ind);
        %if max(abs(L*q-1i*omega0*q)) > 1e-10
        %    error('eigenvector q not correct')
        %end

        [ev,ea]=eig(transpose(L));
        ea = diag(ea);
        [~,ind] = min( abs(ea+1i*omega0) );
        p = ev(:,ind);
        %if max(abs(transpose(L)*p+1i*omega0*p)) > 1e-10
        %    error('eigenvector p not correct')
        %end
        p = p / conj(p'*q);

        AinvB = L\[Tpp*q(3)*q(3)';0;0];
        val = (2*1i*omega0*eye(3,3)-L)\[Tpp*q(3)^2;0;0];
        l10 = real( ...
                                 p'*[Tppp*(q(3)^2*q(3)');0;0] ...
                              -2*p'*[Tpp*q(3)*AinvB(3);0;0] ...
                                +p'*[Tpp*q(3)'*val(3);0;0]...
                             )/2/omega0;
end

function [T,Tp,Tpp,Tppp] = T_fun(Btot,Atot,gtot,Kd,K1,K2,rho)
% evaluated the transcription function as well as its derivatives with
% respect to Btot

    [Aqss,Apqss,Appqss,Apppqss] = Aqssfun(Btot,Atot,gtot,Kd,K1,K2);
    T    = K1/((1+K2)/Aqss+K1) * (1-rho*K2/(1+K2)) + rho*K2/(1+K2);

    Tp   = K1*(1+K2-rho*K2)/(1+K2+K1*Aqss)^2*Apqss;
    Tpp  = K1*(1+K2-rho*K2)/(1+K2+K1*Aqss)^2*Appqss - 2*K1^2*(1+K2-rho*K2)/(1+K2+K1*Aqss)^3*Apqss^2 ;
    Tppp = (K1*(K2 - K2*rho + 1)*(2*K2*Apppqss + 6*K1^2*Apqss^3 + K2^2*Apppqss - 6*K1*Apqss*Appqss + K1^2*Aqss^2*Apppqss + 2*K1*Aqss*Apppqss - 6*K1^2*Aqss*Apqss*Appqss + 2*K1*K2*Aqss*Apppqss - 6*K1*K2*Apqss*Appqss + Apppqss))/(K2 + K1*Aqss + 1)^4; % TODO: simplify this
    
    %syms K1 K2 rho Btot Aqss(Btot)
    %T    = K1/((1+K2)/Aqss(Btot)+K1) * (1-rho*K2/(1+K2)) + rho*K2/(1+K2);
    %for n=1:3
    %    simplify(diff(T,Btot,n))
    %end
end

function val = Phifun(Aqss,Btot,Atot,gtot,Kd,K1,K2) 
    % evaluates the cubic for Aqss
    
    a = K1;
    b = (1+K2) + (gtot+Btot-Atot+Kd)*K1;
    c = (1+K2)*(Btot-Atot+Kd)+ (gtot-Atot)*Kd*K1;
    d = -(1+K2)*Atot*Kd;
    
    val = polyval([a,b,c,d],Aqss);
    
end

function [Aqss,Apqss,Appqss,Apppqss] = Aqssfun(Btot,Atot,gtot,Kd,K1,K2) 
    % function to solve the cubic for Aqss
    % also returns the derivative of Aqq with respect to Btot
    
    a = K1;
    b = (1+K2) + (gtot+Btot-Atot+Kd)*K1;
    c = (1+K2)*(Btot-Atot+Kd)+ (gtot-Atot)*Kd*K1;
    d = -(1+K2)*Atot*Kd;
    
    Aqss = roots( [a,b,c,d] );
    Aqss = Aqss( isreal(Aqss) & real(Aqss)>0 );
    
    if nargout > 1
        bp = K1;
        cp = 1+K2;
        Apqss = -(cp*Aqss+bp*Aqss^2)/(3*a*Aqss^2+2*b*Aqss+c);
    end
    if nargout > 2
        Appqss = -( (2*b+6*a*Aqss)*Apqss^2 + (2*cp+4*bp*Aqss)*Apqss ) /(3*a*Aqss^2+2*b*Aqss+c);
    end
    if nargout > 3
        Apppqss = -( 6*a*Apqss^3 + 6*bp*Apqss^2 + 6*(b+3*a*Aqss)*Appqss*Apqss + 3*(cp+2*bp*Aqss)*Appqss ) /(3*a*Aqss^2+2*b*Aqss+c);
    end
    
    % symbolically compute these derivatives
    %syms B a b(B)  d A(B)
    %syms bp cp b0 c0
    %b = B*bp + b0; % linear dependence bp = K1
    %c = B*cp + c0; % linear dependence cp = 1+K2
    %Phi = a*A^3+b*A^2+c*A+d;
    %for n=1:3
    %    collect( simplify(diff(Phi,B,n))==0, diff(A(B), B, n) )
    %end
end

function value = nth_output(N,fcn,varargin)
  [value{1:N}] = fcn(varargin{:});
  value = value{N};
end