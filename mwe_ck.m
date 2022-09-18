% pdepe solution check pannel

% hard-coded parameters %
kae=24.9194653795352e+000;
kah=23.2581676875662e+000;
Cn1=773.634541436672e-006;
Cp1=773.634541436672e-006;
h1 =7.54293677900755e+000;
te =61.7330028639972e-030;
th =66.1425030685685e-030;
fswpNc =3.23975816372672e+021;
fswpNv =3.02377428614494e+021;
% --------------------- %

figure

for i=2:size(sol,1)

    xq = x;

    %1. df/dx
        [u1i,dudx1i] = pdeval(m,x,sol(i,:,1),xq);
        [u2i,dudx2i] = pdeval(m,x,sol(i,:,2),xq);
    
        f1 = Cn1.^-1.*( (f_ext(u1i) + kae.*u1i).*dudx1i - kah.*u1i.*dudx2i );
        f2 = Cp1.^-1.*( (f_ext(u2i) + kah.*u2i).*dudx2i - kae.*u2i.*dudx1i );
    
        [~,df1dx] = pdeval(m,x,f1,xq);
        [~,df2dx] = pdeval(m,x,f2,xq);
    
        dfdx = [df1dx;df2dx];

    %2. source term
        u=[u1i;u2i];    
        s=u;
        for j=1:length(x)
            F=x;
            F(j) = REM(u(1,j),u(2,j),-h1,te,th);
            s(:,j) = [-fswpNc.^-1.*F(j);
                      -fswpNv.^-1.*F(j)];
        end
        dfdx_plus_s = dfdx + s;

    %3. du/dt via backward difference formula
        dudt_bwd = [ (t(i)-t(i-1)).^-1.*(u1(i,:)-u1(i-1,:));
                     (t(i)-t(i-1)).^-1.*(u2(i,:)-u2(i-1,:))];

    % compare du/dt and df/dx + s
    for z=1:2
    subplot(1,2,z)
    title(['$i=$ ',num2str(z),', $\tau = $',num2str(t(i))],"Interpreter","latex",'FontSize',14)
    xlabel('$\lambda$','Interpreter','latex','FontSize',14)
    yyaxis left
    plot(x,dudt_bwd(z,:),'-')%du/dt
    if z==1
        ylabel('$\partial_{\tau}{u_i}$','Interpreter','latex','FontSize',14)
    end
    yyaxis right
    plot(x,dfdx_plus_s(z,:),'*-')%df/dx + s
    if z==2
        ylabel('$c_i^{-1}\partial_{\lambda}{f_i}+\alpha_i^{-1}R$','Interpreter','latex','FontSize',14)
    end
    xlim([0.01,0.99])
    end

    pause(0.05)

end


%% function define
% copy from the main file
function y=f_ext(x)

    %            t
    % y =  -------------
    %       1 - exp(-t)

    abs_tol=5e-7;
    
    Bn=[1 1/2 1/12 -1/720 1/30240 -1/1209600];
    f1=@(t) Bn(1)+Bn(2).*t;
    f_o=@(t) t./(1-exp(-t));
    
    y=f_o(x);
    x_id=abs(x)<abs_tol;
    x_sp=x(x_id);
    if sum(x_id)~=0
        y(x_id)=f1(x_sp);
    end

end

function r=REM(u,v,h,mu1,mu2)

    %                       uv - Ff(u)Ff(v)exp(2h) 
    % r =  -------------------------------------------------------
    %       mu2 * ( u + Ff(u)exp(h) ) + mu1 * ( v + Ff(v)exp(h) )

    % recombination 
    r=(u.*v-Ff(u).*Ff(v).*exp(2.*h)).*...
    (df(mu2,u,h)+df(mu1,v,h)).^-1;

    function y=df(mu0,u0,h)
        y=mu0.*(u0+Ff(u0).*exp(h));
    end
    
    function y=Ff(x)
        y=x.*(exp(x)-1).^-1;
        idx=(abs(x)<1e-3);
        y(idx)=1-x(idx)./2+x(idx).^2./12;
    end

end
