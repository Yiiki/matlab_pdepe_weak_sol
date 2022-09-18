% minimum work example
% stimulus settings
Vamp=0.5;
vdsfx=1.5;
vgsf=@(t) Vamp.*t.^0;
% vdsf=@(t) vdsfx.*0.5.*(1.001+Eapp(t));
vdsf=@(t) vdsfx.*0.5.*(1.001+sin(2.*pi.*t));

% solution settings
x = linspace(0,1,101);
t = linspace(0,2,201);
m = 0;

% pdepe solver
pdeic = @(x) pdeic_func(x,vgsf(0),vdsf(0));
pdebc = @(xl,ul,xr,ur,t) pdebc_func(xl,ul,xr,ur,t,vgsf(t),vdsf(t));
sol = pdepe(m,@pdefun,pdeic,pdebc,x,t);
u1 = sol(:,:,1);% u
u2 = sol(:,:,2);% v


%% function define

function [c,f,s] = pdefun(~,~,u,dudx)
    
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
    
    c = [1; 1];
    
    f1 = Cn1.^-1.*((f_ext(u(1)) + kae.*u(1)).*dudx(1) - kah.*u(1).*dudx(2));
    f2 = Cp1.^-1.*((f_ext(u(2)) + kah.*u(2)).*dudx(2) - kae.*u(2).*dudx(1));
    
    f = [f1;f2];
    
    F = REM(u(1),u(2),-h1,te,th);
    
    s = [-fswpNc.^-1.*F;
         -fswpNv.^-1.*F];

end


function uv_mat = pdeic_func(x,vgsf,vdsf)

    % - -  - - - - - -- - - - - - +
    % hard-coded parameters %   % |
    pn=39.0000000000000e-003;   % +
    pp=351.000000000000e-003;   % |
    VT=25.8519997864355e-003;   % +
    kae=24.9194653795352e+000;  % |
    kah=23.2581676875662e+000;  % +
    % --------------------------+ |
    
    vdsl=0;
    
    vdsr=vdsf;
    
    Vgs=vgsf;
    
    uv_mat=kron([0 0]',x);
    
    for i=1:size(x,2) % x should be a row vector
    
        xi=x(i);
    
        vdsx=interp1([0,1],[vdsl,vdsr],xi);
        
        xi_mat=([vdsx+pn;...          % xi_n
                 vdsx-pp]-Vgs)./VT;   % xi_p
                 
        xi_vec=sol_bnd(kae,kah,xi_mat(1),xi_mat(2));
        
        uv_vec=[log_exp_plus(xi_vec-xi_mat(1));... % u
                log_exp_plus(xi_mat(2)-xi_vec)];...   % v
        
        uv_mat(:,i)=uv_vec;
    end

end


function [pl,ql,pr,qr]=pdebc_func(~,ul,~,ur,~,vgsf,vdsf)
    
    % - -  - - - - - -- - - - - - +
    % hard-coded parameters %   % |
    pn=39.0000000000000e-003;   % +
    pp=351.000000000000e-003;   % |
    VT=25.8519997864355e-003;   % +
    kae=24.9194653795352e+000;  % |
    kah=23.2581676875662e+000;  % +
    % --------------------------+ |


    %     | x=0    x=1   |
    xi_mat=([pn,  vdsf+pn;...          % xi_n
             -pp, vdsf-pp]-vgsf)./VT;   % xi_p

    xi0=sol_bnd(kae,kah,xi_mat(1,1),xi_mat(2,1));
    
    psi0=xi0.*VT + vgsf;

    u0 = log_exp_plus((psi0-pn).*VT.^-1);
    v0 = log_exp_plus((-psi0-pp).*VT.^-1);

    xi1=sol_bnd(kae,kah,xi_mat(1,2),xi_mat(2,2));
    
    psi=xi1.*VT + vgsf;

    u1 = log_exp_plus((psi-vdsf-pn).*VT.^-1);
    v1 = log_exp_plus((-psi+vdsf-pp).*VT.^-1);
    
    ql = [0; 0];
    qr = [0; 0];
    pl = [0; 0];
    pr = [0; 0];
    pl(1)=ul(1)-u0;
    pl(2)=ul(2)-v0;
    pr(1)=ur(1)-u1;
    pr(2)=ur(2)-v1;

end


function Es=Eapp(tau)

    % ------------------------+
    %         1               |
    %          /\    1        |
    %    -----/--\--/----> t  |
    %        0    \/          |
    %              -1         |
    % ------------------------+

    Tobs_mod=mod(tau + 0.25,1);
    Es = 1 - abs(1- (Tobs_mod.*4 - 1));

end


function y=log_exp_plus(x)

    % y = log(1+exp(x))

    bdn=-6;
    y=(x<bdn).*log_exp_plus_expand(x)+(x>=bdn).*log(exp(x)+1);
    y(isnan(y))=x(isnan(y));
    
    function y=log_exp_plus_expand(x)
        y=exp(x)-exp(2.*x)./2+exp(3.*x)./3-exp(4.*x)./4+...
            exp(5.*x)./5-exp(6.*x)./6+exp(7.*x)./7-exp(8.*x)./8;
    end

end



function xi=sol_bnd(a,b,xi_n,xi_p)
    xi=fzero(@(x)bdeq(a,b,xi_n,xi_p,x),ap_sol(a,b,xi_n,xi_p));

    function res=bdeq(a,b,xi_n,xi_p,xi)
        res=a.*log_exp_plus(xi-xi_n)-b.*log_exp_plus(xi_p-xi)+xi;
    end
    
    function xi=ap_sol(a,b,xi_n,xi_p)
        xi=((max(0,xi_p)>xi_n).*a.*xi_n+(min(0,xi_n)<xi_p).*b.*xi_p)...
            ./(1+(max(0,xi_p)>xi_n).*a+(min(0,xi_n)<xi_p).*b);
    end

end

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
