delay=0.0001;

figure
filename = 'dynamic_iv.gif';
% filename = 'MonteCarlo_20.gif';

% figure
% for i=t_lst2:10:length(t)-10
% plot([vdsf(t(i)),vdsf(t(i+10))],[jtot(i,1),jtot(i+10,1)],'b-')
% hold on

% outw=out.*(pa.ch.*pa.wl).^-1;

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




      drawnow
      frame = getframe(1);
      im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);
      if i == 2
          imwrite(imind,cm,filename,'gif', 'DelayTime', delay, 'Loopcount',inf);
      else
          imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', delay);
      end

end
