function euler_post(cfg, m, sol)
% EULER_POST  Final plots: Mach, Cp, entropy error, vorticity, velocity,
% surface Cp distribution and convergence / force histories.

Lref = m.Lref;
xl = [m.xb - 3*Lref, m.xb + 12*Lref];
yl = [m.yb - 3*Lref, m.yb + 3*Lref];
win = find(m.xc >= xl(1) & m.xc <= xl(2) & m.yc >= yl(1) & m.yc <= yl(2));

cellplot(m,win,sol.Mach,'Mach Number',xl,yl,[]);
cellplot(m,win,sol.Cp,'Pressure Coefficient',xl,yl,[]);
cellplot(m,win,sol.dS,'Entropy Error  s/s_\infty - 1',xl,yl,[]);
vm = max(abs(sol.vort(win)));
cellplot(m,win,sol.vort,'Vorticity  \omega L_{ref}/U_\infty',xl,yl,[-0.3*vm 0.3*vm]);

% Velocity field on a uniform grid, masked inside the body
figure;
[Xq,Yq] = meshgrid(linspace(xl(1),xl(2),70),linspace(yl(1),yl(2),35));
uq = griddata(m.xc(win),m.yc(win),sol.W(win,2),Xq,Yq);
vq = griddata(m.xc(win),m.yc(win),sol.W(win,3),Xq,Yq);
in = inpolygon(Xq,Yq,m.bodyX,m.bodyY);
uq(in) = NaN;  vq(in) = NaN;
quiver(Xq,Yq,uq,vq);  hold on;
line(m.wx,m.wy,'Color','k','LineWidth',1.5);
axis equal;  xlim(xl);  ylim(yl);  xlabel('x');  ylabel('y');  title('Velocity Field');

% Surface pressure distribution (exact wall faces)
figure;
plot(sol.wallX,sol.wallCp,'.','MarkerSize',8);
set(gca,'YDir','reverse');  grid on;
xlabel('x');  ylabel('C_p');
title(sprintf('Surface C_p   (C_l = %.4f, C_d = %.4f)',sol.Cl,sol.Cd));

% Histories
figure;
if sol.transient
    xh = sol.tHist*sol.U_inf/Lref;   xlab = 't^* = t U_\infty / L_{ref}';
else
    xh = (1:sol.iter)';              xlab = 'Iteration';
end
subplot(3,1,1);  semilogy(xh,sol.resHist,'LineWidth',1.5);
grid on;  xlabel(xlab);  ylabel('Rel. density residual');  title('Convergence History');
subplot(3,1,2);  plot(xh,sol.ClHist,'LineWidth',1.5);  hold on;  plot(xh,sol.CdHist,'LineWidth',1.5);
grid on;  xlabel(xlab);  ylabel('Coefficient');  legend('C_l','C_d');
subplot(3,1,3);  plot(xh,sol.vpHist,'LineWidth',1.5);
grid on;  xlabel(xlab);  ylabel('v at wake probe [m/s]');

fprintf('\nAll plots generated.\n');
end

function cellplot(m,win,fld,ttl,xl,yl,clim)
figure;
patch('Faces',m.cells(win,:),'Vertices',[m.xn m.yn],'FaceVertexCData',fld(win), ...
      'FaceColor','flat','EdgeColor','none');
hold on;  line(m.wx,m.wy,'Color','k','LineWidth',1.5);
axis equal;  xlim(xl);  ylim(yl);  colorbar;
if ~isempty(clim), caxis(clim); end
xlabel('x');  ylabel('y');  title(ttl);
end
