function PlotSTAR(ShatClean,RClean,thetaTarget,phiTarget,sty)

xSTAR = RClean.*sind(thetaTarget).*cosd(phiTarget); ySTAR = RClean.*sind(thetaTarget).*sind(phiTarget); zSTAR = RClean.*cosd(thetaTarget);

figure(1)
subplot(1,3,1)
hold on
plot3(xSTAR, ySTAR, 20*log10(ShatClean), sty,'MarkerSize',20,'LineWidth',2)
xlabel('X')
zlabel('Shat (dB)')
set(gca,'FontSize',20)
view(0,0)

subplot(1,3,2)
hold on
plot3(xSTAR, ySTAR, 20*log10(ShatClean), sty,'MarkerSize',20,'LineWidth',2)
ylabel('Y')
zlabel('Shat (dB)')
set(gca,'FontSize',20)
view(90,0)

subplot(1,3,3)
hold on
plot3(xSTAR, zSTAR, 20*log10(ShatClean), sty,'MarkerSize',20,'LineWidth',2)
ylabel('Z')
zlabel('Shat (dB)')
set(gca,'FontSize',20)
view(90,0)
% w = waitforbuttonpress;

end

