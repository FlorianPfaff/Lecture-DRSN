function visualizeCovEllipsoidsAndBounds(Cxx, Cyy, rmin, rmax, nExamples, nBounds)
% @author Florian Pfaff
% @date 2011-2021
arguments
    Cxx (1,1) double = 9
    Cyy (1,1) double = 4
    rmin (1,1) double {mustBeInRange(rmin,-1,1)} = -0.999 
    rmax (1,1) double {mustBeInRange(rmax,-1,1)} = 0.999
    nExamples (1,1) double {mustBeInteger, mustBePositive} = 40 
    nBounds (1,1) double {mustBeInteger, mustBePositive} = 20
end
addpath ellipsoidalToolboxEssentials

figure(1),clf,axis([-1.1*sqrt(Cxx),1.1*sqrt(Cxx),-1.1*sqrt(Cyy),1.1*sqrt(Cyy)]),hold on

plotExamples(Cxx, Cyy, rmin, rmax, nExamples, false)
axis([-10,10,-10,10])
drawBounds(Cxx, Cyy, rmax, nBounds, false)
end

function plotExamples(Cxx, Cyy, rmin, rmax, nExamples, savePlots)
arguments
    Cxx (1,1) double
    Cyy (1,1) double
    rmin (1,1) double {mustBeInRange(rmin,-1,1)}
    rmax (1,1) double {mustBeInRange(rmax,-1,1)}
    nExamples (1,1) double {mustBeInteger, mustBePositive}
    savePlots (1,1) logical = true
end
if rmin==-1.0 
    warning('Formula will not work for rmin=-1, increasing it slightly.')
    rmin = rmin + 0.001;
end
if rmax==1.0 
    warning('Formula will not work for rmax=1, reducing it slightly.')
    rmax = rmax - 0.001;
end
i = 0;
for r = linspace(rmin, rmax, nExamples)
    plot(ellipsoid([0;0],[Cxx,r*sqrt(Cxx*Cyy);r*sqrt(Cxx*Cyy),Cyy]),'r');
    drawnow
    pause(0.1)
    if savePlots
        set(gcf,'color',[1,1,1]);
        exportgraphics(gcf,sprintf('example%02d.png',i));
    end
    i = i+1;
end
end

function drawBounds(Cxx, Cyy, rabsmax, nBounds, savePlots)
arguments
    Cxx (1,1) double
    Cyy (1,1) double
    rabsmax (1,1) double
    nBounds (1,1) double {mustBeInteger, mustBePositive} = 20
    savePlots (1,1) logical = true
end
i = 0;
eta=@(kappa)(1-sqrt(rabsmax^2+kappa^2*(1-rabsmax^2)^2))/(1-rabsmax^2);
for kappa = linspace(-0.49, 0.49, nBounds)
    Ctotal=[(1/(eta(kappa)-kappa))*Cxx,0;0,(1/(eta(kappa)+kappa))*Cyy];
    plot(ellipsoid([0;0],Ctotal),'b')
    drawnow
    pause(0.1)
    if savePlots
        set(gcf,'color',[1,1,1]);
        exportgraphics(gcf,sprintf('bounds%02d.png',i));
    end
    i = i+1;
end
end