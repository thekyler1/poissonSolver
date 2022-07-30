
k = linspace(-0.1,0.3,100   );

for j = 1:100

    
    aa = simCat;
    aa.setSimParameters(k(j), 90);

    aa.newLayer;
    aa.layerProp(3e-6,1e16,'n','mct',0.3);
    aa.endLayer;

    % aa.newLayer;
    % aa.layerProp(0.2e-6,1e18,'n','mct',0.35);
    % aa.endLayer;

    aa.newLayer;
    aa.layerProp(1e-6,2e17,'p','mct',0.5);
    aa.endLayer;

    aa.simHandle;

    aa.deriveScaling;
    aa.initConditions;

    aa.startSim(1000);

    aa.sim1.chargesDerivative;
    aa.sim1.findTotCurrent;

    x(j) = aa.sim1.totCurrent

end












