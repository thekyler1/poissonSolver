
k = linspace(-0.2,0.1,26);
len = length(k);
out = zeros(1,len);

for j = 1:len   

    mySim = simCat;
    mySim.setSimParameters(k(j), 90);

    mySim.newLayer;
    mySim.layerProp(3e-6,1e16,'n','mct',0.3);
    mySim.endLayer;

    % mySim.newLayer;
    % mySim.layerProp(0.2e-6,1e18,'n','mct',0.35);
    % mySim.endLayer;

    mySim.newLayer;
    mySim.layerProp(1e-6,2e17,'p','mct',0.5);
    mySim.endLayer;

    mySim.simHandle;

    mySim.deriveScaling;
    mySim.initConditions;

    mySim.startSim(1000);

    mySim.sim1.chargesDerivative;
    mySim.sim1.findTotCurrent;

    out(1,j) = mySim.sim1.totCurrent;

end












