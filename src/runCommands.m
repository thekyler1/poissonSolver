clear all;

mySim = numSimPoisson;
mySim.setSimParameters(-0.1, 90);
mySim.newLayer;
mySim.layerProp(3e-6,1e16,'n','mct',0.3);
mySim.endLayer;

% mySim.newLayer;
% mySim.layerProp(0.5e-6,1e18,'n','mct',0.35);
% mySim.endLayer; 

mySim.newLayer;
mySim.layerProp(1e-6,2e17,'p','mct',0.5);
mySim.endLayer;

mySim.deriveScaling;
mySim.initConditions;


%%

mySim.nonBiasSim(1000);
mySim.biasClosing(1000);
mySim.chargesDerivative;
mySim.findTotCurrent;

mySim.totCurrent

