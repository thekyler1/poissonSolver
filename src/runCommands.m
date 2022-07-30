clear all;

aa = numSimPoisson;
aa.setSimParameters(-0.1, 90);
aa.newLayer;
aa.layerProp(3e-6,1e16,'n','mct',0.3);
aa.endLayer;

% aa.newLayer;
% aa.layerProp(0.5e-6,1e18,'n','mct',0.35);
% aa.endLayer; 

aa.newLayer;
aa.layerProp(1e-6,2e17,'p','mct',0.5);
aa.endLayer;

aa.deriveScaling;
aa.initConditions;


%%

aa.nonBiasSim(1000);
aa.biasClosing(1000);
aa.chargesDerivative;
aa.findTotCurrent;

aa.totCurrent

