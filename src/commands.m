bb = simCat;
bb.setSimParameters(-0.1, 90);

bb.newLayer;
bb.layerProp(3e-6,1e16,'n','mct',0.3);
bb.endLayer;

% bb.newLayer;
% bb.layerProp(0.5e-6,1e18,'n','mct',0.35);
% bb.endLayer;

bb.newLayer;
bb.layerProp(1e-6,2e17,'p','mct',0.5);
bb.endLayer;

bb.simHandle;

bb.deriveScaling;
bb.initConditions;

bb.startSim(1000);

%%

bb.sim1.chargesDerivative;
bb.sim1.findTotCurrent;

bb.sim1.totCurrent

bb.sim1.currentDens;

plot(bb.sim1.dummyVar)
hold on
plot(bb.sim1.totderCurrDensBiased)



% plot(aa.sim1.np_intrinsic_diff)
