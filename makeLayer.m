
classdef makeLayer < handle
    
    properties
        Temp = 300;
        layerLength = 1e-6;
        dopingDensity = 0;
        dopingType = 'n';
        materialType = 'mct';
        materialComp = 0.0 ;
        IDofLayer = 0;
        permittivity = 0;
        intrinsic = 0;
        eMobility = 0;
        hMobility = 0;
        bandGap = 0;
        affinity = 0;
        isWellDefined = false;
        startPos = 0;
        endPos = 0; 
    end
    
    methods
        function setTemperature(obj, temp_in)
            obj.Temp = temp_in;
        end
        
        function setLayerLength(obj, length_input)
            obj.layerLength = 100*length_input;
        end
        
        function setDopingDensity(obj, density_input, type)
            obj.dopingType = type;
            if obj.dopingType == 'n'
                obj.dopingDensity = density_input;
            else
                if obj.dopingType == 'p'
                obj.dopingDensity = -1 * density_input;
                else
                    fprintf('Unrecognized type of doping please try using "n" or "p"');
                end
            end
        end
        
        function setMaterialType(obj, mat_input)
            obj.materialType = mat_input;
        end
        
        function setMaterialComp(obj, mat_comp)
            obj.materialComp = mat_comp;
        end
        
        function setID(obj, ID)
            obj.IDofLayer = ID;
        end
        
        function deriveParam(obj)
            obj.permittivity = getPermittivity(obj.materialComp);
            obj.intrinsic = getIntrinsicConcentration(obj.materialComp, obj.Temp);
            obj.eMobility = getElectronMobility(obj.materialComp, obj.Temp);
            obj.hMobility = getHoleMobility(obj.materialComp, obj.Temp);
            obj.bandGap = getBandgap(obj.materialComp, obj.Temp);
            obj.affinity = getAffinity(obj.bandGap);
            obj.isWellDefined = true;
        end
        
     
    end
    
        

    
end







%%


function [ tmp ] = getBandgap(x , Temp) % (eV) Bandgap creator function for different x compositions and different temperatures
    
    tmp = -0.302 + 1.93 .* x + 5.35 .* 10^(-4) .* (1-2.*x) .* (-1882 + Temp.^3) ./ (255.2 + Temp.^2);
    
end

function [ tmp ] = getIntrinsicConcentration(x , Temp) % (cm-3) Intrinsic carrrier concentration for different x composition and different temperatures 

    bandGap = getBandgap(x, Temp);
    n_i_initial = 5.24256 - 3.5729 .* x - 4.74019 .* Temp .* 10^(-4) + 1.25942 .* 10^(-2) .* x .* Temp - 5.77046 .* x .* x - 4.24123 .* 10^(-6) .* Temp .* Temp;
    tmp = n_i_initial .* 10^(14) .* bandGap.^(0.75) .* Temp.^(1.5) .* exp( (-bandGap .* 1.60217662 * 10^(-19)) ./ (2 .* 1.38064852 * 10^(-23) .* Temp) );

end

function [ tmp ] = getAffinity(bandGap) % (eV) Affinity generator function

    tmp = 4.23 - 0.813 .* (bandGap - 0.083);
    
end

function [ tmp ] = getPermittivity(x) % (Farad/cm) Permittivity generator function

    vacuum_permittivity = 8.85 * 10^(-14);
    tmp = (20.5 - 15.6 .* x + 5.7 .* x .* x ) .* vacuum_permittivity;
    
end

function [ tmp ] = getElectronMobility(x , Temp) % (cm2.V-1.s-1) Estimated Electron mobility with given x mole fraction and temperature

    a=(0.2./x).^0.6;
    b=(0.2./x).^7.5;
    
    if Temp > 50
        tmp = (9 .* (10^8) .* b ) ./ ( Temp .^ (2.*a)) ;
    else
        tmp = (9 .* (10^8) .* b ) ./ ( ((1.18 .^ 10^5 ) ./ (2600 - abs((Temp-35).^2.07))) .^ (2.*a)) ;
    end
    
end

function [ tmp ] = getHoleMobility(x, Temp) % (cm2.V-1.s-1) eMobility/100 estimation

    tmp = getElectronMobility(x , Temp);
    tmp = tmp ./ 100;

end

function [ tmp ] = getBuiltinVoltage(x1 , x2 , N1 , N2 , Temp , junction_name) % (eV) For any hypothetical pn, np, nn junctions built-in voltages
    
    gap1 = getBandgap(x1,Temp); % (eV) Finding bandgap for x1
    gap2 = getBandgap(x2,Temp); % (eV) Finding bandgap for x2
    
    if junction_name == "pn"
        acceptor_doping = N1; % (cm-3) Defining acceptor doping as N1
        donor_doping = N2;    % (cm-3) Defining donor doping as N2
 
    else
        if juntion_name == "np" % may need to be corrected later on
        acceptor_doping = N2; % (cm-3) Defining acceptor doping as N2
        donor_doping = N1;    % (cm-3) Defining donor doping as N1
        end
    end
    
    int_conc1 = getIntrinsicConcentration(x1, Temp); % (cm-3) Find Intrinsic Concentration for x1
    int_conc2 = getIntrinsicConcentration(x2, Temp); % (cm-3) Find Intrinsic Concentration for x2
    
    aff1 = getAffinity(gap1); % (eV) Affinity of x1
    aff2 = getAffinity(gap2); % (eV) Affinity of x2
    delta_cond = aff2 - aff1 ; % (eV) Conduction band differences
    
    ktq = 1.38064852 .* 10^(-23) .* Temp ./ (1.6 .* 10^(-19)); % (eV) Thermal voltage for given Temperatures
     
    Nvaln1 = int_conc1 .* exp(gap1 ./ (2.*ktq)); % (cm-3) Valence band carrier concentration of x1
    Ncond2 = int_conc2 .* exp(gap2 ./ (2.*ktq)); % (cm-3) Conduction band carrier concentration of x2
    
    tmp = gap1 - delta_cond + ktq .* log( (acceptor_doping * donor_doping) ./ (Nvaln1 .* Ncond2) ); % Finding built-in voltage
  

end

function [ depletion_n, depletion_p, depletion_total ] = getDepletionWidth(x1, x2, N1, N2, Temp, junction_name, bias) % (um) Depletion layer widths

    Vbi = getBuiltinVoltage(x1, x2, N1, N2, Temp, junction_name);
    
    if junction_name == "pn"
        acceptor_doping = N1; % (cm-3) Defining acceptor doping as N1
        donor_doping = N2;    % (cm-3) Defining donor doping as N2
 
    else
        if juntion_name == "np" % may need to be corrected later on
        acceptor_doping = N2; % (cm-3) Defining acceptor doping as N2
        donor_doping = N1;    % (cm-3) Defining donor doping as N1
        end
    end
    
    perm1 = getPermittivity(x1);
    perm2 = getPermittivity(x2);
    
    depletion_p = sqrt( 2.*(perm1.*perm2.*donor_doping.*(Vbi-bias)) ./ (1.6e-19 .* acceptor_doping .* (perm1.*acceptor_doping + perm2.*donor_doping) ));
    depletion_n = sqrt( 2.*(perm1.*perm2.*acceptor_doping.*(Vbi-bias)) ./ (1.6e-19 .* donor_doping .* (perm1.*acceptor_doping + perm2.*donor_doping) ));
    depletion_p = depletion_p .* 10^4; % Converting to um
    depletion_n = depletion_n .* 10^4; % Converting to um
    
    depletion_total = depletion_p + depletion_n ;

end


