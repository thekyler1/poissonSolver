

classdef numSimPoisson < handle
    
    properties (Constant)
        electron_charge = 1.6e-19; % Coulomb - charge of electron
        boltzmann_const = 1.38064852e-23; % m2.kg.s-2.K-1 - Boltzmann constant
        vacuum_permittivity = 8.85 * 10^(-14); % Farad/cm - Vacuum permittivity
        electron_mass = 9.10938e-31; % kg - mass of electron
        mct_electron_mass = 0.006937 * 9.10938e-31; % kg - electron effective mass
        mct_hole_mass = 0.55 * 9.10938e-31; % kg - hole effective mass
    end
    
    properties % Layer Definition and Environment
        V_bias = 0.0; % Default biasing value
        Temp = 0; % 0K defined for user to submit
        numofLayers = 0; % Initial number of layers    
        layerID = 0;
        layerNames;
        lastIteration;
        
        layer1;
        layer2;
        layer3;
        
        sim1;
        sim2;
        sim3;
        
        
    end
    
    properties % Device Properties
        device_length = 0;
        ktq = 0;
        concentration_scaling = 0;
        permittivity_scaling = 0;
        debye_length = 0;
        mobility_scaling = 0;
        currentDens_scaling = 0;
        length_scaling = 0;
        time_scaling = 0;
        isScaled = false;
        total_mesh_points = 0;
        dopDensProf;
        permittivity;
        intrinsic;
        voltage;
        voltageBiased;
        eMobility;
        hMobility;
        bandGap; % eV
        
    end
    
    properties % Matrix Definitions
        A_matrix;
        B_matrix;
        Bsub_matrix;
        X_matrix;
        S_matrix_upper;
        U_matrix_upper;
        Y_matrix_upper;
        lambda_upper;
        S_matrix_lower;
        U_matrix_lower;
        Y_matrix_lower;
        lambda_lower;
        lambda;
        error;
    end
    
    properties % Simulation Parameters    
        nCons;
        pCons;
        charges;
        secondDerVolt;
        totCurrent;
        
        nConsBiased;
        pConsBiased;
        chargesBiased;
        firstDerVoltBiased;
        secondDerVoltBiased;
        
        v_diff;
        n_diff;
        p_diff;
        v_diffBiased;
        n_diffBiased;
        p_diffBiased;
        
        np_intrinsic_diff;
        np_intrinsic_diffBiased;
        NcurrentDens;
        PcurrentDens;
        derNcurrentDens;
        derPcurrentDens;
        NcurrentDensBiased;
        PcurrentDensBiased;
        derNcurrentDensBiased;
        derPcurrentDensBiased;
        totderCurrDensBiased;
        
        paramTotCurrent;

        dummyVar;
        
        dummyVar1;
        dummyVar2;
        dummyVar3;
        
        tau_radiative;
        R_radiative;
        R_srh;
        
        
    end
    
    methods
        
        function setTemperature(obj, temp_in)
            if temp_in < 90
                fprintf('Enter temperatures higher than 90K!');
            end
            obj.Temp = temp_in;
        end 
        
        function setBias(obj, bias_in)
            obj.V_bias = bias_in;
        end
        
        function setSimParameters(obj, bias_in, temp_in)
            obj.V_bias = bias_in;
            obj.Temp = temp_in;
        end
        
        function y = getBiasingVoltage(obj)
            y = obj.V_bias;
        end
        
        function newLayer(obj)
            if obj.numofLayers == 0
                obj.layer1 = makeLayer;
                obj.layer1.setID(obj.numofLayers+1);
                obj.layerNames = {'layer1'};
                obj.layer1.setTemperature(obj.Temp);
            end
            if obj.numofLayers == 1
                obj.layer2 = makeLayer;
                obj.layer2.setID(obj.numofLayers+1);
                obj.layerNames = {'layer1' 'layer2'};
                obj.layer2.setTemperature(obj.Temp);
            end
            if obj.numofLayers == 2
                obj.layer3 = makeLayer;
                obj.layer3.setID(obj.numofLayers+1);
                obj.layerNames = {'layer1' 'layer2' 'layer3'};
                obj.layer3.setTemperature(obj.Temp);
            end
            obj.numofLayers = obj.numofLayers+1;
            obj.layerID = obj.numofLayers;
        end
        
        function layerProp(obj, layerLength, dopingDensity, dopingType, materialType, materialComp)
            if obj.numofLayers == 0
                disp('There is no layer created')
            end
            if obj.layerID == 1
                obj.layer1.setLayerLength(layerLength);
                obj.layer1.setDopingDensity(dopingDensity, dopingType);
                obj.layer1.setMaterialType(materialType);
                obj.layer1.setMaterialComp(materialComp);
            end
            if obj.layerID == 2                
                obj.layer2.setLayerLength(layerLength);
                obj.layer2.setDopingDensity(dopingDensity, dopingType);
                obj.layer2.setMaterialType(materialType);
                obj.layer2.setMaterialComp(materialComp);
            end
            if obj.layerID == 3                
                obj.layer3.setLayerLength(layerLength);
                obj.layer3.setDopingDensity(dopingDensity, dopingType);
                obj.layer3.setMaterialType(materialType);
                obj.layer3.setMaterialComp(materialComp);
            end
        end
        
        function endLayer(obj)
            if obj.numofLayers == 0
                fprintf('No layers defined\n')
            end
            if obj.Temp == 0
                fprintf('Set temperature first!\n')
            else
                obj.endLayerProp;
            end
        end
        
        function endLayerProp(obj)
            if obj.layerID == 1
                obj.layer1.setTemperature(obj.Temp);
                obj.layer1.deriveParam;
            end
            if obj.layerID == 2
                obj.layer2.setTemperature(obj.Temp);
                obj.layer2.deriveParam;
            end
            if obj.layerID == 3
                obj.layer3.setTemperature(obj.Temp);
                obj.layer3.deriveParam;
            end
        end
        
        function deriveScaling(obj)
            if obj.Temp == 0
                fprintf('You need to define temperature first!\n');
                return
            else
                obj.ktq = obj.boltzmann_const * obj.Temp / obj.electron_charge;
            end
            
            if isa(obj.layer1, 'makeLayer')
                if obj.layer1.isWellDefined
                    ly1_dd = obj.layer1.dopingDensity;
                    ly1_perm = obj.layer1.permittivity;
                    ly1_int = obj.layer1.intrinsic;
                    ly1_eMob = obj.layer1.eMobility;
                    ly1_hMob = obj.layer1.hMobility;   
                else
                    fprintf('You need to complete the layer1 first! call endLayer\n');
                    return
                end
            else 
                fprintf('You need to define at least two layers\n');
                return
            end
            
            if isa(obj.layer2, 'makeLayer')
                if obj.layer2.isWellDefined
                    ly2_dd = obj.layer2.dopingDensity;
                    ly2_perm = obj.layer2.permittivity;
                    ly2_int = obj.layer2.intrinsic;
                    ly2_eMob = obj.layer2.eMobility;
                    ly2_hMob = obj.layer2.hMobility;
                else
                    fprintf('You need to complete the layer2 before calling this function! call endLayer');
                    return
                end
            else 
                fprintf('You need to define at least two layers');
                return
            end
            
            if isa(obj.layer3, 'makeLayer')
                if obj.layer3.isWellDefined
                    ly3_dd = obj.layer3.dopingDensity;
                    ly3_perm = obj.layer3.permittivity;
                    ly3_int = obj.layer3.intrinsic;
                    ly3_eMob = obj.layer3.eMobility;
                    ly3_hMob = obj.layer3.hMobility;
                else
                    fprintf('Complete layer3 before calling this function! call endLayer');
                    return
                end
            else
                ly3_dd = 0;
                ly3_perm = 0;
                ly3_int = 0;
                ly3_eMob = 0;
                ly3_hMob = 0;
            end
            
            obj.concentration_scaling = max([abs(ly1_dd) abs(ly2_dd) abs(ly3_dd) ly1_int ly2_int ly3_int]);
            obj.permittivity_scaling = max([ly1_perm ly2_perm ly3_perm]);
            obj.mobility_scaling = max([ly1_eMob ly2_eMob ly3_eMob ly1_hMob ly2_hMob ly3_hMob]);
            obj.debye_length = sqrt(obj.permittivity_scaling * obj.ktq / (obj.electron_charge * obj.concentration_scaling));
            obj.length_scaling = obj.debye_length;
            obj.currentDens_scaling = obj.electron_charge * obj.mobility_scaling * obj.concentration_scaling * obj.ktq / obj.length_scaling;
            obj.time_scaling = obj.length_scaling .* obj.length_scaling ./ obj.mobility_scaling ./ obj.ktq;
            obj.isScaled = true;
            
            
        end
        
        function setLayerPos(obj)

            if isa(obj.layer3, 'makeLayer') && obj.layer3.isWellDefined
                obj.device_length = obj.layer1.layerLength + obj.layer2.layerLength + obj.layer3.layerLength;
                obj.total_mesh_points = obj.device_length / obj.length_scaling;
                obj.total_mesh_points = int16(obj.total_mesh_points);
                obj.layer1.startPos = 0;
                obj.layer1.endPos = obj.layer1.layerLength / obj.device_length;
                obj.layer2.startPos = obj.layer1.endPos;
                obj.layer2.endPos = obj.layer2.startPos + obj.layer2.layerLength / obj.device_length;
                obj.layer3.startPos = obj.layer2.endPos;
                obj.layer3.endPos = obj.layer3.startPos + obj.layer3.layerLength / obj.device_length;
                
                obj.layer1.startPos = int16(1+obj.layer1.startPos * obj.total_mesh_points);
                obj.layer1.endPos = int16(1+obj.layer1.endPos * obj.total_mesh_points-1);
                obj.layer2.startPos = int16(1+obj.layer2.startPos * obj.total_mesh_points);
                obj.layer2.endPos = int16(1+obj.layer2.endPos * obj.total_mesh_points-1);
                obj.layer3.startPos = int16(1+obj.layer3.startPos * obj.total_mesh_points);
                obj.layer3.endPos = int16(1+obj.layer3.endPos * obj.total_mesh_points-1);
                
            else
                obj.device_length = obj.layer1.layerLength + obj.layer2.layerLength;
                obj.total_mesh_points = obj.device_length / obj.length_scaling;
                obj.total_mesh_points = int16(obj.total_mesh_points);
                obj.layer1.startPos = 0;
                obj.layer1.endPos = obj.layer1.layerLength / obj.device_length;
                obj.layer2.startPos = obj.layer1.endPos;
                obj.layer2.endPos = obj.layer2.startPos + obj.layer2.layerLength / obj.device_length;
                
                obj.layer1.startPos = int16(1+obj.layer1.startPos * obj.total_mesh_points);
                obj.layer1.endPos = int16(1+obj.layer1.endPos * obj.total_mesh_points-1);
                obj.layer2.startPos = int16(1+obj.layer2.startPos * obj.total_mesh_points);
                obj.layer2.endPos = int16(1+obj.layer2.endPos * obj.total_mesh_points-1);
            end    
        end
        
        function setZeroMatrices(obj)
            obj.dopDensProf = zeros(1, obj.total_mesh_points+2);
            obj.permittivity = zeros(1, obj.total_mesh_points+2);
            obj.intrinsic = zeros(1, obj.total_mesh_points+2);
            obj.voltage = zeros(1, obj.total_mesh_points+2);
            obj.eMobility = zeros(1, obj.total_mesh_points+2);
            obj.hMobility = zeros(1, obj.total_mesh_points+2);
            
            obj.A_matrix = zeros(1, obj.total_mesh_points+2);
            obj.B_matrix = zeros(1, obj.total_mesh_points+2);
            obj.Bsub_matrix = zeros(1, obj.total_mesh_points+2);
            obj.X_matrix = zeros(1, obj.total_mesh_points+2);
            
            obj.S_matrix_upper = zeros(1, obj.total_mesh_points+2);
            obj.U_matrix_upper = zeros(1, obj.total_mesh_points+2);
            obj.Y_matrix_upper = zeros(1, obj.total_mesh_points+2);
            obj.lambda_upper = zeros(1, obj.total_mesh_points+2);
            
            obj.S_matrix_lower = zeros(1, obj.total_mesh_points+2);
            obj.U_matrix_lower = zeros(1, obj.total_mesh_points+2);
            obj.Y_matrix_lower = zeros(1, obj.total_mesh_points+2);
            obj.lambda_lower = zeros(1, obj.total_mesh_points+2);
            
            obj.lambda = zeros(1, obj.total_mesh_points+2);
            
        end
        
        function initConditions(obj)
            
            obj.setLayerPos;
            obj.setZeroMatrices;
            
            
            if isa(obj.layer3, 'makeLayer') && obj.layer3.isWellDefined
                obj.dopDensProf(1,1) = obj.layer1.dopingDensity / obj.concentration_scaling;
                obj.dopDensProf(1,1+obj.layer1.startPos:1+obj.layer1.endPos) = obj.layer1.dopingDensity / obj.concentration_scaling;
                obj.dopDensProf(1,1+obj.layer2.startPos:1+obj.layer2.endPos) = obj.layer2.dopingDensity / obj.concentration_scaling;
                obj.dopDensProf(1,1+obj.layer3.startPos:2+obj.layer3.endPos) = obj.layer3.dopingDensity / obj.concentration_scaling;
                
                obj.permittivity(1,1) = obj.layer1.permittivity / obj.permittivity_scaling;
                obj.permittivity(1,1+obj.layer1.startPos:1+obj.layer1.endPos) = obj.layer1.permittivity / obj.permittivity_scaling;
                obj.permittivity(1,1+obj.layer2.startPos:1+obj.layer2.endPos) = obj.layer2.permittivity / obj.permittivity_scaling;
                obj.permittivity(1,1+obj.layer3.startPos:2+obj.layer3.endPos) = obj.layer3.permittivity / obj.permittivity_scaling;
                                
                obj.intrinsic(1,1) = obj.layer1.intrinsic / obj.concentration_scaling;
                obj.intrinsic(1,1+obj.layer1.startPos:1+obj.layer1.endPos) = obj.layer1.intrinsic / obj.concentration_scaling;
                obj.intrinsic(1,1+obj.layer2.startPos:1+obj.layer2.endPos) = obj.layer2.intrinsic / obj.concentration_scaling;
                obj.intrinsic(1,1+obj.layer3.startPos:2+obj.layer3.endPos) = obj.layer3.intrinsic / obj.concentration_scaling;
                                
                obj.eMobility(1,1) = obj.layer1.eMobility / obj.mobility_scaling;
                obj.eMobility(1,1+obj.layer1.startPos:1+obj.layer1.endPos) = obj.layer1.eMobility / obj.mobility_scaling;
                obj.eMobility(1,1+obj.layer2.startPos:1+obj.layer2.endPos) = obj.layer2.eMobility / obj.mobility_scaling;
                obj.eMobility(1,1+obj.layer3.startPos:2+obj.layer3.endPos) = obj.layer3.eMobility / obj.mobility_scaling;
                                
                obj.hMobility(1,1) = obj.layer1.hMobility / obj.mobility_scaling;
                obj.hMobility(1,1+obj.layer1.startPos:1+obj.layer1.endPos) = obj.layer1.hMobility / obj.mobility_scaling;
                obj.hMobility(1,1+obj.layer2.startPos:1+obj.layer2.endPos) = obj.layer2.hMobility / obj.mobility_scaling;
                obj.hMobility(1,1+obj.layer3.startPos:2+obj.layer3.endPos) = obj.layer3.hMobility / obj.mobility_scaling;
                
                obj.bandGap(1,1) = obj.layer1.bandGap;
                obj.bandGap(1,1+obj.layer1.startPos:1+obj.layer1.endPos) = obj.layer1.bandGap;
                obj.bandGap(1,1+obj.layer2.startPos:1+obj.layer2.endPos) = obj.layer2.bandGap;
                obj.bandGap(1,1+obj.layer3.startPos:2+obj.layer3.endPos) = obj.layer3.bandGap;
                
                
                if obj.layer1.dopingType == 'n'
                    obj.voltage(1,1) = log(obj.layer1.dopingDensity/obj.layer1.intrinsic) - 0.5*obj.V_bias/obj.ktq;
                    obj.voltage(1,1+obj.layer1.startPos:1+obj.layer1.endPos) = log(obj.layer1.dopingDensity/obj.layer1.intrinsic) - 0.5*obj.V_bias/obj.ktq;
                else
                    if obj.layer1.dopingType == 'p'
                        obj.voltage(1,1) = -1*log(-1*obj.layer1.dopingDensity/obj.layer1.intrinsic) + 0.5*obj.V_bias/obj.ktq;
                        obj.voltage(1,1+obj.layer1.startPos:1+obj.layer1.endPos) = -1*log(-1*obj.layer1.dopingDensity/obj.layer1.intrinsic) + 0.5*obj.V_bias/obj.ktq;
                    end
                end
                %-------
                if obj.layer2.dopingType == 'n'
                    obj.voltage(1,1+obj.layer2.startPos:1+obj.layer2.endPos) = log(obj.layer2.dopingDensity/obj.layer2.intrinsic) - 0.5*obj.V_bias/obj.ktq;
                else
                    if obj.layer2.dopingType == 'p'
                        obj.voltage(1,1+obj.layer2.startPos:1+obj.layer2.endPos) = -1*log(-1*obj.layer2.dopingDensity/obj.layer2.intrinsic) + 0.5*obj.V_bias/obj.ktq;
                    end
                end
                %----------
                if obj.layer3.dopingType == 'n'
                    obj.voltage(1,1+obj.layer3.startPos:2+obj.layer3.endPos) = log(obj.layer3.dopingDensity/obj.layer3.intrinsic) - 0.5*obj.V_bias/obj.ktq;
                else
                    if obj.layer3.dopingType == 'p'
                        obj.voltage(1,1+obj.layer3.startPos:2+obj.layer3.endPos) = -1*log(-1*obj.layer3.dopingDensity/obj.layer3.intrinsic) + 0.5*obj.V_bias/obj.ktq;
                    end
                end

            else
                obj.dopDensProf(1,1) = obj.layer1.dopingDensity / obj.concentration_scaling;
                obj.dopDensProf(1,1+obj.layer1.startPos:1+obj.layer1.endPos) = obj.layer1.dopingDensity / obj.concentration_scaling;
                obj.dopDensProf(1,1+obj.layer2.startPos:2+obj.layer2.endPos) = obj.layer2.dopingDensity / obj.concentration_scaling;
                                
                obj.permittivity(1,1) = obj.layer1.permittivity / obj.permittivity_scaling;
                obj.permittivity(1,1+obj.layer1.startPos:1+obj.layer1.endPos) = obj.layer1.permittivity / obj.permittivity_scaling;
                obj.permittivity(1,1+obj.layer2.startPos:2+obj.layer2.endPos) = obj.layer2.permittivity / obj.permittivity_scaling;
                
                obj.intrinsic(1,1) = obj.layer1.intrinsic / obj.concentration_scaling;
                obj.intrinsic(1,1+obj.layer1.startPos:1+obj.layer1.endPos) = obj.layer1.intrinsic / obj.concentration_scaling;
                obj.intrinsic(1,1+obj.layer2.startPos:2+obj.layer2.endPos) = obj.layer2.intrinsic / obj.concentration_scaling;
                
                obj.eMobility(1,1) = obj.layer1.eMobility / obj.mobility_scaling;
                obj.eMobility(1,1+obj.layer1.startPos:1+obj.layer1.endPos) = obj.layer1.eMobility / obj.mobility_scaling;
                obj.eMobility(1,1+obj.layer2.startPos:2+obj.layer2.endPos) = obj.layer2.eMobility / obj.mobility_scaling;
                                
                obj.hMobility(1,1) = obj.layer1.hMobility / obj.mobility_scaling;
                obj.hMobility(1,1+obj.layer1.startPos:1+obj.layer1.endPos) = obj.layer1.hMobility / obj.mobility_scaling;
                obj.hMobility(1,1+obj.layer2.startPos:2+obj.layer2.endPos) = obj.layer2.hMobility / obj.mobility_scaling;
                
                obj.bandGap(1,1) = obj.layer1.bandGap;
                obj.bandGap(1,1+obj.layer1.startPos:1+obj.layer1.endPos) = obj.layer1.bandGap;
                obj.bandGap(1,1+obj.layer2.startPos:2+obj.layer2.endPos) = obj.layer2.bandGap;
                
                
                if obj.layer1.dopingType == 'n'
                    obj.voltage(1,1) = log(obj.layer1.dopingDensity/obj.layer1.intrinsic) - 0.5*obj.V_bias/obj.ktq;
                    obj.voltage(1,1+obj.layer1.startPos:1+obj.layer1.endPos) = log(obj.layer1.dopingDensity/obj.layer1.intrinsic) - 0.5*obj.V_bias/obj.ktq;
                else
                    if obj.layer1.dopingType == 'p'
                        obj.voltage(1,1) = -1*log(obj.layer1.dopingDensity/obj.layer1.intrinsic) + 0.5*obj.V_bias/obj.ktq;
                        obj.voltage(1,1+obj.layer1.startPos:1+obj.layer1.endPos) = -1*log(-1*obj.layer1.dopingDensity/obj.layer1.intrinsic) + 0.5*obj.V_bias/obj.ktq;
                    end
                end
                %--------------
                if obj.layer2.dopingType == 'n'
                    obj.voltage(1,1+obj.layer2.startPos:2+obj.layer2.endPos) = log(obj.layer2.dopingDensity/obj.layer2.intrinsic) - 0.5*obj.V_bias/obj.ktq;
                else
                    if obj.layer2.dopingType == 'p'
                        obj.voltage(1,1+obj.layer2.startPos:2+obj.layer2.endPos) = -1*log(-1*obj.layer2.dopingDensity/obj.layer2.intrinsic) + 0.5*obj.V_bias/obj.ktq;
                    end
                end
            end 
        end
        
        function startSim(obj, iteration)
            obj.lastIteration = iteration;
            kk = 1;
            err_counter = 0;
            
            
                while (kk < iteration+1)
                    
            obj.A_matrix = (exp(obj.voltage).*exp(0.5*obj.V_bias/obj.ktq) .* obj.intrinsic + exp(-1*obj.voltage) .* exp(0.5*obj.V_bias/obj.ktq) .* obj.intrinsic) .* obj.length_scaling .* obj.length_scaling ./ (obj.debye_length .* obj.debye_length);
            obj.Bsub_matrix = (obj.dopDensProf - exp(obj.voltage).* exp(0.5*obj.V_bias/obj.ktq) .* obj.intrinsic + exp(-1*obj.voltage).* exp(0.5*obj.V_bias/obj.ktq) .* obj.intrinsic) .* obj.length_scaling .* obj.length_scaling ./ (obj.debye_length .* obj.debye_length);  
            obj.B_matrix = obj.Bsub_matrix;  
            obj.B_matrix(1,2:obj.total_mesh_points+1) = obj.B_matrix(1,2:obj.total_mesh_points+1) + (obj.permittivity(1,3:obj.total_mesh_points+2) - obj.permittivity(1,2:obj.total_mesh_points+1)) .* (obj.voltage(1,3:obj.total_mesh_points+2) - obj.voltage(1,2:obj.total_mesh_points+1)) + obj.permittivity(1,2:obj.total_mesh_points+1) .* (obj.voltage(1,1:obj.total_mesh_points) + obj.voltage(1,3:obj.total_mesh_points+2) - 2*obj.voltage(1,2:obj.total_mesh_points+1));

            obj.X_matrix(1,1:obj.total_mesh_points+1) = obj.A_matrix(1,1:obj.total_mesh_points+1) + obj.permittivity(1,1:obj.total_mesh_points+1) + obj.permittivity(1,2:obj.total_mesh_points+2);
            obj.X_matrix(1,obj.total_mesh_points+2) = obj.A_matrix(1,obj.total_mesh_points+2)+ 2*obj.permittivity(1,obj.total_mesh_points+2);

            obj.S_matrix_upper = obj.X_matrix;
            obj.U_matrix_upper = obj.B_matrix;


            for j = obj.total_mesh_points+1:-1:1
                obj.S_matrix_upper(1,j) = obj.X_matrix(1,j) - obj.permittivity(1,j)*obj.permittivity(1,j) / obj.S_matrix_upper(1,j+1);
            end

            for j = obj.total_mesh_points+1:-1:1
                obj.U_matrix_upper(1,j) = obj.B_matrix(1,j) + obj.permittivity(1,j) * obj.U_matrix_upper(1,j+1) / obj.S_matrix_upper(1,j+1);
            end

            obj.Y_matrix_upper(1,1) = obj.U_matrix_upper(1,1);

            for j = 2:obj.total_mesh_points+2
                obj.Y_matrix_upper(1,j) = obj.U_matrix_upper(1,j) + obj.permittivity(1,j) * obj.U_matrix_upper(1,j-1) / obj.S_matrix_upper(1,j-1);
            end

            for j = 2:obj.total_mesh_points+2
                obj.lambda_upper(1,j) = obj.Y_matrix_upper(1,j) / obj.S_matrix_upper(1,j);
            end

            obj.S_matrix_lower(1,1) = obj.X_matrix(1,1); % Initial condition for S_matrix
            obj.U_matrix_lower(1,1) = obj.B_matrix(1,1); % Initial condition for U_matrix


            for j = 2:obj.total_mesh_points+2
                obj.S_matrix_lower(1,j) = obj.X_matrix(1,j) - obj.permittivity(1,j) * obj.permittivity(1,j) / obj.S_matrix_lower(1,j-1);
            end

            for j = 2:obj.total_mesh_points+2
                obj.U_matrix_lower(1,j) = obj.B_matrix(1,j) + obj.permittivity(1,j) * obj.U_matrix_lower(1,j-1) / obj.S_matrix_lower(1,j-1);
            end

            obj.Y_matrix_lower(1, obj.total_mesh_points+2) = obj.U_matrix_lower(1, obj.total_mesh_points+2);

            for j = obj.total_mesh_points+1:-1:1
                obj.Y_matrix_lower(1,j) = obj.U_matrix_lower(1,j) + obj.permittivity(1,j) * obj.U_matrix_lower(1,j+1) / obj.S_matrix_lower(1,j+1);
            end

            for j = 2:obj.total_mesh_points+2
                obj.lambda_lower(1,j) = obj.Y_matrix_lower(1,j) / obj.S_matrix_lower(1,j);
            end

            obj.lambda = (obj.lambda_lower + obj.lambda_upper) / 2;

            
            
            error_prev = obj.error;
            
            damping = 1;
            obj.voltage = obj.voltage + damping * obj.lambda;
            obj.error = norm(obj.lambda_upper,2) / sqrt(double(obj.total_mesh_points+2));


            if error_prev < obj.error
                err_counter = err_counter + 1;
                if err_counter == 2
                    fprintf("Simulation converged with error = %e", obj.error);
                    fprintf(" at iteration #%d\n", kk);
                    break
                end
            end
            
            
            kk = kk + 1;
            
                end
                
                if kk == iteration + 1 
                    kk = kk - 1;
                    fprintf("Simulation ended with error = %e", obj.error);
                    fprintf(" at iteration #%d\n", kk);
                end
        end
        
        function nonBiasSim(obj, iteration)
            obj.lastIteration = iteration;
            tmp = obj.V_bias;
            tmpVolt = obj.voltage;
            obj.V_bias = 0.0;
            obj.initConditions;
            obj.startSim(obj.lastIteration);
            obj.voltageBiased = tmpVolt;
            obj.V_bias = tmp;
        end
        
        function zeroBiasSim(obj, iteration)
            obj.lastIteration = iteration;
            tmp_bias = obj.V_bias;
            obj.V_bias = 0.0;
            obj.initConditions;
            obj.startSim(obj.lastIteration);
            obj.V_bias = tmp_bias;            
        end
        
        function biasSim(obj, iteration)
            obj.lastIteration = iteration;
            tmp_voltage = obj.voltage;
            obj.voltage(1,1) = log(obj.layer1.dopingDensity/obj.layer1.intrinsic) - 0.5*obj.V_bias/obj.ktq;
            obj.voltage(1,1+obj.layer2.endPos) = -1*log(-1*obj.layer2.dopingDensity/obj.layer2.intrinsic) + 0.5*obj.V_bias/obj.ktq;
            obj.startSim(obj.lastIteration);
            obj.voltageBiased = obj.voltage;
            obj.voltage = tmp_voltage;

            
                        
        end
        
        function biasClosing(obj, iteration)
            obj.lastIteration = iteration;
            tmp_zeroVoltage = obj.voltage;
            last_bias = obj.V_bias;
            
            start = -0.01;
            step = -0.01;
            
            if obj.V_bias > 0
                start = 0.01;
                step = 0.01;
            end
            
            for v = start:step:last_bias

%                 lineLength = strlength("Biasing voltage: -2.000000e-01,     Simulation converged with error = 1.615055e-14 at iteration #202");
%                 fprintf(repmat('\b',1,lineLength));
                obj.V_bias = v;
                fprintf("Biasing voltage: %e,    ", v);
                obj.voltage(1,1) = log(obj.layer1.dopingDensity/obj.layer1.intrinsic) - 0.5*v/obj.ktq;
            if isa(obj.layer3, 'makeLayer') && obj.layer3.isWellDefined
                obj.voltage(1,2+obj.layer3.endPos) = -1*log(-1*obj.layer3.dopingDensity/obj.layer3.intrinsic) + 0.5*v/obj.ktq;
            else
                obj.voltage(1,1+obj.layer2.endPos) = -1*log(-1*obj.layer2.dopingDensity/obj.layer2.intrinsic) + 0.5*v/obj.ktq;
            end
            
                obj.startSim(obj.lastIteration);
                obj.voltageBiased = obj.voltage;

            end
            
            obj.voltage = tmp_zeroVoltage;
            
        
        
        end
        
        function currentDens2Voltage(obj)
            
            obj.nCons = obj.nCons ./ obj.concentration_scaling;
            obj.pCons = obj.pCons ./ obj.concentration_scaling;
            
            obj.v_diff = obj.voltage(1,2:obj.total_mesh_points+1) - obj.voltage(1,1:obj.total_mesh_points);
            obj.n_diff = obj.nCons(1,2:obj.total_mesh_points+1) - obj.nCons(1,1:obj.total_mesh_points);
            obj.p_diff = obj.pCons(1,2:obj.total_mesh_points+1) - obj.pCons(1,1:obj.total_mesh_points);
                        
            obj.NcurrentDens = obj.eMobility(1,2:obj.total_mesh_points+1) .* (obj.n_diff - obj.nCons(1,2:obj.total_mesh_points+1) .* obj.v_diff );
            obj.PcurrentDens = obj.hMobility(1,2:obj.total_mesh_points+1) .* (-1 .* obj.p_diff - obj.pCons(1,2:obj.total_mesh_points+1) .* obj.v_diff );
                    
            obj.derNcurrentDens = obj.NcurrentDens(1,2:obj.total_mesh_points) - obj.NcurrentDens(1,1:obj.total_mesh_points-1);
            obj.derPcurrentDens = obj.PcurrentDens(1,2:obj.total_mesh_points) - obj.PcurrentDens(1,1:obj.total_mesh_points-1); 
            
        end
        
        function radiativeRecombination(obj)

            bigB = 5.8e-13 .* sqrt(100 .* obj.permittivity .* obj.permittivity_scaling .* (obj.electron_mass ./ (obj.mct_electron_mass+obj.mct_hole_mass))) .* (obj.electron_mass ./ (obj.mct_electron_mass+obj.mct_hole_mass));
            bigB = bigB .* (1 + obj.electron_mass ./ obj.mct_electron_mass + obj.electron_mass ./ obj.mct_hole_mass);
            bigB = bigB .* (300/obj.Temp) .* sqrt(300/obj.Temp);
            bigB = bigB .* (obj.bandGap .* obj.bandGap + 3 .* obj.ktq .* obj.bandGap + 3.75 .* obj.ktq .* obj.ktq);
            bigB = bigB .* obj.time_scaling ./ ((obj.length_scaling).^3);
            
            obj.R_radiative = bigB .* (obj.nConsBiased .* obj.pConsBiased - obj.intrinsic .* obj.intrinsic .* obj.concentration_scaling .* obj.concentration_scaling);

%             obj.tau_radiative = tau ./ obj.time_scaling;
%             obj.R_radiative = (obj.nCons - obj.pCons - obj.intrinsic .* obj.intrinsic) .* obj.intrinsic .* obj.intrinsic ./ (obj.tau_radiative .* (obj.nCons + obj.pCons));
% 
%             obj.R_radiative = (obj.nConsBiased - obj.pConsBiased - obj.intrinsic .* obj.intrinsic) .* obj.intrinsic .* obj.intrinsic ./ (obj.tau_radiative .* (obj.nConsBiased + obj.pConsBiased));
%             
        end
        
        function srhRecombination(obj)
           
            
            
        end
        
        
        function chargesDerivative(obj)
            obj.nConsBiased = exp(obj.voltageBiased) .* exp(0.5*obj.V_bias/obj.ktq) .* obj.concentration_scaling .* obj.intrinsic;
            obj.pConsBiased = exp(-obj.voltageBiased) .* exp(0.5*obj.V_bias/obj.ktq) .* obj.concentration_scaling .* obj.intrinsic;
            obj.chargesBiased = obj.pConsBiased - obj.nConsBiased + obj.dopDensProf .* obj.concentration_scaling;
            obj.chargesBiased = -1 * obj.chargesBiased ./ obj.concentration_scaling;
            obj.secondDerVoltBiased = obj.voltageBiased(1,1:obj.total_mesh_points) + obj.voltageBiased(1,3:obj.total_mesh_points+2) - 2*obj.voltageBiased(1,2:obj.total_mesh_points+1);
            obj.secondDerVoltBiased(1,2:obj.total_mesh_points-1) = obj.secondDerVoltBiased(1,1:obj.total_mesh_points-2);
            
            tmp = obj.V_bias;
            obj.V_bias = 0.0;
            obj.nCons = exp(obj.voltage) .* exp(0.5*obj.V_bias/obj.ktq) .* obj.concentration_scaling .* obj.intrinsic;
            obj.pCons = exp(-obj.voltage) .* exp(0.5*obj.V_bias/obj.ktq) .* obj.concentration_scaling .* obj.intrinsic;
            obj.charges = obj.pCons - obj.nCons + obj.dopDensProf .* obj.concentration_scaling;
            obj.charges = -1 * obj.charges ./ obj.concentration_scaling;
            obj.secondDerVolt = obj.voltage(1,1:obj.total_mesh_points) + obj.voltage(1,3:obj.total_mesh_points+2) - 2*obj.voltage(1,2:obj.total_mesh_points+1);
            obj.secondDerVolt(1,2:obj.total_mesh_points-1) = obj.secondDerVolt(1,1:obj.total_mesh_points-2);
            obj.V_bias = tmp;
            obj.np_intrinsic_diffBiased = obj.nConsBiased .* obj.pConsBiased - obj.intrinsic .* obj.intrinsic .* obj.concentration_scaling .* obj.concentration_scaling;
            obj.np_intrinsic_diff = obj.nCons .* obj.pCons - obj.intrinsic .* obj.intrinsic .* obj.concentration_scaling .* obj.concentration_scaling;

            
        end
        
        function findTotCurrent(obj)
            
            charge_difference = obj.chargesBiased - obj.charges;
            charge_difference = charge_difference * obj.concentration_scaling;
    
            obj.dummyVar = charge_difference;
            
            mean_charge_diff = mean(charge_difference);
            mean_charge_diff = obj.electron_charge * mean_charge_diff;
            mean_charge_diff = mean_charge_diff * obj.currentDens_scaling;
            obj.totCurrent = mean_charge_diff;
        end
        
        function plotChargesDerv(obj)
            plot(obj.secondDerVolt)
            hold on
            plot(obj.charges)
        end
        
        function eField(obj)
            obj.firstDerVoltBiased = obj.voltageBiased(1,1:obj.total_mesh_points) - obj.voltageBiased(1,2:obj.total_mesh_points+1);
        end
        
        function currentDens(obj)
        
            obj.v_diff = obj.voltage(1,2:obj.total_mesh_points+1) - obj.voltage(1,1:obj.total_mesh_points);
            obj.n_diff = obj.nCons(1,2:obj.total_mesh_points+1) - obj.nCons(1,1:obj.total_mesh_points);
            obj.p_diff = obj.pCons(1,2:obj.total_mesh_points+1) - obj.pCons(1,1:obj.total_mesh_points);
            
            obj.v_diffBiased = obj.voltageBiased(1,2:obj.total_mesh_points+1) - obj.voltageBiased(1,1:obj.total_mesh_points);
            obj.n_diffBiased = obj.nConsBiased(1,2:obj.total_mesh_points+1) - obj.nConsBiased(1,1:obj.total_mesh_points);
            obj.p_diffBiased = obj.pConsBiased(1,2:obj.total_mesh_points+1) - obj.pConsBiased(1,1:obj.total_mesh_points);
            
            obj.NcurrentDens = obj.electron_charge .* obj.eMobility(1,2:obj.total_mesh_points+1) .* obj.mobility_scaling ./ obj.length_scaling .* ( obj.ktq .* obj.n_diff  - obj.nCons(1,2:obj.total_mesh_points+1) .* obj.v_diff ); 
            obj.PcurrentDens = obj.electron_charge .* obj.hMobility(1,2:obj.total_mesh_points+1) .* obj.mobility_scaling ./ obj.length_scaling .* ( -1 .* obj.ktq .* obj.p_diff  - obj.pCons(1,2:obj.total_mesh_points+1) .* obj.v_diff ); 
            
            obj.derNcurrentDens = obj.NcurrentDens(1,2:obj.total_mesh_points) - obj.NcurrentDens(1,1:obj.total_mesh_points-1);
            obj.derPcurrentDens = obj.PcurrentDens(1,2:obj.total_mesh_points) - obj.PcurrentDens(1,1:obj.total_mesh_points-1); 
            
            obj.NcurrentDensBiased = obj.electron_charge .* obj.eMobility(1,2:obj.total_mesh_points+1) .* obj.mobility_scaling ./ obj.length_scaling .* ( obj.ktq .* obj.n_diffBiased  - obj.nConsBiased(1,2:obj.total_mesh_points+1) .* obj.v_diffBiased ); 
            obj.PcurrentDensBiased = obj.electron_charge .* obj.hMobility(1,2:obj.total_mesh_points+1) .* obj.mobility_scaling ./ obj.length_scaling .* ( -1 .* obj.ktq .* obj.p_diffBiased  - obj.pConsBiased(1,2:obj.total_mesh_points+1) .* obj.v_diffBiased ); 

            obj.derNcurrentDensBiased = obj.NcurrentDensBiased(1,2:obj.total_mesh_points) - obj.NcurrentDensBiased(1,1:obj.total_mesh_points-1);
            obj.derPcurrentDensBiased = obj.PcurrentDensBiased(1,2:obj.total_mesh_points) - obj.PcurrentDensBiased(1,1:obj.total_mesh_points-1);
            
            obj.totderCurrDensBiased = obj.NcurrentDensBiased + obj.PcurrentDensBiased;
            
            obj.dummyVar = obj.nConsBiased(1,1:obj.total_mesh_points) .* obj.p_diffBiased + obj.pConsBiased(1,1:obj.total_mesh_points) .* obj.n_diffBiased + obj.p_diffBiased .* obj.n_diffBiased;
            obj.dummyVar1 = obj.nConsBiased(1,1:obj.total_mesh_points) .* obj.p_diffBiased;
            obj.dummyVar2 = obj.pConsBiased(1,1:obj.total_mesh_points) .* obj.n_diffBiased;
            obj.dummyVar3 = obj.p_diffBiased .* obj.n_diffBiased;
            
            
        end
        
        
        
%         function startTempParamSim(obj, temp_start, temp_end, pointNumber)
%             i = 1;
%             
%             for i:pointNumber
%                 
%                 
%                 
%             end
%             
%             
%         end
%         
        
 
    end
end








