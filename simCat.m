
classdef simCat < handle
   
  properties (Constant)
        electron_charge = 1.6e-19; % Coulomb - charge of electron
        boltzmann_const = 1.38064852e-23; % m2.kg.s-2.K-1 - Boltzmann constant
        vacuum_permittivity = 8.85 * 10^(-14); % Farad/cm - Vacuum permittivity
  end
    
    properties % Layer Definition and Environment
        V_bias = 0.0; % Default biasing value
        Temp = 0; % 0K defined for user to submit
        numofLayers = 0; % Initial number of layers  
        numofSims = 0; % Initial number of simulations
        
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
        isScaled = false;
        total_mesh_points = 0;
        dopDensProf;
        permittivity;
        intrinsic;
        voltage;
        voltageBiased;
        voltagePosition;
        eMobility;
        hMobility;
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
        
        paramTotCurrent;

        dummyVar;
        
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
        
        function simHandle(obj)
           
            if obj.numofLayers == 3
                obj.numofSims = 2;
                % sim1 start ===========
                obj.sim1 = numSimPoisson;
                obj.sim1.setSimParameters(obj.V_bias, obj.Temp);
                % =============================
                obj.sim1.newLayer;
                obj.sim1.layer1 = obj.layer1;
                obj.sim1.endLayer;
                % ============================
                obj.sim1.newLayer;
                obj.sim1.layer2 = obj.layer2;
                obj.sim1.endLayer;
                % sim1 end ===============
                
                
                % sim2 start ==================
                obj.sim2 = numSimPoisson;
                obj.sim2.setSimParameters(obj.V_bias, obj.Temp);
                % ===============================
                obj.sim2.newLayer;
                obj.sim2.layer1 = obj.layer2;
                obj.sim2.endLayer;
                % ===============================
                obj.sim2.newLayer;
                obj.sim2.layer2 = obj.layer3;
                obj.sim2.endLayer;
                % sim2 end ====================
                
            else
                if obj.numofLayers == 2
                    obj.numofSims = 1;
                    obj.sim1 = numSimPoisson;
                    
                    obj.sim1.setSimParameters(obj.V_bias, obj.Temp);
                    obj.sim1.newLayer;
                    obj.sim1.layer1 = obj.layer1;
                    obj.sim1.endLayer;
                    
                    obj.sim1.newLayer;
                    obj.sim1.layer2 = obj.layer2;
                    obj.sim1.endLayer;


            
                end
                
                
            end
            

            
        end
        
        function deriveScaling(obj)

            obj.sim1.deriveScaling;
            if obj.numofSims == 2
                obj.sim2.deriveScaling;
            end
            
        end
        
        function initConditions(obj)
            
            obj.sim1.initConditions;
            if obj.numofSims == 2
                obj.sim2.initConditions;
            end
            
        end
        
        function startSim(obj, iteration)
            
            obj.lastIteration = iteration;
            
            obj.sim1.zeroBiasSim(iteration);
            obj.sim1.biasClosing(iteration);
            
            if obj.numofSims == 2
                obj.sim2.nonBiasSim;
                obj.sim2.startSim(obj.lastIteration);
            end
            
            
        end
        
        function voltageConcatanate(obj)
            
            
            
        end
        
        
    end
end