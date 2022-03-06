classdef Receiver
    
    properties
        
        %The Usual
        pi_data;
        c_light;
        R_earth;
        sidereal;
        
        %Constructed Satellite, Xi, Xs, ts
        s = Satellite();
        
        %Keep track of the satellites
        numberOfSatellites;
        satelliteIndices;
        satellitePositions;
        satelliteTimes;
        
        %An initial guess obtained from the vehicle data
        vehValues;
    end
    
    
    methods
       
        function obj = Receiver()
        
            %Define Properties
            fidd = fopen('data.dat','r');
            formatSpec = '%f %*[^\n]';
            A = textscan(fidd,formatSpec);
            fclose(fidd);
            ddat = A{1};
            
            obj.pi_data = ddat(1);
            obj.c_light = ddat(2);
            obj.R_earth = ddat(3);
            obj.sidereal = ddat(4);
            
            %Read in Satellite Data
            fids = fopen('satellite_temp.log', 'r');
            C = fscanf(fids, '%f');
            fclose(fids);
            C = C';
            obj.vehValues = C(1:10);
            C = C(11:length(C));
            obj.numberOfSatellites = length(C)/5;
            obj.satelliteIndices = zeros(obj.numberOfSatellites, 1);
            obj.satelliteTimes = zeros(obj.numberOfSatellites, 1);
            satelliteX = zeros(obj.numberOfSatellites, 1);
            satelliteY = zeros(obj.numberOfSatellites, 1);
            satelliteZ = zeros(obj.numberOfSatellites, 1);
            
            j = 1;
            for k = 1:obj.numberOfSatellites
                obj.satelliteIndices(k) = C(j);
                obj.satelliteTimes(k) = C(j+1);
                satelliteX(k) = C(j+2);
                satelliteY(k) = C(j+3);
                satelliteZ(k) = C(j+4);
                obj.satellitePositions(k, :) = [satelliteX(k) satelliteY(k) satelliteZ(k)];
                j = j + 5;
            end           
        end
        
        %A function to convert from radians to DMS
        function DMS = convertFromRadians(Obj,radians)
            DMS(4) = 1;
            if (radians < 0)
                DMS(4) = -1;
                radians = -radians;
            end
            decimalDeg = radians*180/Obj.pi_data;
            degrees = floor(decimalDeg);
            minutes = mod(floor((decimalDeg - degrees) * 60),60);
            seconds = round((decimalDeg - degrees - (minutes / 60)) * 3600, 2);
            if (seconds ~= 0 && mod(seconds, 60) == 0) 
                minutes = minutes + 1;
                seconds = 0;
            end
            if (minutes == 60)
                degrees = degrees + 1;
                minutes = 0;
            end
            DMS(1) = degrees;
            DMS(2) = minutes;
            DMS(3) = seconds;
        end
        
                %A function to convert DMS to Radians
        function radians = convertToRadians(Satellite,degrees, minutes, seconds, sign)
            deg = degrees + (minutes/60) + (seconds/3600);
            radians = deg*(Satellite.pi_data)*sign/180;
        end
        
        %A function to convert cartesian coordinates at some time t to DMS
        %coordinates in longitude and latitude
        function values = convertCartesianToDMS(Obj, xyzt)
            x = round(xyzt(1), 2);
            y = round(xyzt(2), 2);
            z = round(xyzt(3), 3);
            longRad = wrapToPi(atan2(y, x) - (2 * Obj.pi_data / Obj.sidereal * xyzt(4))); 
            xyztRadius = sqrt(x^2 + y^2 + z^2);
            altitude = round(xyztRadius - Obj.R_earth, 2);
            latRad = asin(z/xyztRadius);
            latDMS = convertFromRadians(Obj,latRad);
            longDMS = convertFromRadians(Obj,longRad);
            time = round(xyzt(4), 11);
            
            values(1) = time;
            values(2) = latDMS(1);
            values(3) = latDMS(2);
            values(4) = latDMS(3);
            values(5) = latDMS(4);
            values(6) = longDMS(1);
            values(7) = longDMS(2);
            values(8) = longDMS(3);
            values(9) = longDMS(4);
            values(10) = altitude;
        end
        
        
        %A function to convert DMS coordinates in latitude and longitude to
        %cartesian coordinates
        function xyz = convertGivenToCartesian(Satellite, values)
            latRad = convertToRadians(Satellite, values(2), values(3), values(4), values(5));
            longRad = convertToRadians(Satellite, values(6), values(7), values (8), values(9));
            
            xyz(1) = cos(latRad) * (Satellite.R_earth + values(10)) * cos(2*Satellite.pi_data/Satellite.sidereal*values(1) + longRad);
            xyz(2) = cos(latRad) * (Satellite.R_earth + values(10)) * sin(2*Satellite.pi_data/Satellite.sidereal*values(1) + longRad);
            xyz(3) = (Satellite.R_earth + values(10)) * sin(latRad);
        end
        
        %If the receiver only has 4 satellites we use this method. Simply
        %grab the 4 satellites from the receiver along with the 4 times. No
        %need to feed in anything besides the object
        function solution = solveWFourSat(obj)
           
            XsMatrix = obj.satellitePositions;
            tsVector = obj.satelliteTimes;
            %Position Vectors
            sat1 = XsMatrix(1,:);
            sat2 = XsMatrix(2,:);
            sat3 = XsMatrix(3,:);
            sat4 = XsMatrix(4,:);
            
            F = zeros(4,1);
            J = zeros(4,4);
            
            %Now We Do Newton's Method For A System.
            posGuess = convertGivenToCartesian(obj, obj.vehValues);
            
            %Our initial guess is a little too good so let's shift it a bit
            posGuess = posGuess + [100 100 100];
            initSol = [posGuess tsVector(1) + 0.1];
            initSol = initSol';
            maxIt = 100;
      
            for k = 1:maxIt
                
                %Extract
                timeEst = initSol(4);
                posEst = [initSol(1) initSol(2) initSol(3)];
               
                %Make F
                F(1) = getValueOfFunction(obj, posEst, sat1, timeEst, tsVector(1));
                F(2) = getValueOfFunction(obj, posEst, sat2, timeEst, tsVector(2));
                F(3) = getValueOfFunction(obj, posEst, sat3, timeEst, tsVector(3));
                F(4) = getValueOfFunction(obj, posEst, sat4, timeEst, tsVector(4));
                
                %Make J
                for j = 1:4
                    xs = XsMatrix(j,:);
                    J(j,1) = partialWithRespectToX(obj, posEst, xs);
                    J(j,2) = partialWithRespectToY(obj, posEst, xs);
                    J(j,3) = partialWithRespectToZ(obj, posEst, xs);
                    J(j,4) = -obj.c_light;
                end
                
                %Solve the system and make our new guess
                y=J\F;
                newGuess = initSol - y;
                
                %Check if the new guess is within a certain tolerance
                if abs(newGuess(4) - initSol(4)) < 10^-11
                    solution = newGuess;
                    return;
                else
                    initSol = newGuess;
                    continue;
                end
            end
            
            %If we can't find a solution, indiciate that to the user in
            %some way
            disp("No solution was found :(");
            solution = newGuess;
        end
        
        %Evaluate the distance between xv and xs
        function value = getDistanceBetweenPositions(obj, xv, xs)
           value = norm(xv - xs); 
        end
        
        %Evaluate the function ||xv - xs|| - c*(tv - ts)
        function value = getValueOfFunction(obj, xv, xs, tv, ts)
            c = obj.c_light;
            norm = getDistanceBetweenPositions(obj, xv, xs);
            secondPart = c*(tv - ts);
            value = norm - secondPart;
        end
        
        %Take the partial derivative with respect to x using finite
        %distance
        function value = partialWithRespectToX(obj, xv, xs)
            norm = getDistanceBetweenPositions(obj, xv, xs);
            x_vehicle = xv(1);
            x_satellite = xs(1);
            value = (x_vehicle - x_satellite)/norm;
        end
        
        %Take the partial derivative with respect to y using finite distance
        function value = partialWithRespectToY(obj, xv, xs)
            norm = getDistanceBetweenPositions(obj, xv, xs);
            y_vehicle = xv(2);
            y_satellite = xs(2);
            value = (y_vehicle - y_satellite)/norm;
        end
        
        %Take the partial derivative with respect to z using finite
        %distance
        function value = partialWithRespectToZ(obj, xv, xs)
            norm = getDistanceBetweenPositions(obj, xv, xs);
            z_vehicle = xv(3);
            z_satellite = xs(3);
            value = (z_vehicle - z_satellite)/norm;
        end
        
        function value = lsqFWTime(obj, posEst, timeEst)
            availSats = obj.numberOfSatellites;
            sum = 0;
            for i = 1:availSats
                sat1 = obj.satellitePositions(i,:);
                time1 = obj.satelliteTimes(i);
                firstNorm = norm(posEst - sat1);
                timeDif = obj.c_light * (timeEst - time1);
                result = (firstNorm - timeDif)^2;
                sum = sum + result;
            end
            value = sum;
        end
        
        
        %Evaluate the first partial of the lsq function above with respect
        %to the variable passed in
        function firstPartial = lsqPartial(obj, posEst,timeEst, variable)
            initSum = 0;
            sats = obj.numberOfSatellites;
            for i=1:sats
                satPos = obj.satellitePositions(i, :);
                dist = norm(posEst - satPos);
                timeDiff = obj.c_light * (timeEst - obj.satelliteTimes(i));
                if (variable == 'x')
                    lastTerm = (posEst(1) - satPos(1))/dist;
                elseif (variable == 'y')
                    lastTerm = (posEst(2) - satPos(2))/dist;
                elseif (variable == 'z')
                    lastTerm = (posEst(3) - satPos(3))/dist;
                elseif (variable == 't')
                    lastTerm = -2*obj.c_light;
                end
                result = (dist - timeDiff)*lastTerm;
                initSum = initSum + result;
            end
            firstPartial = initSum;
        end
        
        %If the receiver has n satellites we use this method. 
        function solution = solveWNSat(obj)

            F = zeros(4,1);
            J = zeros(4,4);
            
            %Now We Do Newton's Method For A System.
            posGuess = convertGivenToCartesian(obj, obj.vehValues);
            
            %Our initial guess is a little too good so let's shift it a bit
            posGuess = posGuess + [100 100 100];
            timeGuess = obj.satelliteTimes(1);
            maxIt = 10;
      
            initGuess = [posGuess timeGuess];
            for k = 1:maxIt
                               
                %Make F
                F(1) = lsqPartial(obj, posGuess, timeGuess, 'x');
                F(2) = lsqPartial(obj, posGuess, timeGuess, 'y');
                F(3) = lsqPartial(obj, posGuess, timeGuess, 'z');
                F(4) = lsqPartial(obj, posGuess, timeGuess, 't');
                
                %Make J
                J(1,1) = secondPartial(obj, posGuess, timeGuess, 'x', 'x');
                J(1,2) = secondPartial(obj, posGuess, timeGuess, 'x', 'y');
                J(1,3) = secondPartial(obj, posGuess, timeGuess, 'x', 'z');
                J(1,4) = secondPartial(obj, posGuess, timeGuess, 'x', 't');
                J(2,1) = secondPartial(obj, posGuess, timeGuess, 'y', 'x');
                J(2,2) = secondPartial(obj, posGuess, timeGuess, 'y', 'y');
                J(2,3) = secondPartial(obj, posGuess, timeGuess, 'y', 'z');
                J(2,4) = secondPartial(obj, posGuess, timeGuess, 'y', 't');
                J(3,1) = secondPartial(obj, posGuess, timeGuess, 'z', 'x');
                J(3,2) = secondPartial(obj, posGuess, timeGuess, 'z', 'y');
                J(3,3) = secondPartial(obj, posGuess, timeGuess, 'z', 'z');
                J(3,4) = secondPartial(obj, posGuess, timeGuess, 'z', 't');
                J(4,1) = secondPartial(obj, posGuess, timeGuess, 't', 'x');
                J(4,2) = secondPartial(obj, posGuess, timeGuess, 't', 'y');
                J(4,3) = secondPartial(obj, posGuess, timeGuess, 't', 'z');
                J(4,4) = secondPartial(obj, posGuess, timeGuess, 't', 't');
                
                %Solve for y
                y=J\F;
                newGuess = initGuess - y';
                
                if norm(newGuess - initGuess) < 10^-2
                    initGuess = newGuess;
                    break;
                else
                    initGuess = newGuess;
                    posGuess = initGuess(1:3);
                    timeGuess = initGuess(4);
                    continue;
                end
            end
            finalGuess = initGuess;
            solution = finalGuess';
            finalValues = convertCartesianToDMS(obj, solution);
        end
        
        %A method to find the second partial derivatives for the Jacobian
        %Matrix. Utilizeds the firstPartial and uses finite difference.
        function value = secondPartial(obj, posEst, timeEst, variable1, variable2)
            %Define a time step for the appropriate variable
            h = 0.001;
            if variable2 == 'x'
                step = [h 0 0];
            elseif variable2 == 'y'
                step = [0 h 0];
            elseif variable2 == 'z'
                step = [0 0 h];
            end
            
            %First partial with respect to x and second with respect to
            %variable2
            if variable1 == 'x' 
                if variable2 ~= 't'
                    value = (lsqPartial(obj, posEst + step, timeEst, 'x') - lsqPartial(obj, posEst, timeEst, 'x')) / h;
                else
                    value = (lsqPartial(obj, posEst, timeEst + h, 'x') - lsqPartial(obj, posEst, timeEst, 'x'))/ h;
                end
            end
            
            %First partial with respect to y and second with respect to
            %variable2
            if variable1 == 'y'
                if variable2 ~= 't'
                    value = (lsqPartial(obj, posEst + step, timeEst, 'y') - lsqPartial(obj, posEst, timeEst, 'y')) / h;
                else
                    value = (lsqPartial(obj, posEst, timeEst + h, 'y') - lsqPartial(obj, posEst, timeEst, 'y'))/ h;
                end
            end
            
            %First partial with respect to z and second with respect to
            %variable2.
            if variable1 == 'z'
                if variable2 ~= 't'
                    value = (lsqPartial(obj, posEst + step, timeEst, 'z') - lsqPartial(obj, posEst, timeEst, 'z')) / h;
                else
                    value = (lsqPartial(obj, posEst, timeEst + h, 'z') - lsqPartial(obj, posEst, timeEst, 'z'))/ h;
                end
            end
            
            %First partial with respect to t and second with respect to
            %variable2
            if variable1 =='t'
                if variable2 ~= 't'
                    value = (lsqPartial(obj, posEst + step, timeEst, 't') - lsqPartial(obj, posEst, timeEst, 't')) / h;
                else
                    value = (lsqPartial(obj, posEst, timeEst + h, 't') - lsqPartial(obj, posEst, timeEst, 't'))/ h;
                end
            end
        end
        
        %Send the necessary data to the receiver_temp.log file
        %Checks how many satellites are available and uses the appropriate
        %method to solve the system
        function receiverLog = sendData(obj)
            availableSatellites = obj.numberOfSatellites;
            if availableSatellites < 4
                disp("Not enough satellites to determine vehicle position");
                return;
            elseif availableSatellites == 4
                vehXYZT = solveWFourSat(obj);
            else
                vehXYZT = solveWNSat(obj);
            end
            vehicleValues = convertCartesianToDMS(obj, vehXYZT);
            fileID = fopen('receiver_temp.log','w');
            fprintf(fileID, '%.11f %d %.2f %.2f %d %d %.2f %.2f %d %.2f\n', vehicleValues);
            loggedData = num2str(vehicleValues, '%.11f %d %.2f %.2f %d %d %.2f %.2f %d %.2f\n');
            loggedData = split(loggedData)'; 
            fclose(fileID);
            receiverLog = loggedData;
        end
    end
    
    
end

