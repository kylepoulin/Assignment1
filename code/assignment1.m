% clear all
clearvars
clearvars -GLOBAL
close all
format shorte

set(0, 'DefaultFigureWindowStyle', 'docked')
global restmass boltzmann vth;
restmass = 9.109*10^-31;
boltzmann = 1.3806*10^-23;
vth = sqrt((boltzmann)*300/(restmass*0.26));



global Q1;
Q1 = false;

global Q2;
Q2 = false;

global Q3;
Q3 = true;

TotalTime = 1000;
global numElectrons;
numElectrons = 1000;
global T;
T = 0;

global Temperature;
global oldTemperature;

global xmax ymax;
xmax = 200*10^-9;
ymax = 100*10^-9;


boundary = zeros(2,4);
checkableboundary = zeros(2,4);
%2 boundaries, 1=xlower, 2=ylower, 3=xhigher, 4=yhigher
numbounds=2;

boundary(1,1)=(xmax/2)-0.125*xmax;
boundary(1,2)=0;
boundary(1,3)=0.125*xmax;
boundary(1,4)=0.4*ymax;

boundary(2,1)=(xmax/2)-0.125*xmax;
boundary(2,2)=ymax-0.4*ymax;
boundary(2,3)=0.125*xmax;
boundary(2,4)=0.4*ymax;

checkableboundary(1,1)=(xmax/2)-0.125*xmax;
checkableboundary(1,2)=0;
checkableboundary(1,3)=boundary(1,1)+0.125*xmax;
checkableboundary(1,4)=boundary(1,2)+0.4*ymax;

checkableboundary(2,1)=(xmax/2)-0.125*xmax;
checkableboundary(2,2)=ymax-0.4*ymax;
checkableboundary(2,3)=boundary(2,1)+0.125*xmax;
checkableboundary(2,4)=boundary(2,1)+0.4*ymax;


global TimeStep;
TimeStep = (xmax/vth)/100;

global ElectronInfo;

ElectronInfo = zeros(numElectrons,4);
%First2  layers = x & y positions
%Layers 3 = velocity
%layers 4 =  angles
if Q3 == true
    for k = 1:numElectrons
        firstX = xmax*rand;
        firstY = ymax*rand;
        while ((firstX > checkableboundary(1,1)&& firstX < checkableboundary(1,3) && firstY>checkableboundary(1,2) && firstY<checkableboundary(1,4)) ...
                || (firstX > checkableboundary(2,1)&& firstX < checkableboundary(2,3) && firstY>checkableboundary(2,2) && firstY<checkableboundary(2,4)))
            firstX = xmax*rand;
            firstY = ymax*rand;
        end
        ElectronInfo(k,1) = firstX;
        ElectronInfo(k,2) = firstY;
    end
else
    ElectronInfo(:,1) = xmax*rand(1,numElectrons);
    ElectronInfo(:,2) = ymax*rand(1,numElectrons);
end

figure(1)

if Q1 == true
    ElectronInfo(:,3) = vth;
    ElectronInfo(:,4) = 2*pi*rand(1,numElectrons);
elseif Q2 == true || Q3 ==true
    for i = 1:numElectrons
        vy = vth*randn*0.8;
        vx = vth*randn*0.8;
        ElectronInfo(i,4) = tan(vy/vx);
        ElectronInfo(i,3) = sqrt(vx^2 + vy^2);
    end
    mean(ElectronInfo(:,3))
end

hist(ElectronInfo(:,3),40);
title('Velocity Distribution');
xlabel('velocity (m/s)');
ylabel('Number in bins');
xlim([0 5*vth]);
pause(0.05);



global electronColours;
electronColours = hsv(10);

hitrightbound = false;
hitleftbound = false;
hitupperbound = false;
hitlowerbound = false;

MFPcalc = zeros(numElectrons,4);
MFPsum = 0;
collisionNum = 0;
collision = false;

pscat = 1 - exp(-TimeStep/(0.2*10^-12));
clearcounter = 0;

collisiontimes = zeros(numElectrons,1);
collisiontimesum = 0;
if Q3 == true
    figure(2)
    for j = 1:numbounds
        rectangle('Position',boundary(j,:));
        hold on;
    end
    xlim([0 xmax]);
    ylim([0 ymax]);
    pause(0.1);
end


while T<(TimeStep*1000)
    T = T + TimeStep
    clearcounter = clearcounter+1;
    figure(2)
    pause(0.05);
    oldTemperature = Temperature;
    
    for i = 1:numElectrons
        %X position update
        if Q2==true  || Q3==true          
            if pscat>rand
                MFPsum = MFPsum + sqrt(((MFPcalc(i,1))^2+(MFPcalc(i,2))^2));
                MFPcalc(i,1) = 0;
                MFPcalc(i,2) = 0;
                collisiontimesum = collisiontimesum + collisiontimes(i,1);
                collisiontimes(i,1) = 0;
                vy = vth*randn*0.8;
                vx = vth*randn*0.8;
                ElectronInfo(i,4) = tan(vy/vx);
                ElectronInfo(i,3) = sqrt(vx^2 + vy^2);
                collisionNum = collisionNum + 1;
            end
        end
        %xposition update
        oldx = ElectronInfo(i,1);
        newX = oldx + TimeStep*ElectronInfo(i,3)*cos(ElectronInfo(i,4));
        
        if newX > xmax
            hitrightbound = true;            
        elseif newX < 0
            hitleftbound = true;            
        end
        

        %Y position update
        oldy = ElectronInfo(i,2);
        oldangle = ElectronInfo(i,4);
        newY = oldy + TimeStep*ElectronInfo(i,3)*sin(ElectronInfo(i,4));
        if newY > ymax
            hitupperbound = true;            
        elseif newY < 0
            hitlowerbound = true;
        end
        
        userbound.triggered = false;
        userbound.which = [1 1];
        if Q3 == true
            for k = 1:2
                if (newX>checkableboundary(k,1))&& (newX<checkableboundary(k,3)) && (newY>checkableboundary(k,2)) && (newY<checkableboundary(k,4))
                   userbound.triggered = true;
                   
                   if oldx < checkableboundary(k,1)
                       hitrighttbound = true
                       userbound.which = [k 1]
                   elseif oldx > checkableboundary(k,3)
                       hitleftbound = true;
                       userbound.which = [k 3]
                   elseif oldy < checkableboundary(k,2)
                       hitupperbound = true;
                       userbound.which = [k 2]
                   elseif oldy > checkableboundary(k,4)
                       hitlowerbound = true;
                       userbound.which = [k 4]
                   end
                end
            end
        end
        
        
        
        if Q2==true || Q3 == true
            MFPcalc(i,1) = MFPcalc(i,1)+(newX-oldx);
            MFPcalc(i,2) = MFPcalc(i,2)+(newY-oldy);
            collisiontimes(i,1) = collisiontimes(i,1)+TimeStep;
        end
        
        if hitrightbound == true            
            if (userbound.triggered == true) && Q3 ==true
                boundaryY = oldy+(checkableboundary(userbound.which(1),userbound.which(2))-oldx)*tan(oldangle);
                newX = checkableboundary(userbound.which(1),userbound.which(2))- (newX-checkableboundary(userbound.which(1),userbound.which(2)));
                ElectronInfo(i,4) = ElectronInfo(i,4)+(90*pi/180);
                userbound.triggered = false;
            else
                newX = newX - xmax;
                boundaryY = oldy+(xmax-oldx)*tan(oldangle);
            end
            if i<10
                plot([oldx xmax], [oldy boundaryY], 'color', electronColours(i,:));
                hold on;
            end
            oldx = 0;
            oldy = boundaryY;
            hitrightbound = false;
            
        elseif hitleftbound == true                        
            if userbound.triggered == true && Q3 ==true
                newX =  2*checkableboundary(userbound.which(1),userbound.which(2))-newX;
                boundaryY = oldy+(oldx - checkableboundary(userbound.which(1),userbound.which(2)))*tan(oldangle);
                ElectronInfo(i,4) = oldangle+(90*pi/180);
                userbound.triggered = false;
            else
                newX = newX + xmax;
                boundaryY = oldy+(oldx)*tan(oldangle);
            end
            if i<10
                plot([oldx 0], [oldy boundaryY], 'color', electronColours(i,:));
                hold on;
            end
            oldx = xmax;
            oldy = boundaryY;
            hitleftbound = false;
            
        elseif hitupperbound == true      
            
            if userbound.triggered == true && Q3 ==true
                boundaryX = oldx+(checkableboundary(userbound.which(1),userbound.which(2))-oldy)/tan(oldangle);
                newY = checkableboundary(userbound.which(1),userbound.which(2)) - (newY-checkableboundary(userbound.which(1),userbound.which(2)));
                ElectronInfo(i,4) = -ElectronInfo(i,4);
                userbound.triggered = false;
            else
                boundaryX = oldx+(ymax-oldy)/tan(oldangle);
                newY = ymax - (newY-ymax);
                ElectronInfo(i,4) = -ElectronInfo(i,4);
            end
            if i<10
                plot([oldx boundaryX], [oldy ymax], 'color', electronColours(i,:));
                hold on;
            end
            oldx = boundaryX;
            oldy = ymax;
            hitupperbound = false;
        elseif hitlowerbound == true
            
            if userbound.triggered == true && Q3 ==true
                newY = 2*checkableboundary(userbound.which(1),userbound.which(2))-newY;
                ElectronInfo(i,4) = -ElectronInfo(i,4);
                boundaryX = oldx+(oldy-checkableboundary(userbound.which(1),userbound.which(2)))/tan(oldangle);
                userbound.triggerd = false;
            else
                newY = -newY;
                ElectronInfo(i,4) = -ElectronInfo(i,4);
                boundaryX = oldx+(oldy)/tan(oldangle);
            end
            if i<10
                plot([oldx boundaryX], [oldy 0], 'color', electronColours(i,:));
                hold on;
            end
            oldx = boundaryX;
            oldy = 0;
            hitlowerbound = false;
        end
        
        %update electron info
        ElectronInfo(i,1) = newX;
        ElectronInfo(i,2) = newY;
        
        if i<10
            plot([oldx newX], [oldy newY], 'color', electronColours(i,:));
            xlim([0 xmax]);  
            ylim([0 ymax]);
            Temperature = ((mean(ElectronInfo(:,3)))^2)*(0.26*restmass)/boltzmann;
            if Q1 ==true
                xlabel(['Temperature: ' num2str(Temperature) ' Kelvin']);
            elseif Q2==true
                xlabel(['Temperature: ' num2str(Temperature) ' Kelvin, MFP: ' num2str(MFPsum/collisionNum) 'm, TBC: ' num2str(collisiontimesum/collisionNum) 's']);
            end
            hold on;
            pause(0.01);
        end
    end
   %{ 
    if clearcounter == 30;
        clf;
        
        if Q3 == true
            figure(2)
            for j = 1:numbounds
                rectangle('Position',boundary(j,:));
                hold on;
            end
            xlim([0 xmax]);
            ylim([0 ymax]);
            pause(0.1);
        end
        
        clearcounter =0;
    end
    %}
    figure(3)
    plot([T-TimeStep T],[oldTemperature Temperature], 'Color', 'Blue');
    title('Temperature');
    ylabel('degrees Kelvin');
    ylim([0 500]);
    hold on;
    pause(0.05);
    
  
    figure(1)
    hist(ElectronInfo(:,3),40);
    title('Velocity Distribution');
    xlabel('velocity (m/s)');
    ylabel('Number in bins');
    xlim([0 5*vth]);
    pause(0.05);

end