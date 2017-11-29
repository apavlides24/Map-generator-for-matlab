% This program creates a continent or island for RPG games using a simple
% form of the Plurigaussian simulation. 
% The resulting png is composed by pixels, with different colors for the 
% various terrain types. Each pixel represents an area that the appropriate
% terrain is dominant. This program -doesn't- place cities or rivers.
% There is also matrix "VAL" with values representing the terrain type for 
% each pixel. The png (or matrix) can be used with different programs to 
% trace or convert them to hex maps.
% Personally, I use the image with the Hexographer program 
% (http://www.hexographer.com) to "convert Underlay" and create hex maps.
%
% The program creates a "base" map first as a trend and uses this to create
% a series of more detailed maps that it then joins together using some
% basic (and crude) smoothing techniques and corrections.
%
% ATTENTION: This program is crudely made, and it suits me and my needs. It
% requires time and RAM to run and it is no way optimized.
% For example, if you want to change any of the parameters you need to get
% in the code and change them as this program is not a function with
% inputs, nargin etc.
%
% DEPENDENCIES: This program requires the statistics toolbox and image 
% processing toolbox of Matlab
%
% This program uses the exponential covariance so that every point on the
% grid has some correlation (even if very small) to all other points. The
% exponential model is the simplest to implement but probably not the best.
% Bibliography suggests that 2D random-walk correlation would be the best
% fit but this program was not created specifically to create terrain maps. 
%
% If you want to use sparse matrices to save on space in order to create
% trully huge maps, then I suggest you use the spherical covariance model
% and triple the correlation lenghts. 
%
% EXAMPLE 1: You can run the program as is and get five maps of a big temperate
% island, without tundras or deserts. 240 x 240 pixels, or 480 x 480 miles.
% It takes my computer 3-4GB of RAM and 14 minutes to make the necessary 
% covariance matrix and once the matrix is calculated then it takes 3-4 
% mins to make the simulations.
% Time varies depending on how many "Attempts" are required to get a good
% simulation.
%
% EXAMPLE 2: You can change SyntheM value to SyntheM=10, change North and
% South to North=1 and South=1, change corl to corl=2 and have 5 large
% continents, comparable to Europe in size with cold climate on top and warm
% climate at the bottom. The resulting map will be 800x800 pixels or 2400 x
% 2400 miles.
% It takes my computer 3-4GB of RAM and about 2 hours to finish the 5
% continents.
% Time varies depending on how many "Attempts" are required to get a good
% simulation.
% 
% PARAMETERS:
% Since the program was made for me, it's not a function with inputs, so to
% change anything you have to change it directly from the code.
% Explanations of the parameters are in comments.

t0=tic;
SyntheM=3; % The final map will be SyntheM x SyntheM basic maps. 
% I.e if you have Nx = 80 pixels and SyntheM=3, your final map will be 240 x 240 pixels

for LOOP=1:5;    

% ******************  PARAMETERS ********************
    
    
Nx=80; % The basic map will be Nx x Nx pixels (i.e. 80 x 80). 
% Nx is combined with the SyntheM to give the final map size. Nx should be
% above 40. Increasing Nx, -very- significantly affects time and RAM
% requirements. I find that Nx=90 is a good compromise and my computer's
% 16GB RAM are barely enough for Nx=120.
imrot=0; % imrot=1 rotates the landmass. This  rotation makes the final 
% image larger
DiffPlain=0; %  DiffPlain=1 makes different plains terrains (steppe and 
% savanna) available
corl=3; % corl multiplies the correlation length.
% Larger values for corl => less change with pixel distance => smaller scale (less rough change)
% Should be 3 for maps with a scale of 2-3 miles per pixel
% corl=1 or less could make the changes in altitude so rapid per pixel that 
% hills may be mostly invisible.

Sea=2; % How much sea would be surrounding the landmass. 
% sea=0 is no sea, sea=1 is low sea, sea=2 (max) is more sea.
StSea=0; % StSea=1 puts patches of sea at two opposing corners of the map 
% to create a less "square" landmass for islands.

North=0; South=0; % North and South allow cold terrain types or hot climate. 
% By default, cold terrain is towards the North and deserts, jungles etc 
% towards the south. If you don't want Deserts and jungles but want a
% colder North, use South=0 and North=1. If you want desert and jungles but
% not tundras, or ice-capped mountains use South=1 and North=0.
COLDF= -2.3; HEATF=-2; % COLDF and HEATF control the climate (in case of 
% North and South). The names ARE NOT indicative of what they do.
%
% Temperature is controlled by a variable called VAL3. The higher VAL3 is
% for a pixel, the higher the temperature.
% Cold climates start at VAL3<2.9 and warm climates start at VAL3>6.9. 
% A temperate zone is between VAL3 = 2.9 and VAL3= 6.9
% COLDF: The higher it is, the Colder the climate. With COLDF=0, VAL3
% starts at 0 at the North and increases linearly towards the South.
% with the default COLDF=-2.3, VAL3 starts at 1.15 and increases linearly.
% HEATF: Controls the range of VAL3. With HEATF=0, the range of VAL3 is 10.
% So with COLDF=0, VAL3 would go from 0 to 10. With HEATF=5 the range of
% VAL3 is 15. So with COLDF=0, VAL3 would go from 0 to 15. More than half
% the map in this case would be in warm climate, while the Northern third
% would be in cold terrain leaving only a small temperate zone.
% Putting COLDF at 5.8 and HEATF=-6 (so a VAL3 range of 2.9-6.9) is the 
% same as deactivating North and South, as it will give only a temperate
% zone.
%
% With the default values, VAL3 would range from 1.15 to 9.15.
% For cold maps use COLDF=-1.5 and HEATF=-8.
% For warm maps use COLDF=-14 and HEATF=-6.

ExtraFor=1; % ExtraFor is used for big maps to not have large areas 
% completely empty from forests unless the climate is dry.
Forechan = 0.5; % Forechan influences the prevelance of forests on plains.
% Default 0.5 Higher value = LESS forests. Values should be between -2 to 2

%% Initialize

Ny=Nx;
calcDM=1;
Ampli = 0.64; Presmooth=1; Mediter=0; 

if Sea>0
    if Mediter>0 && Sea>1; Sea=1; end
    if Sea>2; Sea=2; end
end

if mod(LOOP,10)~=1 % Recalculate the matrix for the axis of anisotropy every 10 loops 
    calcDM=0;
end

if Nx<40; Nx=40; end
if Ny<40; Ny=40; end
if corl<0.5; corl=0.5; end
if corl>5; corl=5; end

N=Nx*Ny;
theta=rand(1)*pi;
theta*180/pi

% Correlation length defines the "range" that two pixels are somewhat
% correlated (36.5% correlation using exponential covariance). 
% Larger values for Correlation length => less change with pixel distance => smaller scale (less rough change)

beta1=[1 16*corl 8*corl]; % [Variance  cor.length.x  cor.length.y] altitude; anisotropic. Use above [1 12 6] to see hills below mountains.
% [1 16 8] is good in my opinion for matrices that go about 3-4 miles per pixel. 
% [1 16*3 8*3] is what I use for matrices that go about 1-2 miles per pixel.  
beta2=[1 8*corl 8*corl];  % [Variance  cor.length.x  cor.length.y] Vegetation; Isotropic. Changes faster than altitude.
 Epan=0; 
 SHM(N,2)=0;
 

 
 
%% Matrices

if calcDM>0

Cexp=inline('expo(1)*(exp(-sqrt( (rx/expo(2)).^2 + (ry/expo(3)).^2 ) ))','expo','rx','ry');

for x=1:Nx
    SHM(Epan*Ny+1:(Epan+1)*Ny,1:2)=[x*ones(Ny,1),(1:Ny)'];
    Epan=Epan+1;
end

R=[cos(theta) -sin(theta); sin(theta) cos(theta)];
SHM=SHM*R;    
    
DX=pdist2(SHM(:,1),SHM(:,1));
DY=pdist2(SHM(:,2),SHM(:,2));

tic
C1=Cexp(beta1,DX,DY); A1=sqrtm(C1); %A1=chol(C1); % Use cholesky decomposition for speed, but it is worse.
clear C1
C2=Cexp(beta2,DX,DY); A2=sqrtm(C2); %A2=chol(C2); 
clear C2 DX DY
toc

save MatrDis SHM R A1 A2

 else
     load MatrDis
    
end

%% Simulation

for IMAG=1:SyntheM^2+1
    if mod(IMAG,10)==0
        clc
        disp('Part')
        IMAG
        toc(t0)
    end
    
load MatrDis

simn=0; goon=0;
t2=tic;
RFS=zeros(N,15);


while goon==0
    if mod(simn,10)==9 
        disp('Sim num')
        simn+1
    end
simn=simn+1;
u1=randn(N,1);
u2=randn(N,1);

if Sea>=2 & IMAG==1;
    u1(1:1*Nx)=1.3 + sqrt(Nx/beta1(2)/5);
    u1(N-1*Nx+1:end)=1.3 + sqrt(Nx/beta1(2)/5);
    u1(1:Ny:end)=1.3 + sqrt(Nx/beta1(2)/5);
    u1(Ny:Ny:end)=1.3 + sqrt(Nx/beta1(2)/5);
    
    Mount=floor(rand(1)*4+3*log(Nx)-4);
    for i=1:Mount
        xm=round(rand(1)*(Nx-8))+4;
        ym=round(rand(1)*(Ny-8))+4;
        u1(xm*Nx+ym-1:xm*Nx+ym)=-0.85*sqrt(Nx/10);
        u1((xm+1)*Nx+ym-1:(xm+1)*Nx+ym)=-0.85*sqrt(Nx/10);
    end
    
    if StSea>0
        ND=max([round(Nx/12.5),1]);
        u1=reshape(u1,Nx,Ny);
        u1(1:ND,1:ND)=1.3 + sqrt(Nx/beta1(2)/5);
        u1(1:ND,end-ND+1:end)=1.3 + sqrt(Nx/beta1(2)/5);
        u1(end-ND+1:end,end-ND+1:end)=1.3 + sqrt(Nx/beta1(2)/5);
        u1(end-ND+1:end,1:ND)=1.3 + sqrt(Nx/beta1(2)/5);
        
        u1(1,1)=u1(1,1)+1;
        u1(1,end)=u1(1,end)+1;
        u1(end,1)=u1(end,1)+1;
        u1(end,end)=u1(end,end)+1;
        
        u1=reshape(u1,N,1);
     
    end
    
end

if Sea==1 & IMAG==1;
    u1(1:1*Nx)=1.1 + sqrt(Nx/beta1(2)/8);
    u1(N-1*Nx+1:end)=1.1 + sqrt(Nx/beta1(2)/8);
    u1(1:Ny:end)=1.1 + sqrt(Nx/beta1(2)/8);
    u1(Ny:Ny:end)=1.1 + sqrt(Nx/beta1(2)/8);
    
    Mount=floor(rand(1)*4+3*log(Nx)-7);
    for i=1:Mount
        xm=round(rand(1)*(Nx-8))+4;
        ym=round(rand(1)*(Ny-8))+4;
        u1(xm*Nx+ym-1:xm*Nx+ym)=-0.85*sqrt(Nx/10);
        u1((xm+1)*Nx+ym-1:(xm+1)*Nx+ym)=-0.85*sqrt(Nx/10);
    end
    
        if StSea>0
        ND=round(Nx/12);
        
        if theta>pi/2+pi/8 & theta< pi-pi/8
            for i=1:ND
                u1( (Ny-ND+i-2)*Nx+2: (Ny-ND+i-2)*Nx+ND)=2.3;
            end
        end
        if theta<pi/2+pi/8 & theta> pi/8
            for i=1:ND
                u1(i*Nx+2:i*Nx+ND-i)=2;
            end
        end            
    end
    
end


if Mediter==1 & Nx>45 & Ny>45 & IMAG==1
        k=linspace(0.2,1.2,Ny);
        u1=reshape(u1,Nx,Ny);
         u1(1:3,1:Ny)=u1(1:3,1:Ny)+[k*2;k*1.5;k];
         u1(end-2:end,1:Ny)=u1(end-2:end,1:Ny)+[k;k*1.5;k*2];

        ND=max([round(Nx/18),1]);
        u1(1:ND,1:ND)= u1(1:ND,1:ND)+0.4;
        u1(1:ND,end-ND+1:end)= u1(1:ND,end-ND+1:end)+0.4;
        u1(end-ND+1:end,end-ND+1:end)= u1(end-ND+1:end,end-ND+1:end)+0.4;
        u1(end-ND+1:end,1:ND) = u1(end-ND+1:end,1:ND)+0.4;

        u1=reshape(u1,N,1);
        clear k
end

RF1=A1*u1;
RF2=A2*u2;

if simn<16
    RFS(:,simn)=RF1;
end

    DM=sum(RF1<-0.85);
    DW=sum(RF1>1.75);
    DPL=sum(RF1>=-0.35 & RF1<=1.75);
    
if Sea>0 
    RF1=A1*u1-1.5*beta1(2)/Nx;
    DM=sum(RF1<-0.65);
    DW=sum(RF1>1.95);
    DPL=sum(RF1>=-0.25 & RF1<=1.95);
end


    RatioM=DM/N;
    RatioW=DW/N;
    RatioML=DM/(N-DW);
    RatioPL=DPL/(N-DW);

if Sea<3    
    DCheck(1)=RatioM>0.06 -Mediter/50;
    DCheck(2)=RatioM<0.14 -Mediter/50; 
    DCheck(3)=RatioML<0.21;
    DCheck(4)=RatioW<0.35+Sea/20 +StSea/12 +Mediter/16;
    DCheck(5)=RatioPL>0.55; % Nearly half are lost as woods
end


if abs(skewness(RF1)<0.4+Sea/4 +StSea/2 +Mediter/2) & sum(DCheck)==5 & abs(mean(RF1))<3.5+Mediter
    goon=1;
end

% ********** Check to stop **************

if IMAG==1 & toc(t2)>(240+Mediter*360) | IMAG>1 & toc(t2)>(30+Mediter*5)
    goon=1;
    disp('error')
   
    DM=sum(RFS<-0.85);
    DW=sum(RFS>1.75);
    
    if Sea>0 & Sea<3
        DM=sum(RFS<-0.65);
        DW=sum(RFS>1.95);
    elseif Sea>=3
        DM=sum(RF1<0.4);
        DW=sum(RF1>1.6);
    end
    
    RatioM=DM/N;
    RatioW=DW/N;
    RatioML=DM/(N-DW);
 
if Sea<3    
    DCheck(1,1:15)=RatioM>0.06 -Mediter/50;
    DCheck(2,1:15)=RatioM<0.14 -Mediter/50; %na exei kai ligo plains... e3ou kai to rixnw me mediteran
    DCheck(3,1:15)=RatioML<0.21;
    DCheck(4,1:15)=RatioW<0.35+Sea/20+StSea/12+Mediter/12;
    DCheck(5,1:15)=RatioPL>0.45;
end
    
    DSUM=sum(DCheck);
    d=max(DSUM);
    DT=find(DSUM==d);
    DFinal=DT(1);
    RF1=RFS(:,DFinal);
    
    DCheck(:,DFinal)
    SK=skewness(RF1)
    DW=sum(RF1>1.95);
    clear DCheck
end

end
clear goon A1 A2 
simn

P1(:,IMAG)=RF1;
P2(:,IMAG)=RF2;

if IMAG==1
   save imag1 RF1 RF2  
end

end

%% Presmooth
% Smooth of first image
Gram=0;

V1T=P1(:,1); V2T=P2(:,1);
V1=reshape(V1T,Nx,Ny); V2=reshape(V2T,Nx,Ny);

h= [0.1250    0.2000    0.2500    0.2000    0.1250
    0.2000    0.5000    1.0000    0.5000    0.2000
    0.2500    1.0000    6.0000    1.0000    0.2500
    0.2000    0.5000    1.0000    0.5000    0.2000
    0.1250    0.2000    0.2500    0.2000    0.1250]/15.1;

V1T=filter2(h,V1); V2T=filter2(h,V2);
P1(:,1)=reshape(V1T,Nx*Ny,1); P2(:,1)=reshape(V2T,Nx*Ny,1);

%% Grids

P1(:,2:SyntheM^2+1)=P1(:,2:SyntheM^2+1)*Ampli;
P2(:,2:SyntheM^2+1)=P2(:,2:SyntheM^2+1)*Ampli; 

Gram=0;
 for i=1:SyntheM:SyntheM*Nx
    for j=1:SyntheM:SyntheM*Ny
        Gram=Gram+1;
        VALB1(i:i+SyntheM-1,j:j+SyntheM-1)=P1(Gram,1);
        VALB2(i:i+SyntheM-1,j:j+SyntheM-1)=P2(Gram,1);
    end
 end

Synthe=2;
for synx=1:SyntheM
for syny=1:SyntheM
 Gram=0; 
 for i=Nx*(synx-1)+1:Nx*synx
    for j=Ny*(syny-1)+1:Ny*syny
        Gram=Gram+1;
        VAL1(i,j)=P1(Gram,Synthe)+VALB1(i,j)*1;
        VAL2(i,j)=P2(Gram,Synthe)+VALB2(i,j)*1;
    end
 end
 Synthe=Synthe+1;
end
end


Nx=SyntheM*Nx; Ny=SyntheM*Ny; N=Nx*Ny;

%% Corrections
% These are arbitary corrections by me to make maps I like and smoothing.
% Some are ruled by values at the initialize part on the beginning.


% ********* adding more forests ********
if ExtraFor>0
    if Nx>120
    Akro=0.35; stpF=round(Nx/6); stpRan=floor(stpF/4)-2;
    for i=round(Akro*Nx)+1:stpF:round((1-Akro)*Nx)
        for j=1:stpF:Ny-stpF
            dfor=VAL2(i:i+stpF,j:j+stpF)>0.35;
            kfor=sum(sum(dfor));
            while kfor<(stpF^2)*0.2
                ranx=floor(rand(1)*stpRan*3);rany=floor(rand(1)*stpRan*3);
                
                VAL2(i+ranx+3:i+ranx+stpRan,j+rany+3:j+rany+stpRan)=VAL2(i+ranx+3:i+ranx+stpRan,j+rany+3:j+rany+stpRan)+rand(stpRan-2,stpRan-2)*0.3;
                               
                dfor=VAL2(i:i+stpF,j:j+stpF)>0.35;
                kfor=sum(sum(dfor))+30;
                clear dfor
            end
            clear kfor
        end
    end
    end
    
end


% ******** Rows ********
h=[1 -2 1]';
Fpluto=imfilter(VAL1,h);
d=abs(Fpluto)>1.2;
nerr=sum(sum(d));
K(nerr,1)=0;

Fpluto2=imfilter(VAL2,h);
d2=abs(Fpluto2)>1.2;
nerr2=sum(sum(d2));
K2(nerr2,1)=0;

[ro,co]=find(d==1);

for i=1:nerr
    if ro(i)>4 & ro(i)<Nx-4
        K(i)= ( VAL1(ro(i)-1,co(i))+VAL1(ro(i),co(i))+VAL1(ro(i)+1,co(i)) )/3;
        VAL1(ro(i)-2:ro(i)+2,co(i))=3*K(i)/4+VAL1(ro(i)-2:ro(i)+2,co(i))/4;
        VAL1(ro(i)-3,co(i))=K(i)/2 +VAL1(ro(i)-3,co(i))/2;
        VAL1(ro(i)+3,co(i))=K(i)/2 +VAL1(ro(i)+3,co(i))/2;
        VAL1(ro(i)-4,co(i))=K(i)/4 +3*VAL1(ro(i)-4,co(i))/2;
        VAL1(ro(i)+4,co(i))=K(i)/4 +3*VAL1(ro(i)+4,co(i))/2;
    end
end
VAL1(d)=K;

[ro,co]=find(d2==1);
for i=1:nerr2
    if ro(i)>4 & ro(i)<Nx-4
        K2(i)= ( VAL2(ro(i)-1,co(i))+VAL2(ro(i),co(i))+VAL2(ro(i)+1,co(i)) )/3;
        VAL2(ro(i)-2:ro(i)+2,co(i))=3*K2(i)/4+VAL2(ro(i)-2:ro(i)+2,co(i))/4;
        VAL2(ro(i)-3,co(i))=K2(i)/2 +VAL2(ro(i)-3,co(i))/2;
        VAL2(ro(i)+3,co(i))=K2(i)/2 +VAL2(ro(i)+3,co(i))/2;
        VAL2(ro(i)-4,co(i))=K2(i)/4 +3*VAL2(ro(i)-4,co(i))/4;
        VAL2(ro(i)+4,co(i))=K2(i)/4 +3*VAL2(ro(i)+4,co(i))/4;
    end
end
VAL2(d2)=K2;

clear K d K2 d2 Fpluto Fpluto2 

% ******** Columns ********
h=[1 -2 1];
Fpluto=imfilter(VAL1,h);
d=abs(Fpluto)>1.2;
nerr=sum(sum(d));
K(nerr,1)=0; 

Fpluto2=imfilter(VAL2,h);
d2=abs(Fpluto2)>1.2;
nerr2=sum(sum(d2));
K2(nerr2,1)=0;

[ro,co]=find(d==1);

for i=1:nerr
    if co(i)>4 & co(i)<Nx-4
        K(i)= ( VAL1(ro(i),co(i)-1)+VAL1(ro(i),co(i))+VAL1(ro(i),co(i)+1) )/3;
        VAL1(ro(i),co(i)-2:co(i)+2)=K(i)/4+3*VAL1(ro(i),co(i)-2:co(i)+2)/4;
        VAL1(ro(i),co(i)-3)=K(i)/2 +VAL1(ro(i),co(i)-3)/2;
        VAL1(ro(i),co(i)+3)=K(i)/2 +VAL1(ro(i),co(i)+3)/2;
        VAL1(ro(i),co(i)-4)=K(i)/4 +3*VAL1(ro(i),co(i)-4)/4;
        VAL1(ro(i),co(i)+4)=K(i)/4 +3*VAL1(ro(i),co(i)+4)/4;
    end
end
VAL1(d)=K;

[ro,co]=find(d2==1);
for i=1:nerr2
    if co(i)>4 & co(i)<Nx-4
        K2(i)= ( VAL2(ro(i),co(i)-1)+VAL2(ro(i),co(i))+VAL2(ro(i),co(i)+1) )/3;
        VAL2(ro(i),co(i)-2:co(i)+2)=K2(i)/4+3*VAL2(ro(i),co(i)-2:co(i)+2)/4;
        VAL2(ro(i),co(i)-3)=K2(i)/2 +VAL2(ro(i),co(i)-3)/2;
        VAL2(ro(i),co(i)+3)=K2(i)/2 +VAL2(ro(i),co(i)+3)/2;
        VAL2(ro(i),co(i)-4)=K2(i)/4 +3*VAL2(ro(i),co(i)-4)/4;
        VAL2(ro(i),co(i)+4)=K2(i)/4 +3*VAL2(ro(i),co(i)+4)/4;
    end
end
VAL2(d2)=K2;

clear K d K2 d2 Fpluto Fpluto2

% ********* adding sea ********
if Sea>0
    VAL1(1,:)=3;
    VAL1(Nx,:)=3;
    VAL1(:,1)=3;
    VAL1(:,Ny)=3;
    
    VAL1(2:4,:)=VAL1(2:4,:)+0.7;
    VAL1(Nx-4:Nx-2,:)=VAL1(2:4,:)+0.7;
    VAL1(:,2:4)=VAL1(:,2:4)+0.7;
    VAL1(:,Ny-2:Ny-4)=VAL1(:,Ny-2:Ny-4)+0.7;
    if Sea==3
        VAL1(2:4,:)=VAL1(2:4,:)+0.7;
        VAL1(Nx-4:Nx-2,:)=VAL1(2:4,:)+0.7;
        VAL1(:,2:4)=VAL1(:,2:4)+0.7;
        VAL1(:,Ny-2:Ny-4)=VAL1(:,Ny-2:Ny-4)+0.7;
    end
    
end

 
%% Values

yk=1:Nx; xk=1:Ny;

figure(1)

axes1 = axes('Parent',figure(1),'DataAspectRatio',[1 1 1],'FontSize',14);

grid(axes1,'on');
hold(axes1,'all');
surf(xk,yk,VAL1)
shading flat
colorbar
colormap('bone')
view(2)

figure(2)
axes1 = axes('Parent',figure(2),'DataAspectRatio',[1 1 1],'FontSize',14);

grid(axes1,'on');
hold(axes1,'all');
surf(xk,yk,VAL2)
shading flat
colorbar
colormap('bone')
view(2)

figure(5)

axes1 = axes('Parent',figure(5),'DataAspectRatio',[1 1 1],'FontSize',14);

grid(axes1,'on');
hold(axes1,'all');
surf(xk,yk,VALB1)
shading flat
colorbar
colormap('bone')
view(2)

figure(6)
axes1 = axes('Parent',figure(6),'DataAspectRatio',[1 1 1],'FontSize',14);

grid(axes1,'on');
hold(axes1,'all');
surf(xk,yk,VALB2)
shading flat
colorbar
colormap('bone')
view(2)

%% TRUNCTUATION 
% Here the values become terrains

%% Trunctuation no sea

VAL(Nx,Ny)=0;
Deik1=VAL1<-0.85;
VAL(Deik1==1)=1; % Mountain

Deik1=VAL1>=-0.85 & VAL1<-0.35;
VAL(Deik1==1)=2; % Hill

Deik1=VAL1>=-0.35 & VAL1<=1.75;
VAL(Deik1==1)=3; % Plains

Deik1=VAL1>1.75;
VAL(Deik1==1)=4; % Water

Deik1=VAL1>=0.35 & VAL1<=1.75;
Deik2=VAL2>0;
VAL(Deik1==1 & Deik2==1)=5; 
Deik1=VAL1>=-0.25 & VAL1<0.35;
Deik2=VAL2>0.5;
VAL(Deik1==1 & Deik2==1)=5; % Forest

Deik1=VAL1>=1 & VAL1<=1.75;
Deik2=VAL2>0.5;
VAL(Deik1==1 & Deik2==1)=6; 
Deik1=VAL1>=-0.15 & VAL1<1;
Deik2=VAL2>1;
VAL(Deik1==1 & Deik2==1)=6; % Dense Forest

Deik1=VAL1>=-0.65 & VAL1<-0.25;
Deik2=VAL2>0.5;
VAL(Deik1==1 & Deik2==1)=7; % Forested Hills

Deik1=VAL1>=1.25 & VAL1<=1.75;
Deik2=VAL2>-0.25 & VAL2<0.75;
VAL(Deik1==1 & Deik2==1)=8; % Swamp



%% Trunctuation with sea

if Sea>0 
    
VAL(Nx,Ny)=0;
% RF1
Deik1=VAL1<-0.65;
VAL(Deik1==1)=1; 

Deik1=VAL1>=-0.65 & VAL1<-0.25;
VAL(Deik1==1)=2; 

Deik1=VAL1>=-0.25 & VAL1<=1.95;
VAL(Deik1==1)=3; 

Deik1=VAL1>1.95;
VAL(Deik1==1)=4;

Deik1=VAL1>=0.45 & VAL1<=1.95;
Deik2=VAL2>0;
VAL(Deik1==1 & Deik2==1)=5; 
Deik1=VAL1>=-0.15 & VAL1<0.45;
Deik2=VAL2>0.5;
VAL(Deik1==1 & Deik2==1)=5; 

Deik1=VAL1>=1.1 & VAL1<=1.95;
Deik2=VAL2>0.5;
VAL(Deik1==1 & Deik2==1)=6; 
Deik1=VAL1>=-0.10 & VAL1<1.1;
Deik2=VAL2>1;
VAL(Deik1==1 & Deik2==1)=6; 

Deik1=VAL1>=-0.55 & VAL1<-0.20;
Deik2=VAL2>0.5;
VAL(Deik1==1 & Deik2==1)=7; 

Deik1=VAL1>=1.65 & VAL1<=1.90;
Deik2=VAL2>-0.10 & VAL2<0.60;
VAL(Deik1==1 & Deik2==1)=8; 
    
    
end


%%  Forest mask
% Converts some plains to forests, DIRECTLY, not using the simulated values

DCH(Nx,Ny)=0;

DFFull=VAL==5 | VAL==6 | VAL==7;
DFHalf=VAL==8 | VAL==4 ; % Sto nero half gia na mhn metra san 0

DF=DFFull+1/2*DFHalf;

FMASK=[0 1 1.5 1 0; 1 2 3 2 1; 1.5 3 -34 3 1.5; 1 2 3 2 1; 0 1 1.5 1 0];

FH=filter2(FMASK,DF,'valid');
DF=FH>19; DRF=FH>15 & FH<=19; 
DPL=VAL==3; DRAN=randn(Nx-4,Ny-4)>Forechan;
DF2=DRF==1 & DRAN==1;
DF=DF+DF2;

DCH(3:Nx-2,3:Ny-2)=DPL(3:Nx-2,3:Ny-2)==1 & DF==1;

VAL(DCH==1)=5;

clear FMASK DF FH DFFull DFHalf DRF DPL DRAN DF2 DCH



%% Rotate

if imrot==1

AVAL=reshape(VAL',Nx*Ny,1);
AVAL1=reshape(VAL1',Nx*Ny,1);
AVAL2=reshape(VAL2',Nx*Ny,1);
Epan=0;
for x=1:Nx
    SHMT(Epan*Ny+1:(Epan+1)*Ny,1:2)=[x*ones(Ny,1),(1:Ny)'];
    Epan=Epan+1;
end

theta2=1/12*pi+rand(1)*4/12*pi;
R2=[cos(theta2) -sin(theta2); sin(theta2) cos(theta2)];
SHMT=SHMT*R2;  
SHMT(:,1)=SHMT(:,1)+abs(min(SHMT(:,1)));
SHMT(:,2)=SHMT(:,2)+abs(min(SHMT(:,2)));
A=SHMT;

Nx=round(Nx*sqrt(2)); Ny=round(Ny*sqrt(2)); 
Epan=0;
for x=1:Nx
    SHMT(Epan*Ny+1:(Epan+1)*Ny,1:2)=[x*ones(Ny,1),(1:Ny)'];
    Epan=Epan+1;
end

t=knnsearch(A,SHMT);
VAL=AVAL(t);
VAL1=AVAL1(t);
VAL2=AVAL2(t);
VAL=reshape(VAL,Nx,Ny)';
VAL1=reshape(VAL1,Nx,Ny)';
VAL2=reshape(VAL2,Nx,Ny)';


clear A SHMT t
end

%% Cold and heat

% Trend
for i=1:Ny
    VAL3(1:Nx,i)=linspace(0,(10+HEATF),Nx)'-COLDF/2;
end
% cold climates start at VAL3<2.9 and warm climates start at VAL3>6.9

    K=linspace(0,-(10+HEATF)/15,floor(Ny/2));
    K2=linspace(-(10+HEATF)/15,0,Ny-floor(Ny/2));
    
for i=1:round(Nx/2)-5
    VAL3(i,1:Ny)=VAL3(i,1:Ny)+[K K2];    
end
for i=round(Nx/2)+5:Nx
    VAL3(i,1:Ny)=VAL3(i,1:Ny)-[K K2];    
end

Tstp=(10+HEATF)/Nx;
kk=max([floor(Nx/90), 1]);
for i=1:kk:Nx-kk+1
    for j=1:kk:Ny-kk+1
        VAL3(i:i+kk-1,j:j+kk-1)=VAL3(i:i+kk-1,j:j+kk-1)+randn(1)*ones(kk,kk)*Tstp*0.99;
    end
end

clear K K2 Tstp

%% ******* Climate ********

% ******* Cold *******
DeikCpp=VAL3<2.9; 
DeikCp=VAL3<2.5;
DeikC=VAL3<2;
DeikCm=VAL3<1;

if North==1
    
% Forest & Cold
Deik1=VAL==5;
Deik2=VAL2>0 & VAL2<0.6;
VAL(Deik1==1 & DeikCpp==1 & Deik2==1)=3; % 3= Plains 
VAL(Deik1==1 & DeikC==1)=9; % 9= Tundra 


% Deep forest & Cold
Deik1=VAL==6;
VAL(Deik1==1 & DeikC==1)=10; % 10= Taiga

% Plains & Cold
Deik1=VAL==3;
VAL(Deik1==1 & DeikC==1)=9; % 9= Tundra

% Hills & Cold+
Deik1=VAL==2;
VAL(Deik1==1 & DeikCp==1)=11; % 11= Snow hills

% Mountains & Cold++
Deik1=VAL==1;
VAL(Deik1==1 & DeikCpp==1)=12; % 12= Snow mountains

% Forest Hills & Cold-
Deik1=VAL==7;
VAL(Deik1==1 & DeikCpp==1)=2; % 2= Hills
VAL(Deik1==1 & DeikCp==1)=11; % 11= Snow hills

% Water & Cold-
Deik1=VAL==4;
VAL(Deik1==1 & DeikCm==1)=20; % 20 = Glacier

end

% ******** warm *******

DeikHpp=VAL3>6.9;
DeikHp=VAL3>7.45; 
DeikH=VAL3>8;
DeikHm=VAL3>8.65; 

VAL2(DeikHp==1)=VAL2(DeikHp==1)-0.15;
VAL2(DeikH==1)=VAL2(DeikH==1)-0.2;
VAL2(DeikHm==1)=VAL2(DeikHm==1)-0.25;

if South==1
    
% Deep forest & Heat

Deik1=VAL==6;
VAL(Deik1==1 & DeikHm==1)=13; % 13= Jungle

% forest & Heat
Deik1=VAL==5;
Deik2=VAL2<0;
VAL(Deik1==1 & Deik2==1 & DeikH==1)=3; % 

% Plains & Heat

Deik1=VAL==3;
Deik2=VAL2<-1.4;
Deik3=VAL2<-1.7;
Deik4=VAL2<-2.1;
Deik5=VAL2<-2.7;
VAL(Deik1==1 & DeikHpp==1 & Deik4==1)=14; % 14= semi-arid
VAL(Deik1==1 & DeikHp==1 & Deik2==1)=14; % 14= semi-arid
VAL(Deik1==1 & DeikH==1 & Deik3==1)=15; % 15= desert
VAL(Deik1==1 & DeikHpp==1 & Deik5==1)=15;
end

% ******** Different Plains *********
if DiffPlain==1
    
    Deik1=VAL==3 | VAL==14;
    Deik2=VAL2<-0.4 & VAL2> -1.2;
    D=Deik1+Deik2;
    
    h=ones(7,7); h(4,4)=50; h=h/98; 
    K=imfilter(D,h)/2;
    dK=K>85/98;
    
    VAL(dK==1 & DeikHp==1) =18; % Savanna
    
    Deik2=VAL2<-0.7 & VAL2> -2;
    D=Deik1+Deik2;
    K=imfilter(D,h)/2;
    dK=K>80/98;
    VAL(dK==1 & VAL3<4 & DeikC==0) =19; % Steppe
    clear K K2 Deik1 Deik2 dK h
    
end

%% Hexographer Compact image
figure(3)

IM(Nx,Ny,3)=0;

DC=VAL==1;
IM(:,:,1)=IM(:,:,1)+DC*178;
IM(:,:,2)=IM(:,:,2)+DC*128;
IM(:,:,3)=IM(:,:,3)+DC*0;

DC=VAL==2;
IM(:,:,1)=IM(:,:,1)+DC*231;
IM(:,:,2)=IM(:,:,2)+DC*206;
IM(:,:,3)=IM(:,:,3)+DC*90;

DC=VAL==3;
IM(:,:,1)=IM(:,:,1)+DC*160;
IM(:,:,2)=IM(:,:,2)+DC*215;
IM(:,:,3)=IM(:,:,3)+DC*107;

DC=VAL==4;
IM(:,:,1)=IM(:,:,1)+DC*153;
IM(:,:,2)=IM(:,:,2)+DC*204;
IM(:,:,3)=IM(:,:,3)+DC*255;

DC=VAL==5;
IM(:,:,1)=IM(:,:,1)+DC*57;
IM(:,:,2)=IM(:,:,2)+DC*145;
IM(:,:,3)=IM(:,:,3)+DC*45;

DC=VAL==6;
IM(:,:,1)=IM(:,:,1)+DC*47;
IM(:,:,2)=IM(:,:,2)+DC*118;
IM(:,:,3)=IM(:,:,3)+DC*33;

DC=VAL==7;
IM(:,:,1)=IM(:,:,1)+DC*72;
IM(:,:,2)=IM(:,:,2)+DC*138;
IM(:,:,3)=IM(:,:,3)+DC*33;

DC=VAL==8;
IM(:,:,1)=IM(:,:,1)+DC*173;
IM(:,:,2)=IM(:,:,2)+DC*222;
IM(:,:,3)=IM(:,:,3)+DC*165;

DC=VAL==19;
IM(:,:,1)=IM(:,:,1)+DC*190;
IM(:,:,2)=IM(:,:,2)+DC*190;
IM(:,:,3)=IM(:,:,3)+DC*190;

if North==1
DC=VAL==9;
IM(:,:,1)=IM(:,:,1)+DC*200;
IM(:,:,2)=IM(:,:,2)+DC*230;
IM(:,:,3)=IM(:,:,3)+DC*195;

DC=VAL==10;
IM(:,:,1)=IM(:,:,1)+DC*160;
IM(:,:,2)=IM(:,:,2)+DC*215;
IM(:,:,3)=IM(:,:,3)+DC*155;

DC=VAL==11;
IM(:,:,1)=IM(:,:,1)+DC*230;
IM(:,:,2)=IM(:,:,2)+DC*245;
IM(:,:,3)=IM(:,:,3)+DC*215;

DC=VAL==12;
IM(:,:,1)=IM(:,:,1)+DC*250;
IM(:,:,2)=IM(:,:,2)+DC*235;
IM(:,:,3)=IM(:,:,3)+DC*225;

DC=VAL==20;
IM(:,:,1)=IM(:,:,1)+DC*250;
IM(:,:,2)=IM(:,:,2)+DC*250;
IM(:,:,3)=IM(:,:,3)+DC*250;
end


if South==1
DC=VAL==13;
IM(:,:,1)=IM(:,:,1)+DC*77;
IM(:,:,2)=IM(:,:,2)+DC*143;
IM(:,:,3)=IM(:,:,3)+DC*90;

DC=VAL==14;
IM(:,:,1)=IM(:,:,1)+DC*180;
IM(:,:,2)=IM(:,:,2)+DC*240;
IM(:,:,3)=IM(:,:,3)+DC*70;

DC=VAL==15;
IM(:,:,1)=IM(:,:,1)+DC*248;
IM(:,:,2)=IM(:,:,2)+DC*240;
IM(:,:,3)=IM(:,:,3)+DC*133;

DC=VAL==18;
IM(:,:,1)=IM(:,:,1)+DC*210;
IM(:,:,2)=IM(:,:,2)+DC*248;
IM(:,:,3)=IM(:,:,3)+DC*150;


end

IM=IM/255;

imshow(IM);

D=VAL==3 | VAL==18 | VAL==19; DW=sum(sum(VAL==4 | VAL==20)); N=Nx*Ny;
RatioPL2=round(sum(sum(D))/(N-DW)*100);

LO=num2str(LOOP); PLR=num2str(RatioPL2); DWL=num2str(N-DW);

TIT=strcat('figs\IM',LO,'PR',PLR,'L',DWL);

if Nx<700
    print(3,'-dpng','-r150',TIT) 
else
    print(3,'-dpng','-r300',TIT) 
end
%export_fig(3,'-pdf','-transparent',TIT)

SynM=SyntheM;
TIT2=strcat('PGSmatr\IM',LO,'PR',PLR,'L',DWL);
save(TIT2, 'VAL1', 'VAL2', 'VAL3', 'VAL', 'COLDF', 'HEATF', 'SynM');

clear VAL IM VAL2 VAL3 VAL1 


end

toc(t0)



