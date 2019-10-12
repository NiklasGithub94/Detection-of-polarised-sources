%% MyProgram
clear all
fname = 'FarFarAway_p.fits';
fifi = fitsread(fname);
fifi(isnan(fifi(:,:,:)))=0;
fname1='FarFaraway.fits';
fif1=fitsread(fname1);
%NEED TO FIND SIZE OF FARADAYCUBE AND POSITION IN SKY TO CORRELATE TO PIXELS

%% IDENTIFYING AND FINDING THE CENTROIDS OF BLOBS/GALAXIES
clc
%PHI = SOME NUMBER
% mesh(fifi(:,:,PHI))
% colorbar
clc

X=fif1; %fif1; %EX: Initialisation 
J = imadjust(X, [0.009 0.01], []); %imadjust(X, [x, y], []); , [0.009,0.01]
% x increases darkness in map and decides about numbers of sources, y increases brigthness 
%J=imadjust(JLite, [0.01,0.3],[]); - OPTIONAL: SECOND FILTER
%  figure(11)
%  imshow(JLite)
% figure(12)
% imshow(J)
BW=J;
figure(1)
imshowpair(J,fif1,'montage')
 %BW = fifi(:,:,10); % FARADAYCUBE OR INTENSITYMAP?
[ii,jj]=find(BW); %FINDS COORDS FOR BLOBS, ii=ycoords, jj=xcoords

CC = bwconncomp(BW); %connects neighbour pixels, identifies them as blobs
s = regionprops(CC,'centroid'); %IMPORTANT STUFF!!
centroids = cat(1, s.Centroid); %IMPORTANT STUFF!!
figure(2)
plot(jj,ii,'.') %MAP

%PLOTS EVERY SOURCE AND ITS CENTROID
stats = regionprops('table',CC,'Centroid',... %ALL RELEVANT INFO, 
   'Area') %SELECT SOURCE, HOW?
Area=table2array(stats(:,1));
Centroid=table2array(stats(:,2));
xcoord=Centroid(:,1);
ycoord=Centroid(:,2);
sprintf('Number of blobs are: ')
NumberOfBlobs=size(Area,1)
hold on
title('Map with centroids')
xlabel('x - axis') 
ylabel('y - axis') 
plot(round(xcoord),round(ycoord),'*') %ROUNDED COORDINATES
%% FIND POLARISATION (COORDs/PIXELs) BASED ON THRESHOLD THEN PRODUCE F(PHI)
clc

x = 1:size(fifi,3);
M=zeros(length(x),sum(Area));
count = 0;
vdistr=zeros(size(fifi,2), size(fifi,1));
vmustd=zeros(size(fifi,2), size(fifi,1));
% img=fifi(:,:,1);
% [xgrid, ygrid] = meshgrid(1:size(img,2), 1:size(img,1));
% values = img(mask)

%HOW TO PRINT # IN MATLAB? 
%fprintf(" # Region file format: DS9 version 4.1 global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
%image")
 for i=1:sum(Area)

   for k=1:length(x) %SUM OVER PHI
       %if(jj(i) > 1 & ii(i)> 1 & jj(i)< 269 & ii(i) < 266)
           
       %M(k,i) = fifi(ii(i),jj(i),x(k)) + fifi(ii(i)+1,jj(i),x(k))+fifi(ii(i),jj(i)+1,x(k))+ fifi(ii(i)-1,jj(i),x(k))+ fifi(ii(i),jj(i)-1,x(k))+ fifi(ii(i)+1,jj(i)-1,x(k))+fifi(ii(i)-1,jj(i)+1,x(k))+fifi(ii(i)-1,jj(i)-1,x(k))+ fifi(ii(i)+1,jj(i)+1,x(k)); 
       
       %M(k,i) = M(k,i)/9;
       %M(k,i) = fifi(ii(i),jj(i),k); - one pixel detection 
       M(k,i) = conv2(fifi(ii(i),jj(i),k), ones(5,5), 'same')./conv2(ones(size(fifi(ii(i),jj(i),k))), ones(5,5), 'same');
      %M(k,i) = medfilt3(fifi(ii(i),jj(i),k), [5 5 1]); - VERY SLOW!!
       % M(k,i)= conv2(fifi(ii(i),jj(i),k), ones(3)/9, 'valid');
       %end

   end
   for k=388:401 %REMOVE INSTRUMENTAL POLARISATION
       if ( M(k,i) > (mean(M(:,i))+3*std(M(:,i))) )
       M(k,i) = mean(M(:,i)); 
       end
   end
  
   for k=1:length(x)
        if (M(k,i) > (mean(M(:,i)) + 8*std(M(:,i))))
         fprintf("circle(" + jj(i) + "," + ii(i) + "," +2.222+ ")" + newline)
         %fprintf("(" + jj(i)+ "," + ii(i) + "," + "phi = " +(x(k)-400)*0.3 + ")" + newline)
            count = count + 1;
         PolarisedSources(1,i)=jj(i);
         PolarisedSources(2,i)=ii(i);
     vdistr(ii(i),jj(i))=(x(k)-400)*0.3;
     %vmustd(ii(i),jj(i)) = mean(M(:,i))/std(M(:,i));
    %     txt = num2str(x(k)-100);
   %text(jj(i),ii(i),txt,'FontSize',7);
   
   %fifi(jj(i),ii(i),x(k));
   hold on
   
        end
     
    end
  
  
 end
fprintf("Number of found polarised sources are: " + count + newline)
% DERIVATIVE/VECTOR FIELD OF RM
% Different Colormaps in one figure
%[DX,DY] = gradient(vdistr);
ax1=axes;
mesh(ax1,fif1)
view(2)
ax2=axes;
contour(ax2,vdistr)
linkaxes([ax1 ax2])
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
cmap=colormap(ax1,gray); %jet - very strong intensity
newmap=imadjust(cmap, [0.005 0.15]);
colormap(ax1,newmap)
colormap(ax2,'Jet')
%set([ax1,ax2],'Position',[.17 .11 .685 .815]);
axis auto
cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);
%hold on
%quiver(DX,DY)
%hold off

% REVERSES THE X-AXIS AND KEEPS THE IMAGE THE SAME
% figure(3)
% hold on 
% im = image('XData', [25 45], 'YData', [45 50], 'CData', fif1, 'CDataMapping', 'scaled');
% im.CData = fliplr(im.CData) 
% colormap(newmap)
% colormap(ax2,'Jet')
% set(gca,'xdir','reverse')
% AxesH = gca;
% axis(AxesH, 'tight');
% hold off



%CDELT NOT TRUE
%YdifferenceFinal=0.00125;
%XdifferenceFinal; %DELTA1/COS(DELTA2); 
% fitsdisp(fname1,'Mode','full');
% info = fitsinfo('M51_0_low_I_corr.fits');
% information=info.PrimaryData;
% Cell=information.Keywords; %CELL WITH ALL INFORMATION FROM THE HEADER
%THE FOLLOWING PIXELS ARE CORRESPONDING TO THE FOLLOWING DEGREES
% format long
% YdifferenceFinal = 0.00125;
% XdifferenceFinal=0;
% CrPix1 = Cell{find(strcmp(Cell, 'CRPIX1')),2};
% CrPix2 = Cell{find(strcmp(Cell, 'CRPIX2')),2};
% CrVal1 = Cell{find(strcmp(Cell, 'CRVAL1')),2};
% CrVal2 = Cell{find(strcmp(Cell, 'CRVAL2')),2};
% vdelta=zeros(size(PolarisedSources,1),size(PolarisedSources,2));
% vDistr=zeros(size(PolarisedSources,1),size(PolarisedSources,2));
% vdelta(:,1)=CrPix1-PolarisedSources(:,1);
% vdelta(:,2)=CrPix2-PolarisedSources(:,2);
% vDistr(:,1)=CrVal1-XdifferenceFinal*vdelta(:,1) ;   
% vDistr(:,2)=CrVal2-YdifferenceFinal*vdelta(:,2);
% vDistr
%% PLOT FARADAY SPECTRUM
for w = 1:length(x)
    O(w) = (fifi(325,70,x(w))); %syntax fifi(y,x,phi)
end
[val,index]=max(O);
figure(2)
bar((x-400)*0.3,O)
for w = 1:length(x)
    O(w) = (fifi(308,45,x(w)));
end
figure(3)
bar((x-400)*0.3,O)
%% CHECK DISTRIBUTION OF PHI - WHICH IS THE MOST COMMON?
%M=mode(vdistr);

nbins=count;
h=histogram(vdistr);
Nbins = morebins(h);
h.Normalization = 'countdensity';
counts=h.Values;
%% MACHINE LEARNING - CLASSIFICATION TASK nnstart
clc
N = length(x); 
limit1=1.5e-4;
limit2=1e-5;
a = 0.3;
b = 0.9;
r = (b-a).*rand(N,1) + a;
LargeMatrix=zeros(N,N*100); %INTIALISING FAKE FARADAY SPECTRA MATRIX
LargeMatrixTarget=zeros(2,N*100); %INTIALISING TARGETS/RESULTS MATRIX

for i=1:(size(LargeMatrix,2)/2) %CREATING FAKE NON-ACCEPTABLE AND TARGETS
    vmach = limit2 + (limit1-limit2).*r;
    LargeMatrix(:,i) = vmach;
    v = linspace(0.1,3.9,500);
    p = randi(length(v));
    random = mean(vmach)+v(p)*std(vmach);
    phi_value = randi(N);
    LargeMatrix(phi_value,i)=random;
    LargeMatrixTarget(1,i)=0; 
    LargeMatrixTarget(2,i)=1;  

end

for i = (size(LargeMatrix,2)/2+1) : size(LargeMatrix,2) % FAKE ACCEPTABLE FARADAYSPECTRA AND TARGETS
   vmach = limit2 + (limit1-limit2).*r;
   LargeMatrix(:,i) = vmach;
   phi_value=randi(N);
   v=4;
   p=randi(length(v));
   random = mean(vmach)+v(p)*std(vmach);
   LargeMatrix(phi_value,i)=random;
   LargeMatrixTarget(1,i)=1; 
   LargeMatrixTarget(2,i)=0;  
end

inputs = LargeMatrix;
targets = LargeMatrixTarget;
%% RUN NETWORK

% Create a Pattern Recognition Network
hiddenLayerSize = [20 20 20]; %20 20 20
net = patternnet(hiddenLayerSize);


% Set up Division of Data for Training, Validation, Testing
net.divideParam.trainRatio = 70/100;
net.divideParam.valRatio = 15/100;
net.divideParam.testRatio = 15/100;


% Train the Network
[net,tr] = train(net,inputs,targets);

% Test the Network
outputs = net(inputs);
errors = gsubtract(targets,outputs);
performance = perform(net,targets,outputs);
% View the Network
%view(net)
figure, plotconfusion(targets,outputs)

%%
%SIMILAR TO USED CODE:
%       NEIGHBOUR = 5;
%  kernel = ones(NEIGHBOUR);
%  kernel = kernel / sum(kernel(:));
%  M(k,i) = conv2(fifi(ii(i),jj(i),k), kernel, 'same');
%AND
% avg = conv2(A, ones(NEIGHBOUR)/9, 'valid');
%AND
% k=[1 1 1; 1 0 1; 1 1 1]/8;
% averageIntensities = conv2(double(yourImage),k,'same');
% NOT PREFERABLE:
%  meanFilter = fspecial('average', [4 4]);
 % toShow = imfilter(A, meanFilter);
 
%figure(3)
% O=zeros(length(x),1);
% for w = 1:length(x) 
%     O(w) = (fifi(ii(7350),jj(7350),x(w))); %syntax fifi(y,x,phi)
% end
% bar(x-100,O)
% xlabel('\phi')
% hold on
% [TF,S1] = ischange(M(:,7350));
% stairs(x-101,S1);