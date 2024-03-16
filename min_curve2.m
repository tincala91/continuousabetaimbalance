function [summary,R2]=min_curve2(b, CSF, PET);

% Estimates individual aggregation and severity scores, 
% as well as a pseudo-R2 for the hyperbolic regression model 
% obtained with the function mincurve1
% 
% FORMAT [summary,R2]=min_curve2(b, CSF, PET);
% 
%         b      - 3x1 array containing the optimized parameters for the modified hyperbolic regression model
%                  obtained from function mincurve
%                  
%         CSF    - 1xN array containing the observed inidividual CSF values
%         
%         PET    - 1xN arrat containing the observed individual PET uptake values 
%         
%         summary- Nx9 matrix containing observed CSF values, observed PET values, 
%                  minimum Euclidean distance from curve (=Aggregation), 
%                  estimated x, estimated y, estimated distance from median
%                  (=Severity)
%
%         R2     - 1x1 arrat containing pseudo-R2
%         
% %--------------------------------------------------------------
% This function derives an Aβ-aggregation score  that reflects the imbalance 
% between soluble and aggregated Aβ, represented by the standardized Euclidean 
% distance of each observed value oij from its predicted value ôij. Here, a 
% positive Aβ-aggregation score denotes a higher PET CL value than expected 
% for a given value of CSF-Aβ42 (i.e., more aggregated relative to soluble Aβ),
% whereas a negative score denotes lower CSF-Aβ42 values than expected for 
% a given value of PET CL (i.e., more soluble relative to aggregated Aβ). 
% In addition, we estimated an Aβ-severity score, represented by the 
% standardized Euclidean distance between each data point’s corresponding 
% predicted value ôij, and the median of all predicted values. A positive 
% Aβ-severity score reflects CSF-Aβ42 and PET CL values that are more severe 
% than the median, while a negative Aβ-severity score reflects CSF-Aβ42 and 
% PET CL values that are less severe than the median. The function also
% derives a pseudo-R2 for model fit.
%
% Inputs: optimized parameters for the hyperbolic model, obtained from mincurve function; 
% observed values for both CSF and PET measures.
%
% Outputs: pseudo-R2, individual aggregation and severity Ab scores
% --------------------------------------------------------------
% Authors: 
% Juan Doming Gispert, BBRC-Fundació Pasqual Maragall
% Arianna Sala, Karolinska Institutet
% Date: 15.03.2024
%
% The code is provided under GPLv3 license 
%
% Any publication based on this code should cite:
% Mastenbroek, Sala et al., A continuous amyloid-β CSF/PET imbalance model
% to capture Alzheimer’s disease heterogeneity (2024). Neurology
% 
% --------------------------------------------------------------
% 
% Ref: 
% Mastenbroek, Sala et al., A continuous amyloid-β CSF/PET imbalance model
% to capture Alzheimer’s disease heterogeneity (2024). Neurology
%
%---------------------------------------------------------------

%A scaling factor is computed so as the balance the weight of each variable on the
%model estimation
ScalingFactor = std(CSF)/std(PET); % To weight x-distance = y-distance

%Set number of points of the fitted curve to be estimated
Npoints = 10000;

%PET is rescaled to match the range of CSF values
PET=PET*ScalingFactor;

%plot data points
%figure; scatter(PET,CSF); hold on

%generate vector of 10000 points (equally spaced) from min to max PET
%value
x=linspace(min(PET),max(PET),Npoints); 
%report here final b coefficient values selected in the curve_estimation
%script
%b = [3673.1 -2924.9 -0.3]; %declares beta coefficients - ref correct
y= b(1)+((b(2).*x./ScalingFactor)./((x/ScalingFactor)-b(3))); %compute y corresponding to x based on curve formula 
tol=0.001;
idx = find(x>(b(3)+tol)*ScalingFactor); %index to select only values in the part of the curve in the upper right quadrant, given condition: x bigger than asymptote
x=x(idx); %select only part of the curve based on condition defined above
y=y(idx); %select only part of the curve based on condition defined above

curvexy = [x;y]'; %store coordinates of the curve

%for all observations, stores distance of point from the curve, and point's corresponding
%estimate on the curve (predicted value)
for i = 1:length(CSF)
  [res(i),est_x(i), est_y(i)]=min_dist_point_curve([PET(i) CSF(i)],curvexy);
end

%computing distance of corresponding estimate from the median of all
%estimated points in the sample. this distance can be use as a measure of where the subject
%"would be" on the curve (if he were on the curve) hence providing a measure of severity. 
est_dot=[est_x',est_y'];

%computing distance from median point 
%computing the median of all estimated points, storing coordinates in ref
ref_x=median(est_x);
ref_y=median(est_y);
ref=[ref_x,ref_y];

%for all estimated points, compute distance from median
for i=1:length(est_dot)
est_d(i)=sqrt(sum((est_dot(i,:)-ref).^2));
%adjust sign
if est_y(i) > ref_y
    est_d(i) = -est_d(i);
end
end

%store useful variables in summary, ie x, y, distance from curve, estimated x, estimated
%y, estimated distance from median
summary=[PET/ScalingFactor,CSF,res',est_x'/ScalingFactor,est_y',est_d'];
%computes Z score for distance from curve
summary(:,7)=((summary(:,3))- mean(summary(:,3)))./std(summary(:,3));
%compute Z score for distance for distance from median
summary(:,8)=((summary(:,6))- mean(summary(:,6)))./std(summary(:,6));
%compute interaction between the two
summary(:,9)=summary(:,7).*summary(:,8);
%save output
save summary.mat summary; 
xlswrite('summary.xlsx',summary); 

%computing pseudo-R squared (using Euclidean distances instead of the standard y-distances)
%RSS PART%
squaredres= power(res,2);
RSS=sum(squaredres)

%TSS PART%
dot=[PET,CSF];
meanCSF=mean(CSF);
meanSUVr=mean(PET);
avg=[meanSUVr,meanCSF];

for i=1:length(dot)
d(i)=sqrt(sum((dot(i,:)-avg).^2));
end

squaredd=power(d,2);
TSS= sum(squaredd)

%compute R^2%
R2=1-(RSS/TSS)

%create various plots;
figure;histfit(res,30);fitdist(res','Normal')
figure;scatter(PET/ScalingFactor,CSF,25,summary(:,7),'filled');hold on;plot(x/ScalingFactor,y);colorbar;hcb=colorbar;title(hcb,'Aggregation');axis([-20 150 0 4000]);xlabel('PET');ylabel('CSF')
figure;scatter(CSF,PET/ScalingFactor,8,summary(:,7),'filled');hold on;plot(y,x/ScalingFactor);colorbar;hcb=colorbar;title(hcb,'Aggregation');axis([-1000 6000 -20 150]);ylabel('PET');xlabel('CSF')
figure;scatter(PET/ScalingFactor,CSF,25,summary(:,8),'filled');hold on;plot(x/ScalingFactor,y);colorbar;colormap(flipud(autumn));hcb=colorbar;title(hcb,'Severity');axis([-50 200 0 4000]);xlabel('PET');ylabel('CSF')
figure;scatter(CSF,PET/ScalingFactor,8,summary(:,8),'filled');hold on;plot(y,x/ScalingFactor);colorbar;hcb=colorbar;title(hcb,'Severity');axis([-1000 6000 -20 150]);ylabel('PET');xlabel('CSF')
%figure;scatter(PET/ScalingFactor,CSF,8,summary(:,9),'filled');hold on;plot(x/ScalingFactor,y);colorbar;hcb=colorbar;title(hcb,'Interaction');axis([-20 150 -1000 6000]);xlabel('PET');ylabel('CSF')
%figure;scatter(CSF,PET/ScalingFactor,8,summary(:,9),'filled');hold on;plot(y,x/ScalingFactor);colorbar;hcb=colorbar;title(hcb,'Interaction');axis([-1000 6000 -20 150]);ylabel('PET');xlabel('CSF')
%figure;histfit(est_d,30);fitdist(res','Normal')
%figure;scatter(PET/ScalingFactor,CSF,8,est_d,'filled');hold on;plot(x/ScalingFactor,y);colorbar;hcb=colorbar;title(hcb,'Severity');axis([-20 150 -1000 6000]);xlabel('PET');ylabel('CSF')


%computing residuals as euclidean distances from the curve, and
%corresponding estimates on the curve for each observations
function [d,est_x,est_y]=min_dist_point_curve(point,curve);

%computing distance of the point to any point on the curve
for i = 1:length(curve)
  d(i)=sqrt(sum((point-curve(i,:)).^2));
end

%getting the minimum distance
[d pos]=min(d);

%getting the coordiantes of the curve corresponding to the estimmte with minimum distance
est_x = curve(pos,1);
est_y = curve(pos,2);

%getting index corresponding to minimum distance
[mn idx] = min(abs(curve(:,1)-point(1)));

%getting sign of the distance based on the position of the point with min
%distance with respect to the curve
curvex = curve(idx,1);
curvey = curve(idx,2);

if point(2) < curvey
    d = -d;
    
end



