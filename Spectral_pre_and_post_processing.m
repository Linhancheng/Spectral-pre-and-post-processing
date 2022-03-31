%% Example_1
% This example illustrates the process of spectral preprocessing.
clc;
clear;
close all;

load('Example_1.mat')

% baseline correction
data_f = data(:,find(wn>900&wn<1800));
data_f_lab = irdata();
data_f_lab.X = data_f;
block = pre_bc_rubber();
block.flag_trim = 0;
data_f_rubber = block.use(data_f_lab);
data_f_rubber = data_f_rubber.X;

data_h = data(:,find(wn>2800&wn<3400));
data_h_lab = irdata();
data_h_lab.X = data_h;
data_h_rubber = block.use(data_h_lab);
data_h_rubber = data_h_rubber.X;

data_rubber = [data_h_rubber,data_f_rubber];

% normalization
data_rubber_vec = normaliz(data_rubber);

% second derivative conversion
data_rubber_vec_sg = savgol(data_rubber_vec,7,3,2);
data_rubber_vec_sg(:,[1:3,154:159,end-2:end]) = 0;

wn_cutted = [wn(find(wn>2800&wn<3400)) wn(find(wn>900&wn<1800))];

clear block data_f data_f_lab data_f_rubber data_h data_h_lab data_h_rubber
save('Example_1_result.mat','data_rubber_vec_sg','wn_cutted');
%% Example_2
% this example illustrates the use of algorithms to remove
% the non-cardiac spectra
clc;
clear;
close all;

load('Example_1_result.mat');

% K-means analysis
pixel = (1:1:size(data_rubber_vec_sg,1))';
cluster = kmeans(data_rubber_vec_sg,10);
pixel_myo = find(cluster(:,1)==mode(cluster));
data_rubber_vec_sg_myo=data_rubber_vec_sg(pixel_myo,:);

% PCA analysis
[~,score,~,~]=pca(data_rubber_vec_sg_myo);
score=roundn(score,-12);

figure
ConfidenceRegion(score(:,1:2),0.001,'norm');

f=findobj('Marker','+','Color','r');
X=get(f,'XData');
Y=get(f,'YData');
for i=1:size(X,2)
    pixel_out_byPCA(i,1)=find(score(:,1)==X(1,i));
end
pixel_myo(pixel_out_byPCA,:)=[];

pixel_outlier = setdiff(pixel, pixel_myo) ;

clear cluster f i pixel_out_byPCA X Y score
save('Example_2_result.mat','pixel_outlier','data_rubber_vec_sg',...
    'pixel','wn_cutted');
%% Example_3
% This example illustrates the use of code to average
% the 8*8 spectra corresponding to each nonoverlapping patch
% of 8*8 pixels (800*800 um2) to generate a single spectrum
clc;
clear;
close all;

load('Example_2_result.mat')

patch_pixel = 8;
data_rubber_vec_sg(pixel_outlier,:) = missing;

col = 47;
row = size(data_rubber_vec_sg,1)/col;

r_col = (col-rem(col,patch_pixel))/patch_pixel;
r_row = (row-rem(row,patch_pixel))/patch_pixel;
   
for j=1:row
    pixel_reshaped(j,1:col)=pixel(col*j-col+1:col*j,1);
end
pixel_reshaped(:,1:rem(col,patch_pixel)) = [];
pixel_reshaped(1:rem(row,patch_pixel),:) = [];

e = mat2cell(pixel_reshaped,ones(r_row,1)*patch_pixel,ones(r_col,1)*patch_pixel);
e = reshape(e,r_col*r_row,1);

for k = 1:r_col*r_row
    position = reshape(e{k},patch_pixel*patch_pixel,1);
    data_rubber_vec_sg_myo_mean(k,:) = nanmean(data_rubber_vec_sg(position,:));
end
TF = ismissing(data_rubber_vec_sg_myo_mean(:,1));
data_rubber_vec_sg_myo_mean(find(TF == 1),:) = [];

clear col e j k r_col r_row row TF pixel patch_pixel pixel_outlier pixel_reshaped position
save('Example_3_result.mat','data_rubber_vec_sg_myo_mean','wn_cutted');
%% Example_4
% this example illustrates the use of PCA toolbox to characterize
% the separation trend between the diabetic and control groups.
%
clc;
clear;
close all;

load('Example_4.mat');
model_pca = pca_model(data_rubb_vect_sg_myo,3,'cent');

figure;
axes1 = axes('Position',[0.13 0.11 0.8 0.82]);
hold(axes1,'on');
plot(model_pca.T(1:571,1),model_pca.T(1:571,3),'DisplayName','Diabetic','MarkerFaceColor',[1 0 0],...
    'MarkerEdgeColor',[0 0 0],...
    'MarkerSize',7,...
    'Marker','o',...
    'LineStyle','none');
plot(model_pca.T(572:end,1),model_pca.T(572:end,3),'DisplayName','Control',...
    'MarkerFaceColor',[0 0.45 0.74],...
    'MarkerEdgeColor',[0 0 0],...
    'MarkerSize',7,...
    'Marker','o',...
    'LineStyle','none');
ylabel('Scores on PC 3 - EV = 13.71%','FontWeight','bold');
xlabel('Scores on PC 1 - EV = 51.96%','FontWeight','bold');
box(axes1,'on');
set(axes1,'FontSize',20,'FontWeight','bold','LineWidth',1.5);
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.4 0.8 0.1 0.1],...
    'FontSize',20,...
    'EdgeColor',[1 1 1]);

%% Example_5
% this example illustrates the use of the Classification toolbox
% to differentiate different stages of diabetic mouse myocardia.
clc;
clear;
close all;

load('Example_5.mat');
res = plsdacompsel(data_rubb_vect_sg_myo_db_07_12,label_db_07_12,'auto','vene',10,'bayes');
cv = plsdacv(data_rubb_vect_sg_myo_db_07_12,label_db_07_12,4,'auto','vene',10,'bayes');
plsda_model_0712 = plsdafit(data_rubb_vect_sg_myo_db_07_12,label_db_07_12,4,'auto','vene','bayes');

figure;
axes1 = axes;
hold(axes1,'on');
plot(plsda_model_0712.T(1:203,1),plsda_model_0712.T(1:203,2),'DisplayName','Diabetic-7w',...
    'MarkerFaceColor',[0 0.48 0.5],...
    'MarkerEdgeColor',[0 0 0],...
    'MarkerSize',7,...
    'Marker','o',...
    'LineStyle','none');
plot(plsda_model_0712.T(204:end,1),plsda_model_0712.T(204:end,2),'DisplayName','Diabetic-12w',...
    'MarkerFaceColor',[1 0.35 0.37],...
    'MarkerEdgeColor',[0 0 0],...
    'MarkerSize',7,...
    'Marker','o',...
    'LineStyle','none');
ylabel('Scores on LV 2 - EV = 9.36%','FontWeight','bold');
xlabel('Scores on LV 1 - EV = 10.64%','FontWeight','bold');
box(axes1,'on');
set(axes1,'DataAspectRatio',[1 1 1],'FontSize',20,'FontWeight','bold',...
    'LineWidth',1.5,'PlotBoxAspectRatio',[30 30 1]);
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.57 0.81 0.14 0.1],...
    'FontSize',20,...
    'EdgeColor',[1 1 1]);

res = plsdacompsel(data_rubb_vect_sg_myo_db_12_21,label_db_12_21,'auto','vene',10,'bayes');
cv = plsdacv(data_rubb_vect_sg_myo_db_12_21,label_db_12_21,4,'auto','vene',10,'bayes');
plsda_model_1221 = plsdafit(data_rubb_vect_sg_myo_db_12_21,label_db_12_21,4,'auto','vene','bayes');

figure;
axes1 = axes;
hold(axes1,'on');
plot(plsda_model_1221.T(1:174,1),plsda_model_1221.T(1:174,2),'DisplayName','Diabetic-12w',...
    'MarkerFaceColor',[1 0.35 0.37],...
    'MarkerEdgeColor',[0 0 0],...
    'MarkerSize',7,...
    'Marker','o',...
    'LineStyle','none');
plot(plsda_model_1221.T(175:end,1),plsda_model_1221.T(175:end,2),'DisplayName','Diabetic-21w',...
    'MarkerFaceColor',[1 0.7 0],...
    'MarkerEdgeColor',[0 0 0],...
    'MarkerSize',7,...
    'Marker','o',...
    'LineStyle','none');
ylabel('Scores on LV 2 - EV = 12.69%','FontWeight','bold');
xlabel('Scores on LV 1 - EV = 17.86%','FontWeight','bold');
box(axes1,'on');
set(axes1,'DataAspectRatio',[1 1 1],'FontSize',20,'FontWeight','bold',...
    'LineWidth',1.5,'PlotBoxAspectRatio',[30 30 1]);
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.57 0.81 0.14 0.1],...
    'FontSize',20,...
    'EdgeColor',[1 1 1]);
%% Example_6
% this example illustrates the use of the RandomForest package
% in combination with PLS-DA to interpretate the biochemical changes in
% diabetic mice's myocardia.
clc;
clear;
close all;

load('Example_6.mat');

% rf_model = classRF_train(data_rubb_vect_sg_myo_db,label_db);

for i = 1:571
[vote(i,1),posi(i,1)]=max(rf_model.votes(i,1:3));
end
a = cell(571,1)
for i = 1:571
    if posi(i) == 1 a{i}='Diabetic-7w'
    elseif posi(i) == 2 a{i}='Diabetic-12w'
    elseif posi(i) == 3 a{i}='Diabetic-21w'
    end
end
predictedLabels = categorical(a)
clear a 

b = cell(571,1)
parfor i = 1:571
    if label_db(i) == 1 b{i}='Diabetic-7w'
    elseif label_db(i) == 2 b{i}='Diabetic-12w'
    elseif label_db(i) == 3 b{i}='Diabetic-21w'    
    end
end
trueLabels = categorical(b)
clear b i 

figure
cm = confusionchart(trueLabels,predictedLabels, ...
    'ColumnSummary','column-normalized', ...
    'RowSummary','row-normalized');

gini = rf_model.importance
gini_sorted = sort(gini,'descend')
for i = 1:19
   var_position(i) = find(gini==gini_sorted(i))
end
data_selected = data_rubb_vect_sg_myo_db(:,var_position)

res = plsdacompsel(data_selected,label_db,'auto','vene',10,'bayes');
plsda_model = plsdafit(data_selected,label_db,4,'auto','vene','bayes');

figure;
axes1 = axes('Position',[0.13 0.11 0.8 0.815]);
hold(axes1,'on');
plot(plsda_model.T(1:203,1),plsda_model.T(1:203,2),'DisplayName','Diabetic-7w',...
    'MarkerFaceColor',[0 0.48 0.5],...
    'MarkerEdgeColor',[0 0 0],...
    'MarkerSize',7,...
    'Marker','o',...
    'LineStyle','none');
plot(plsda_model.T(204:377,1),plsda_model.T(204:377,2),'DisplayName','Diabetic-12w',...
    'MarkerFaceColor',[1 0.35 0.37],...
    'MarkerEdgeColor',[0 0 0],...
    'MarkerSize',7,...
    'Marker','o',...
    'LineStyle','none');
plot(plsda_model.T(378:end,1),plsda_model.T(378:end,2),'DisplayName','Diabetic-21w',...
    'MarkerFaceColor',[1 0.7 0],...
    'MarkerEdgeColor',[0 0 0],...
    'MarkerSize',7,...
    'Marker','o',...
    'LineStyle','none');
ylabel('Scores on LV 2 - EV = 17.36%','FontWeight','bold');
xlabel('Scores on LV 1 - EV = 22.48%','FontWeight','bold');
box(axes1,'on');
set(axes1,'FontSize',20,'FontWeight','bold','LineWidth',1.5);
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.12 0.77 0.14 0.15],...
    'FontSize',20,...
    'EdgeColor',[1 1 1]);
clear cm data_selected gini gini_sorted i posi res trueLabels var_position vote
clear predictedLabels


%% Example_7
% this example illustrates the use of the RandomForest package
% to develop the classification model for predicting DbCM.
clc;
clear;
close all;

load('Example_7.mat');

% rf_model = classRF_train(data_train,label_train) 

for i = 1:size(label_train,1)
[vote(i,1),posi(i,1)]=max(rf_model.votes(i,1:2));
end

a=cell (size(label_train,1),1)
parfor i = 1:size(label_train,1)
    if posi(i) == 1 a{i}='Diabetic'
    elseif posi(i) == 2 a{i}='Control'
    
    end
end
predictedLabels = categorical(a)
clear a 

b = cell(size(label_train,1),1)
parfor i = 1:size(label_train,1)
    if label_train(i) == 1 b{i}='Diabetic'
    elseif label_train(i) == 2 b{i}='Control'
    
    end
end
trueLabels = categorical(b)
clear b i 
figure
cm = confusionchart(trueLabels,predictedLabels, ...
    'ColumnSummary','column-normalized', ...
    'RowSummary','row-normalized');

% external validation
[y votes] = classRF_predict(data_test,rf_model)

a=cell (size(label_test,1),1)
parfor i = 1:size(label_test,1)
    if y(i) == 1 a{i}='Diabetic'
    elseif y(i) == 2 a{i}='Control'
    
    end
end
predictedLabels = categorical(a)
clear a 

b = cell(size(label_test,1),1)
parfor i = 1:size(label_test,1)
    if label_test(i) == 1 b{i}='Diabetic'
    elseif label_test(i) == 2 b{i}='Control'
    
    end
end
trueLabels = categorical(b)
clear b i 
figure
cm = confusionchart(trueLabels,predictedLabels, ...
    'ColumnSummary','column-normalized', ...
    'RowSummary','row-normalized');
clear cm posi vote votes y trueLabels predictedLabels
