celltype = 'EndothelialCells'; % 'Fibroblasts' / 'Macrophages' / 'EndothelialCells'
wt_threshold = 0.5; %0.3; 0.4; % Above are specialists
gen_threshold = 0; 0.3; % Below are generalists -> remove generalists from figure.
switch celltype
    case 'EndothelialCells'
        arcOrigFile = 'C:\Users\amorp\Documents\Work\Adler\Projects\Supportive_Archetypes\ParTi_Endo_Manual_Muscle_6arcs\Endo_6_arcOrig.mat';
        geneExpression = table2array(readtable('C:\Users\amorp\Documents\Work\Adler\Projects\Supportive_Archetypes\ParTi_Endo_Manual_Muscle_6arcs\geneExpression.csv'));
        meta = readtable('C:\Users\amorp\Documents\Work\Adler\Projects\Supportive_Archetypes\Endo_Meta.csv','ReadRowNames',true);

    case 'Fibroblasts'
        arcOrigFile = 'C:\Users\amorp\Large_Files\archetypes_geneExpressionMats\Fibro_arcOrig.mat';
        geneExpression = table2array(readtable('C:\Users\amorp\Large_Files\archetypes_geneExpressionMats\fibro.csv'));
        meta = readtable('C:\Users\amorp\Large_Files\archetypes_geneExpressionMats\fibro_metadata.csv');

    case 'Macrophages'
        arcOrigFile = 'C:\Users\amorp\Large_Files\archetypes_geneExpressionMats\Mac_rev_muscle_stomech_k5_arcOrig.mat';
        geneExpression_table = readtable('C:\Users\amorp\Large_Files\archetypes_geneExpressionMats\mac_rev_muscle_normalized.csv');
        geneExpression = table2array(geneExpression_table(:,~strcmp(geneExpression_table.Properties.VariableNames,'cell')));
        clear geneExpression_table
        meta = readtable('C:\Users\amorp\Large_Files\archetypes_geneExpressionMats\Mac_rev_muscle_stomech_meta.csv');
        
end
% geneExpression_2 = geneExpression(1:2:end,:);
geneExpression_2 = geneExpression;
% tissues = meta.tissue_in_publication(1:2:end);
tissues = meta.tissue_in_publication;

savefolder = 'C:\Users\amorp\Documents\Work\Adler\Projects\Supportive_Archetypes\tsne_with_arcs';
load(arcOrigFile,'arcOrig');
tissue_legend = {'Bladder','Fat','Heart','Kidney','Large intestine','Liver',...
    'Lung','Pancreas','Skin','Small intestine','Spleen','Thymus','Tongue',...
    'Trachea','Stomach','Muscle'};
tissue_legend_clr = {'#3953A4','#ED2024','#6ABD45','#0F1031' ,'#E03F97' ,...
     '#1F5429','#FED304','#4996D2' ,'#9A4D42','#DF7D26' ,'#714EA0' ,...
    '#1E9698','#E3AFD1','#A93493','#C49A6C','#5C8A00'};
ARCHETYPE_PALETTE = hex2rgb({'#FFFF33','#377EB8','#4DAF4A',...
    '#E41A1C','#984EA3','#FF7F00'});
Generalist_clr = hex2rgb('#808080');
no_arc_clr = [0.8 0.8 0.8];
mrksz = 2;
%% run TSNE
if ~exist([savefolder filesep 'tsne_' celltype '.mat'],'file')
    Y_tsne = tsne(geneExpression_2,'Algorithm','exact','Distance','cosine');
    save([savefolder filesep 'tsne_' celltype '.mat'],'Y_tsne');
else
    load([savefolder filesep 'tsne_' celltype '.mat'],'Y_tsne');
end
%% Calc weights

% Preallocate result
num_points = size(geneExpression_2, 1);
archetype_weights = zeros(num_points, size(arcOrig, 1));

% Loop over each point
for i = 1:num_points
    point = geneExpression_2(i, :);
    weights = calculate_archetype_weights(point, arcOrig);
    archetype_weights(i, :) = weights;
end
%% Plot

figure(1); clf; 
% set(gcf,'Units','normalized','Position',[0.05 0.2 0.85 0.5],'name',celltype)
set(gcf,'Units','centimeters','Position',[2 2 18 7],'name',celltype)

ax_tsne = axes('Units','normalized','Position',[0.02 0.07 0.38 0.9]);

hold on
% arc_clrs = [1 1 0 ; 0 0 1 ; 0 1 0; 1 0 0; 0.8 0 0.8; 1 0.5 0];
arc_clrs = ARCHETYPE_PALETTE;
cell_color = zeros(size(archetype_weights,1),3);
gen_count = 0;
cell_order = zeros(size(archetype_weights,1),1);
for cl = 1:size(archetype_weights,1)
    [max_arc,max_arc_ind] = max(archetype_weights(cl,:));
    if max_arc>=wt_threshold
        cell_color(cl,:) = arc_clrs(max_arc_ind,:)*max_arc;
        cell_order(cl) = 1;
    elseif max_arc<=gen_threshold
        cell_color(cl,:) = Generalist_clr;
        gen_count = gen_count+1;
    else
        cell_color(cl,:) = no_arc_clr;
    end
    
    % scatter(ax_tsne,Y_tsne(cl,1), Y_tsne(cl,2),mrksz,cell_color(cl,:),'filled')

end

[~,plot_order] = sort(cell_order);
scatter(ax_tsne,Y_tsne(plot_order,1), Y_tsne(plot_order,2),mrksz,cell_color(plot_order,:),'filled')

fprintf('%d generalists\n',gen_count)
set(gca,'XTick',[],'YTick',[])
% Arcs Colorbar:
for arcs = 1:size(archetype_weights,2)
    ax_lgd = axes('Units','normalized','Position',[0.43 0.075+(arcs-1)*0.117 0.01 0.09]);

    imagesc(ax_lgd,(wt_threshold:0.05:1)'); 
    set(ax_lgd,'Colormap',(wt_threshold:0.05:1)'.*arc_clrs(arcs,:),'YDir','normal','YTick',[],'XTick',[]);
    ylabel(['Arc #' num2str(arcs)])
    set(gca,'FontSize',6)
end
% Add generalists and unassigned:
ax_lgd_unass = axes('Units','normalized','Position',[0.43 0.075+(arcs)*0.117 0.01 0.09]);
imagesc(ax_lgd_unass,1);
set(ax_lgd_unass,'Colormap',no_arc_clr,'YDir','normal','YTick',[],'XTick',[]);
ylabel('Unasgd.')
set(gca,'FontSize',6)

% ax_lgd_gen = axes('Units','normalized','Position',[0.43 0.075+(arcs)*0.117 0.01 0.09]);
% imagesc(ax_lgd_gen,1);
% set(ax_lgd_gen,'Colormap',Generalist_clr,'YDir','normal','YTick',[],'XTick',[]);
% ylabel('Gen.')
% set(gca,'FontSize',6)
% ax_lgd_unass = axes('Units','normalized','Position',[0.43 0.075+(arcs+1)*0.117 0.01 0.09]);
% imagesc(ax_lgd_unass,1);
% set(ax_lgd_unass,'Colormap',no_arc_clr,'YDir','normal','YTick',[],'XTick',[]);
% ylabel('Unasgd.')
% set(gca,'FontSize',6)

ax_tsne_tis = axes('Units','normalized','Position',[0.47 0.072 0.38 0.9]);
hold on
[a,b] = ismember(replace(lower(tissues),' ','_'),replace(lower(tissue_legend)',' ','_'));
cell_color = tissue_legend_clr(b);
scatter(ax_tsne_tis,Y_tsne(:,1), Y_tsne(:,2),mrksz,hex2rgb(cell_color),'filled')
set(gca,'XTick',[],'YTick',[])

ax_lgd_tis = axes('Units','normalized','Position',[0.86 0.072 0.014 0.9]);
hold on
space = 1.12;
for ts = 1:numel(tissue_legend)
    plot(ax_lgd_tis,0,0.1+(numel(tissue_legend)-ts+1)*space,'o','MarkerFaceColor',hex2rgb([tissue_legend_clr{ts}]),...
        'MarkerEdgeColor',hex2rgb([tissue_legend_clr{ts}]),'MarkerSize',4)
    text(ax_lgd_tis,2,0.1+(numel(tissue_legend)-ts+1)*space,tissue_legend{ts},'FontSize',8); 
end
set(gca,'Visible','off')
savefig([ savefolder filesep 'tsne_arc_tissues_' celltype])
exportgraphics(figure(1),[ savefolder filesep 'tsne_arc_tissues_' celltype '.pdf'],'ContentType','vector');
%%
function weights_norm = calculate_archetype_weights(point, arcOrig)
    % Check dimensionality
    if size(arcOrig,2) ~= length(point)
        error('Mismatch in dimensions: archetypes have %d features, but point has %d features.', ...
            size(arcOrig,2), length(point));
    end

    % Extend archetypes and point for homogeneous coordinates
    archetypes_hmg = [arcOrig, ones(size(arcOrig,1), 1)];
    point_hmg = [point, 1];

    % Compute A and b
    A = archetypes_hmg * archetypes_hmg';
    b = archetypes_hmg * point_hmg';

    % Solve for weights
    weights = A \ b;

    % Normalize
    weights_norm = weights / sum(weights);
end
