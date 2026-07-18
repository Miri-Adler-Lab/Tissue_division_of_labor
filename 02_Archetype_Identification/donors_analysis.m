%% load data:
celltype = 'Fibro'; % Endo / Macro / Fibro
dim = 8;
algNum = 5;
binSize = 0.05;

work_folder = 'C:\Users\amorp\Documents\Work\Adler\Projects\Supportive_Archetypes\donor_enrichment\';
celltype_folder = [work_folder filesep celltype];
if ~exist(celltype_folder,'dir')
    mkdir(celltype_folder)
end
donor_table_file = [work_folder filesep 'donor_info_added_TSP9_TSP20.xlsx'];
donor_table = readtable(donor_table_file,'Sheet','Manually_Curated_metadata');
fprintf('Loading %s data...\n',celltype)
switch celltype
    case 'Macro'
        meta_tbl = readtable('C:\Users\amorp\Large_Files\archetypes_geneExpressionMats\Mac_rev_muscle_stomech_meta.csv','ReadRowNames',true);
        geneExpression_table = readtable('C:\Users\amorp\Large_Files\archetypes_geneExpressionMats\mac_rev_muscle_normalized.csv');
        geneNamesStrcut = load('C:\Users\amorp\Documents\Work\Adler\Projects\Supportive_Archetypes\Figure_human_mouse_cmp\Macrophages\human_geneNames.mat','geneNames');
        geneNames = geneNamesStrcut.geneNames;
        geneExpression = table2array(geneExpression_table(:,~strcmp(geneExpression_table.Properties.VariableNames,'cell')));
        load('C:\Users\amorp\Large_Files\archetypes_geneExpressionMats\Mac_rev_muscle_stomech_k5_arc.mat','arc')
        load('C:\Users\amorp\Large_Files\archetypes_geneExpressionMats\Mac_rev_muscle_stomech_k5_arcOrig.mat','arcOrig')
    case 'Fibro'
        meta_tbl = readtable('C:\Users\amorp\Large_Files\archetypes_geneExpressionMats\fibro_metadata.csv','ReadRowNames',true);
        geneExpression = table2array(readtable('C:\Users\amorp\Large_Files\archetypes_geneExpressionMats\fibro.csv'));
        geneNamesStrcut = load('C:\Users\amorp\Documents\Work\Adler\Projects\Supportive_Archetypes\Figure_human_mouse_cmp\Fibroblasts\human_geneNames.mat','geneNames');
        geneNames = geneNamesStrcut.geneNames;

        load('C:\Users\amorp\Large_Files\archetypes_geneExpressionMats\Fibro_arc.mat','arc')
        load('C:\Users\amorp\Large_Files\archetypes_geneExpressionMats\Fibro_arcOrig.mat','arcOrig')
    case 'Endo'
        meta_tbl = readtable('C:\Users\amorp\Documents\Work\Adler\Projects\Supportive_Archetypes\Endo_Meta.csv','ReadRowNames',true);
        geneExpression = table2array(readtable('C:\Users\amorp\Documents\Work\Adler\Projects\Supportive_Archetypes\ParTi_Endo_Manual_Muscle_6arcs\geneExpression.csv'));
        geneNamesStrcut = load('C:\Users\amorp\Documents\Work\Adler\Projects\Supportive_Archetypes\Figure_human_mouse_cmp\EndothelialCells\human_geneNames.mat','geneNames');
        geneNames = geneNamesStrcut.geneNames;

        load('C:\Users\amorp\Documents\Work\Adler\Projects\Supportive_Archetypes\ParTi_Endo_Manual_Muscle_6arcs\Endo_6_arc.mat','arc')
        load('C:\Users\amorp\Documents\Work\Adler\Projects\Supportive_Archetypes\ParTi_Endo_Manual_Muscle_6arcs\Endo_6_arcOrig.mat','arcOrig')

end
deathCauseShort = {'Stroke','Anoxi','Head trauma','Pneumonia',...
    'Intracranial hemorrhage','Respiratory failure','NA'};
ethnicityNames = { 'Asian','Black','Hispanic','White'};         
donor_clrs = hex2rgb({'#3953A4','#ED2024','#6ABD45','#B3B3B3' ,'#E03F97' ,...
     '#1F5429','#FED304','#4996D2' ,'#9A4D42','#DF7D26' ,'#714EA0' ,...
    '#1E9698','#E3AFD1','#A93493','#75C3CD'});
tissue_legend = {'Bladder','Fat','Heart','Kidney','Large intestine','Liver',...
    'Lung','Pancreas','Skin','Small intestine','Spleen','Thymus','Tongue',...
    'Trachea','Stomach','Muscle'};
tissue_legend_clr = hex2rgb({'#3953A4','#ED2024','#6ABD45','#0F1031' ,'#E03F97' ,...
     '#1F5429','#FED304','#4996D2' ,'#9A4D42','#DF7D26' ,'#714EA0' ,...
    '#1E9698','#E3AFD1','#A93493','#C49A6C','#5C8A00'});
ARCHETYPE_PALETTE = hex2rgb({'#FFFF33','#377EB8','#4DAF4A',...
    '#E41A1C','#984EA3','#FF7F00'});
%%
fprintf('Updating metadata per cell...\n')
switch celltype
    case {'Endo','Macro'}
        donor_per_cell = meta_tbl.donor_id;
        udonors = unique(donor_per_cell,'stable');
    case 'Fibro'
        donor_per_cell = cellfun(@(y) y{1},cellfun(@(x) strsplit(x,'_'),meta_tbl.donor_id,'UniformOutput',false),'UniformOutput',false);
        udonors = unique(donor_per_cell,'stable');
end

tissues_per_cel = meta_tbl.tissue_in_publication;
utissues = unique(meta_tbl.tissue_in_publication,'stable');
don_tis_mat = zeros(numel(udonors),numel(utissues));
deathCause =  zeros(numel(udonors),numel(deathCauseShort));
deathCauseMat = zeros(size(meta_tbl,1),numel(deathCauseShort));
ethnicity = zeros(numel(udonors),numel(ethnicityNames));
ethnicityMat = zeros(size(meta_tbl,1),numel(ethnicityNames));
age = [];
sex = [];
bmi = [];
alcohol = [];
drugs = [];
diabetes = [];
cardio = [];

for d = numel(udonors):-1:1
    for t = 1:numel(utissues)
        don_tis_mat(d,t) = sum(strcmp(donor_per_cell,udonors{d}) & strcmp(meta_tbl.tissue_in_publication,utissues{t}));
    end
    dn_num = str2double(erase(udonors{d},'TSP'));
    dn_ln = find(str2double(erase(donor_table.Donor,'Donor '))==dn_num);
    if ~isempty(dn_ln)
        age(d) = donor_table.Age(dn_ln);
        sex(d) = strcmp(donor_table.Sex(dn_ln),'M'); % 1-M, 0-F
        bmi(d) = donor_table.BMI_kg_m2_(dn_ln);
        % alcohol(d) = strcmp(donor_table.Alcohol(dn_ln),'Y');
        % drugs(d) = strcmp(donor_table.IVDrugAbuse(dn_ln),'Y');
        % diabetes(d) = strcmp(donor_table.Diabetes(dn_ln),'Y');
        alcohol(d) = find(strcmp(donor_table.Alcohol(dn_ln),{'NA','N','Y'}))-2;  % -1:NA, 0:N, 1:Y
        drugs(d) = find(strcmp(donor_table.IVDrugAbuse(dn_ln),{'NA','N','Y'}))-2;
        diabetes(d) = find(strcmp(donor_table.Diabetes(dn_ln),{'NA','N','Y'}))-2;
        cardio(d) = strcmp(donor_table.Hypertension(dn_ln),'Y') | strcmp(donor_table.CoronaryArteryDisease(dn_ln),'Y');
        ethnicity(d,:) = cellfun(@(x) contains(donor_table.Ethnicity{dn_ln},x,'IgnoreCase',true),ethnicityNames);
        deathCause(d,:) = cellfun(@(x) contains(donor_table.CauseOfDeath{dn_ln},x,'IgnoreCase',true),deathCauseShort);
    end
    lines_in_meta = find(strcmp(donor_per_cell,udonors{d}));
    meta_tbl.donor_age(lines_in_meta) = age(d);
    meta_tbl.donor_bmi(lines_in_meta) = bmi(d);
    meta_tbl.donor_sex(lines_in_meta) = sex(d);
    meta_tbl.donor_alcohol_usage(lines_in_meta) = alcohol(d);
    meta_tbl.donor_drugs_usage(lines_in_meta) = drugs(d);
    meta_tbl.donor_diabetes(lines_in_meta) = diabetes(d);
    meta_tbl.donor_cardio(lines_in_meta) = cardio(d);
    for cause = 1:numel(deathCauseShort)
        meta_tbl.(['donor_deathCause_' deathCauseShort{cause}])(lines_in_meta) = deathCause(d,cause);
    end
    for etn = 1:numel(ethnicityNames)
        meta_tbl.(['donor_Ethnicity_' ethnicityNames{etn}])(lines_in_meta) = ethnicity(d,etn);
    end
end
for etn = 1:numel(ethnicityNames)
    ethnicityMat(:,etn) = meta_tbl.(['donor_Ethnicity_' ethnicityNames{etn}]);
end
for cause = 1:numel(deathCauseShort)
    deathCauseMat(:,cause) = meta_tbl.(['donor_deathCause_' deathCauseShort{cause}]);
end

%% Calculate all data enrichment
ParTi_MetaEnrich_FileName_csv = [celltype_folder filesep celltype '_Enrich_All'];
meta_cont = [meta_tbl.donor_age meta_tbl.donor_bmi];
titles_cont = {'Age','BMI'};
meta_discrete = [meta_tbl.donor_sex ethnicityMat meta_tbl.donor_alcohol_usage meta_tbl.donor_drugs_usage meta_tbl.donor_diabetes ...
        meta_tbl.donor_cardio deathCauseMat];
titles_discrete = ['Sex',ethnicityNames,'Alcohol','Drugs','Diabities','Cardio',deathCauseShort];
if ~exist([ParTi_MetaEnrich_FileName_csv '_continuous_All.csv'],'file')   
    ParTI_lite(geneExpression,algNum,dim,titles_discrete,...
            meta_discrete,0,titles_cont,meta_cont,[],binSize,ParTi_MetaEnrich_FileName_csv,arcOrig);
else
    fprintf('ParTI data exists.\n')
end
close all
%%
[~,i_donor] = sort(str2double(erase(udonors,'TSP')));
wd = 0.017;
mrg = 0.015;
bot = 0.16;
left = 0.04;
hgt = 0.75;
mat_w = 0.4;
titangle = -25;
titfont = 12;

figure(1); clf; 
set(gcf,'units','normalized','position',[0 0.045 1 0.95])
ax_mat = axes('units','normalized','position',[left bot mat_w hgt]);
imagesc(don_tis_mat(i_donor,:)./sum(don_tis_mat(i_donor,:)))
xticks(1:numel(utissues))
xticklabels(replace(utissues,'_',' '))
yticks(1:numel(udonors))
yticklabels(erase(udonors(i_donor),'TSP'))
ylabel('Donor')
xlabel('Tissue')
title('Cell Origin')
set(ax_mat,'Colormap',brewermap(100,'RdPu'),'FontSize',13)

ax_age = axes('units','normalized','position',[left+mat_w+mrg bot wd hgt]);
imagesc(ax_age,age(i_donor)')
text(ones(numel(i_donor),1),1.2:numel(age(i_donor))+0.2,...
    num2str(age(i_donor)'),'Rotation',90,'FontSize',12)
title('Age','Rotation',titangle,'FontSize',titfont)
set(ax_age,'Colormap',brewermap(20,'YlOrBr'),'YTick',[],'XTick',[])

ax_bmi = axes('units','normalized','position',[left+mat_w+mrg+1*(mrg+wd) bot wd hgt]);
imagesc(ax_bmi,bmi(i_donor)')
text(ones(numel(i_donor),1),1.2:numel(bmi(i_donor))+0.2,...
    num2str(round(bmi(i_donor))'),'Rotation',90,'FontSize',12)
title('BMI','Rotation',titangle,'FontSize',titfont)
set(ax_bmi,'Colormap',brewermap(20,'YlOrBr'),'YTick',[],'XTick',[])

ax_sex = axes('units','normalized','position',[left+mat_w+mrg+2*(mrg+wd) bot wd hgt]);
imagesc(sex(i_donor)')
title('Sex','Rotation',titangle,'FontSize',titfont)
set(ax_sex,'Colormap',brewermap(2,'PiYG'),'YTick',[],'XTick',[])

y_n_na_clrmp = [1 1 1; 0.656    0.777    0.656 ;0.9373    0.5412    0.3843];
y_n_clrmp = [0.656    0.777    0.656 ;0.9373    0.5412    0.3843];

ax_alch = axes('units','normalized','position',[left+mat_w+mrg+3*(mrg+wd) bot wd hgt]);
imagesc(alcohol(i_donor)')
title('Alcohol','Rotation',titangle,'FontSize',titfont)
% set(ax5,'Colormap',flipud(brewermap(3,'RdGy')),'YTick',[],'XTick',[])
if ismember(-1,alcohol(i_donor)')
    clrmp = y_n_na_clrmp;
else
    clrmp = y_n_clrmp;
end
set(ax_alch,'Colormap',clrmp,'YTick',[],'XTick',[])

ax_drug = axes('units','normalized','position',[left+mat_w+mrg+4*(mrg+wd) bot wd hgt]);
imagesc(drugs(i_donor)')
title('Drugs','Rotation',titangle,'FontSize',titfont)
if ismember(-1,drugs(i_donor)')
    clrmp = y_n_na_clrmp;
else
    clrmp = y_n_clrmp;
end
set(ax_drug,'Colormap',clrmp,'YTick',[],'XTick',[])

ax_diab = axes('units','normalized','position',[left+mat_w+mrg+5*(mrg+wd) bot wd hgt]);
imagesc(diabetes(i_donor)')
title('Diabetes','Rotation',titangle,'FontSize',titfont)
if ismember(-1,diabetes(i_donor)')
    clrmp = y_n_na_clrmp;
else
    clrmp = y_n_clrmp;
end
set(ax_diab,'Colormap',clrmp,'YTick',[],'XTick',[])

ax_card = axes('units','normalized','position',[left+mat_w+mrg+6*(mrg+wd) bot wd hgt]);
imagesc(cardio(i_donor)')
title('Cardio','Rotation',titangle,'FontSize',titfont)
if ismember(-1,cardio(i_donor)')
    clrmp = y_n_na_clrmp;
else
    clrmp = y_n_clrmp;
end
set(ax_card,'Colormap',clrmp,'YTick',[],'XTick',[])

ax_ethn = axes('units','normalized','position',[left+mat_w+mrg+7.1*(mrg+wd) bot wd hgt]);
[~,ethnicity_vec] = max(ethnicity(i_donor,:),[],2);
imagesc(ethnicity_vec)
title('Ethnicity','Rotation',titangle,'FontSize',titfont)
clrmap_ethnicity = flipud(brewermap(numel(ethnicityNames),'Set2'));
for et = 1:numel(ethnicityNames)
    [~,et_clr] = ismember(ethnicityNames{et},ethnicityNames);
    text(1.7,0.6+0.6*et,ethnicityNames{et},'Color',clrmap_ethnicity(et_clr,:),...
        'FontSize',13,'FontWeight','bold')
end
set(ax_ethn,'Colormap',clrmap_ethnicity,'YTick',[],'XTick',[])

ax_death = axes('units','normalized','position',[left+mat_w+2*mrg+9.2*(mrg+wd) bot wd hgt]);
[~,deathCause_vec] = max(deathCause(i_donor,:),[],2);
imagesc(deathCause_vec)
title({'Cause of' ;'Death'},'FontSize',titfont,'Rotation',titangle)
uniq_deathCause_vec = unique(deathCause_vec,'stable');
% clrmap_deathcause = flipud(brewermap(numel(uniq_deathCause_vec),'Set2'));
clrmap_deathcause = brewermap(numel(deathCauseShort),'Set2');
set(ax_death,'Colormap',clrmap_deathcause(sort(uniq_deathCause_vec),:),'YTick',[],'XTick',[])

for dc = 1:numel(unique(deathCause_vec,'stable'))
    [~,d_color] = ismember(deathCauseShort{uniq_deathCause_vec(dc)},deathCauseShort);
    text(1.7,0.6+0.6*dc,deathCauseShort{uniq_deathCause_vec(dc)},...
        'Color',clrmap_deathcause(d_color,:),'FontSize',13,'FontWeight','bold')
end

exportgraphics(figure(1),[celltype_folder filesep celltype '_cell_origin_meta_distribution.pdf'],'ContentType','vector');
close all

%% check donor bias on PCA plot
if ~exist([celltype_folder filesep 'pca.mat'],'file')
    [~,score] = pca(geneExpression);
    save([celltype_folder filesep 'pca.mat'],'score')
else
    load([celltype_folder filesep 'pca.mat'],'score')
end
%% Plot each feature with archetypes
if ~exist([celltype_folder filesep 'feature_plots'],'dir')
    mkdir([celltype_folder filesep 'feature_plots'])
end

if ~exist([celltype_folder filesep 'feature_plots'],'dir')
    mkdir([celltype_folder filesep 'feature_plots'])
end
g = str2double(erase(donor_per_cell,'TSP'));
[~,~,gi] = unique(g,'stable');
c = donor_clrs(gi,:);
[a,b] = ismember(replace(lower(tissues_per_cel),'_',' '),lower(tissue_legend));
tc = tissue_legend_clr(b,:);
sig_cont_tbl = readtable([ParTi_MetaEnrich_FileName_csv '_continuous_significant.csv']);
sig_disc_tbl = readtable([ParTi_MetaEnrich_FileName_csv '_discrete_significant.csv']);

%% write table and plot true enrichment cases
close all
% Set bins: 
numOfBins = round(1 / binSize);
[DataPointsInd, distances] = sortDataByDistance(score(:,1:size(arc,2)),arc);
[Numarchs, numDataPoints] = size(DataPointsInd);
breakPoints = floor(linspace(0.5, numDataPoints + 0.5, numOfBins+1));
numPointInBin = diff(breakPoints);
breakPoints = breakPoints(2:end);
numFeatures = 1;
enrich_tbl_cont = table();
spec_th = 0.7;
PoQ_th_disc = 1.5;
PoQ_th_cont = 1.2;
if ~exist([celltype_folder filesep 'parti'],'dir')
    mkdir([celltype_folder filesep 'parti'])
end
for arcft = 1:size(sig_cont_tbl,1)
    enr_feat = sig_cont_tbl.FeatureName{arcft};
    enrArcNum = sig_cont_tbl.archetype_(arcft);
    enrich_tbl_cont.Arc(arcft) = enrArcNum;
    enrich_tbl_cont.Feature{arcft} = enr_feat;
    EnMatCont = meta_cont(:,strcmp(enr_feat,titles_cont)); % metadata values
    tempEnrich =  EnMatCont(DataPointsInd(enrArcNum,:),:);    
    Binned = mat2cell(tempEnrich,numPointInBin, numFeatures);
    DataPointsIndBinned = mat2cell(DataPointsInd(enrArcNum,:)',numPointInBin, numFeatures);
    BinnedDist = mat2cell(distances(enrArcNum,:)',numPointInBin, numFeatures);
    % [r,p] = corrcoef([BinnedDist{1}; BinnedDist{2}],[Binned{1} ; Binned{2}]);
    % Calculate Q: The global baseline median BMI of the entire dataset
    Q = nanmedian(EnMatCont);
    [Aid,Bid] = groupcounts(donor_per_cell(DataPointsIndBinned{1}(EnMatCont(DataPointsIndBinned{1})>Q)));
    [Ats,Bts] = groupcounts(meta_tbl.tissue_in_publication(DataPointsIndBinned{1}(EnMatCont(DataPointsIndBinned{1})>Q)));

    enrich_tbl_cont.MaxTissueAbundance(arcft) = max(Ats)./sum(Ats);
    enrich_tbl_cont.MaxDonorAbundance(arcft) = max(Aid)./sum(Aid);
    enrich_tbl_cont.DonorsWithFeature(arcft) = numel(unique(donor_per_cell(EnMatCont>Q)));

    PoverQ_continuous = zeros(Numarchs,numOfBins);
    for arch = 1:Numarchs
        % Divide metadata to bins
        tempEnrich =  EnMatCont(DataPointsInd(arch,:),:);    
        Binned = mat2cell(tempEnrich,numPointInBin, numFeatures);
        DataPointsIndBinned = mat2cell(DataPointsInd(arch,:)',numPointInBin, numFeatures);
        BinnedDist = mat2cell(distances(arch,:)',numPointInBin, numFeatures);
        % Calculate P: The median BMI for each bin
        P = cellfun(@(x) nanmedian(x), Binned, 'UniformOutput', true);     
        % Calculate PoverQ: Fold change of feature relative to the baseline
        PoverQ_continuous(arch,:) = P ./ Q;
    end
    % Plot enrichment
    figure(arcft); clf
    set(gcf,'Units','centimeters','Position',[5 5 10 3])
    hold on
    for arch = 1:Numarchs
        if arch==enrArcNum
            lw = 2.2;
            clr = ARCHETYPE_PALETTE(arch,:);
        else
            lw = 0.8;
            clr = [0.7 0.7 0.7];
        end
        plot(1:numOfBins,PoverQ_continuous(arch,:),'LineWidth',lw,'Color',clr)
    end
    xlabel('# Bin')
    ylabel('Enrichment')
    legend(strcat('Arc #' ,num2str((1:Numarchs)')), 'Location','eastoutside')
    set(gca,'fontsize',8)
    exportgraphics(figure(arcft),[celltype_folder filesep 'feature_plots' filesep ...
        enr_feat '_' num2str(enrArcNum) '.png'])
    close(arcft)
    enrich_tbl_cont.PoverQ(arcft) = PoverQ_continuous(enrArcNum,1);
    diffPoverQ = diff(PoverQ_continuous(enrArcNum,:));
    enrich_tbl_cont.StableTopBins(arcft) = all(diffPoverQ(1:2)<=0);
    
    % Plot case if "real" and run Parti_lite without the positive cells in the first bin:
    if enrich_tbl_cont.MaxTissueAbundance(arcft) < spec_th &&...
            enrich_tbl_cont.MaxDonorAbundance(arcft) < spec_th &&...
            enrich_tbl_cont.DonorsWithFeature(arcft) > 1 &&...
            enrich_tbl_cont.PoverQ(arcft) > PoQ_th_cont &&... 
            enrich_tbl_cont.StableTopBins(arcft)
        enrich_tbl_cont.PossibleEnrichment(arcft) = 1;
        DataPointsIndBinned = mat2cell(DataPointsInd(enrArcNum,:)',numPointInBin, numFeatures);
        pos_feat = find(EnMatCont>Q);
        neg_feat = find(EnMatCont<=Q);
        figure(100+arcft); clf;
        hold on
        scatter3(score(pos_feat,1),score(pos_feat,2),score(pos_feat,3),...
            3,'r','filled')
        scatter3(score(neg_feat,1),score(neg_feat,2),score(neg_feat,3),...
            3,'k','filled')
        pos_feat_bin1 = intersect(pos_feat,DataPointsIndBinned{1});
        neg_feat_bin1 = intersect(neg_feat,DataPointsIndBinned{1});
        scatter3(score(pos_feat_bin1,1),score(pos_feat_bin1,2),score(pos_feat_bin1,3),...
            12,'r','filled')
        scatter3(score(neg_feat_bin1,1),score(neg_feat_bin1,2),score(neg_feat_bin1,3),...
            12,'k','filled')
        plot3(arc(:,1),arc(:,2),arc(:,3),'d','MarkerSize',4,...
            'MarkerFaceColor','k','MarkerEdgeColor','k')        
        xlabel('PC 1')
        ylabel('PC 2')
        zlabel('PC 3')
        shpParti = alphaShape(arc(:,1),arc(:,2),arc(:,3));
        plot(shpParti,'FaceAlpha',0.1,'FaceColor',[0.5 0.5 0.5])
        text(arc(:,1)+3,arc(:,2)+3,arc(:,3)+3,num2str((1:size(arc,1))'),"FontSize",12);
        title([enr_feat ' Enriched at ' num2str(enrArcNum)]);
        savefig([celltype_folder filesep 'feature_plots' filesep ...
        'PCA_' enr_feat '_' num2str(enrArcNum) '.fig'])

    else
        enrich_tbl_cont.PossibleEnrichment(arcft) = 0;
    end
end

enrich_tbl_disc = table();
arc_enrich = zeros([size(sig_disc_tbl,1) size(arc)]);
for arcft = 1:size(sig_disc_tbl,1)
    enr_feat_full = sig_disc_tbl.FeatureName{arcft};
    if contains(enr_feat_full,':')
        enr_feat = strip(enr_feat_full(1:strfind(enr_feat_full,':')-1));
        enr_direction = strip(enr_feat_full(strfind(enr_feat_full,':')+1:end));
    end
    enrArcNum = sig_disc_tbl.archetype_(arcft);
    switch enr_feat
        case 'NA'
            continue
    end
    switch enr_direction
        case '1'
            EnMatDis = meta_discrete(:,strcmp(enr_feat,titles_discrete)); % metadata values
        case '0'
            EnMatDis = 1-meta_discrete(:,strcmp(enr_feat,titles_discrete)); % metadata values
        case '-1'
            continue
    end
    % Divide metadata to bins
    tempEnrich =  cumsum(EnMatDis(DataPointsInd(enrArcNum,:),:));
    binnedEnrichment = diff([zeros(1,numFeatures) ; tempEnrich(breakPoints,:) ]);
    DataPointsIndBinned = mat2cell(DataPointsInd(enrArcNum,:)',numPointInBin, numFeatures);
    [Aid,Bid] = groupcounts(donor_per_cell(DataPointsIndBinned{1}(EnMatDis(DataPointsIndBinned{1})==1)));
    [Ats,Bts] = groupcounts(meta_tbl.tissue_in_publication(DataPointsIndBinned{1}(EnMatDis(DataPointsIndBinned{1})==1)));
    enrich_tbl_disc.Arc(arcft) = enrArcNum;
    enrich_tbl_disc.Feature{arcft} = enr_feat_full;
    enrich_tbl_disc.MaxTissueAbundance(arcft) = max(Ats./sum(Ats));
    enrich_tbl_disc.MaxDonorAbundance(arcft) = max(Aid./sum(Aid));
    enrich_tbl_disc.DonorsWithFeature(arcft) = numel(unique(donor_per_cell(EnMatDis==1)));
    PoverQ = zeros(Numarchs,numOfBins);
    for arch = 1:Numarchs
        tempEnrich =  cumsum(EnMatDis(DataPointsInd(arch,:),:));
        binnedEnrichment = diff([zeros(1,numFeatures) ; tempEnrich(breakPoints,:) ]);
        DataPointsIndBinned = mat2cell(DataPointsInd(arch,:)',numPointInBin, numFeatures);
        BinnedDist = mat2cell(distances(arch,:)',numPointInBin, numFeatures);
        
        %EnPQ - Prepare enrichment graphs
        P = (bsxfun(@rdivide,(binnedEnrichment'),numPointInBin))'; %%round(mean(numPointInBin)));
        Q = tempEnrich(end,:)./numDataPoints;
        PoverQ(arch,:) = bsxfun(@rdivide, P,Q);
    end
    % Plot enrichment
    figure(arcft); clf
    set(gcf,'Units','centimeters','Position',[5 5 10 3])
    hold on
    for arch = 1:Numarchs
        if arch==enrArcNum
            lw = 2.2;
            clr = ARCHETYPE_PALETTE(arch,:);
        else
            lw = 0.8;
            clr = [0.7 0.7 0.7];
        end
        plot(1:numOfBins,PoverQ(arch,:),'LineWidth',lw,'Color',clr)
    end
    xlabel('# Bin')
    ylabel('Enrichment')
    legend(strcat('Arc #' ,num2str((1:Numarchs)')), 'Location','eastoutside')
    set(gca,'fontsize',8)
    exportgraphics(figure(arcft),[celltype_folder filesep 'feature_plots' filesep ...
        replace(erase(enr_feat_full,' '),':','-') '_' num2str(enrArcNum) '.png'])
    close(arcft)
    enrich_tbl_disc.PoverQ(arcft) = PoverQ(enrArcNum,1);
    diffPoverQ = diff(PoverQ(enrArcNum,:));
    enrich_tbl_disc.StableTopBins(arcft) = all(diffPoverQ(1:2)<=0);
    % plot case if "real":
    if enrich_tbl_disc.MaxTissueAbundance(arcft) < spec_th &&...
            enrich_tbl_disc.MaxDonorAbundance(arcft) < spec_th &&...
            enrich_tbl_disc.DonorsWithFeature(arcft) > 1 &&...
            enrich_tbl_disc.PoverQ(arcft) > PoQ_th_disc &&... 
            enrich_tbl_disc.StableTopBins(arcft)
        enrich_tbl_disc.PossibleEnrichment(arcft) = 1;
        DataPointsIndBinned = mat2cell(DataPointsInd(enrArcNum,:)',numPointInBin, numFeatures);
        pos_feat = find(EnMatDis==1);
        neg_feat = find(EnMatDis==0);
        figure(100+arcft); clf;
        set(gcf,'units','centimeters','position',[3 3 20 15])
        hold on
        scatter3(score(pos_feat,1),score(pos_feat,2),score(pos_feat,3),...
            3,'r','filled')
        scatter3(score(neg_feat,1),score(neg_feat,2),score(neg_feat,3),...
            3,'k','filled')
        pos_feat_bin1 = intersect(pos_feat,DataPointsIndBinned{1});
        neg_feat_bin1 = intersect(neg_feat,DataPointsIndBinned{1});
        scatter3(score(pos_feat_bin1,1),score(pos_feat_bin1,2),score(pos_feat_bin1,3),...
            12,'r','filled')
        scatter3(score(neg_feat_bin1,1),score(neg_feat_bin1,2),score(neg_feat_bin1,3),...
            12,'k','filled')
        plot3(arc(:,1),arc(:,2),arc(:,3),'d','MarkerSize',4,...
            'MarkerFaceColor','k','MarkerEdgeColor','k')
        xlabel('PC 1')
        ylabel('PC 2')
        zlabel('PC 3')
        shpParti = alphaShape(arc(:,1),arc(:,2),arc(:,3));
        plot(shpParti,'FaceAlpha',0.1,'FaceColor',[0.5 0.5 0.5])
        text(arc(:,1)+3,arc(:,2)+3,arc(:,3)+3,num2str((1:size(arc,1))'),"FontSize",12);
        title([replace(enr_feat_full,'Anoxi','Anoxia') ' Enriched at ' num2str(enrArcNum)]);
        set(gca,'fontsize',16)
        savefig([celltype_folder filesep 'feature_plots' filesep ...
        'PCA_' replace(erase(enr_feat_full,' '),':','-') '_' num2str(enrArcNum) '.fig'])
        
        % % Remove enriched cells in first bin:
        % geneExpression_sub = geneExpression;
        % geneExpression_sub(pos_feat_bin1,:) = [];
        % % Run parti:
        % OutputFileNamePrefix = [celltype_folder filesep 'parti' filesep 'ParTi_res_' replace(erase(enr_feat_full,' '),':','-') '_' num2str(enrArcNum) ];
        % if ~exist(strcat(OutputFileNamePrefix,'.mat'),'file')
        %     fprintf('Arc %s Running ParTi Lite on %s...\n',num2str(enrArcNum),enr_feat_full)
        %     % arc_temp = ParTI_lite(geneExpression_sub,5,...
        %     %     8,[],[],0,geneNames,geneExpression_sub,[],0.05,ParTiLiteFileName_csv,[],size(arc,1));
        %     arc_enrich = ParTI_lite(geneExpression_sub,5,...
        %         8,[],[],0,[],[],[],0.05,[],[],size(arc,1));
        %      [R2_global,R2_per_archetype]  = archetype_R2(arc,arc_enrich);
        %     save(strcat(OutputFileNamePrefix,'.mat'),...
        %         'arc_enrich','R2_global','R2_per_archetype')
        %     clear arc_enrich R2_global R2_per_archetype
        % end
    else
        enrich_tbl_disc.PossibleEnrichment(arcft) = 0;

    end
end
enrich_tbl_disc(enrich_tbl_disc.Arc==0,:) = [];
enrich_tbl = [ enrich_tbl_cont; enrich_tbl_disc];
writetable(enrich_tbl,[celltype_folder filesep celltype '_enrichment_table.xlsx'])


