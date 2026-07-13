% Download the three cell type folder with data from https://drive.google.com/drive/folders/1qMSm7LHTBpLo_I2-KqWs5KjF9H9ox4ij?usp=drive_link
% and change the pathes below, including the savefig_path folder: 
%% Gene enrichment tables and gene name tables used for ParTI:
celltypes = {'Macro','Fibro','Endo'};
human_weight_th = 0.5;
mouse_weight_th = 0.5;

% Path to gene enrichment data:
savefig_path = 'C:\Users\amorp\Documents\Work\Adler\Projects\Supportive_Archetypes\Figure_human_mouse_cmp';
geneEnrich_Hu_tblFile{1} = [savefig_path filesep 'Macrophages\human_Mac_rev_muscle_stomech_k5_continuous_significant.csv'];
geneEnrich_Ms_tblFile{1} = [savefig_path filesep  'Macrophages\mice_macro_5_continuous_significant.csv'];
geneEnrich_Hu_tblFile{2} = [savefig_path filesep  'Fibroblasts\human_ParTI_lite_gene_enrichment_results_continuous_significant.csv'];
geneEnrich_Ms_tblFile{2} = [savefig_path filesep  'Fibroblasts\mouse_ParTI_lite_gene_enrichment_results_continuous_significant.csv'];
geneEnrich_Hu_tblFile{3} = [savefig_path filesep  'EndothelialCells\human_ParTi_Lite_Out_Endo_6_continuous_significant.csv'];
geneEnrich_Ms_tblFile{3} = [savefig_path filesep  'EndothelialCells\mouse_ParTi_Lite_Out_EndoMs_continuous_significant.csv'];
% Path to geneNames:
geneNames_Hu_matFile{1} = [savefig_path filesep  'Macrophages\human_geneNames.mat'];
geneNames_Ms_matFile{1} = [savefig_path filesep  'Macrophages\mouse_geneNames.mat'];
geneNames_Hu_matFile{2} = [savefig_path filesep  'Fibroblasts\human_geneNames.mat'];
geneNames_Ms_matFile{2} = [savefig_path filesep  'Fibroblasts\mouse_geneNames.mat'];
geneNames_Hu_matFile{3} = [savefig_path filesep  'EndothelialCells\human_geneNames.mat'];
geneNames_Ms_matFile{3} = [savefig_path filesep  'EndothelialCells\mouse_geneNames.mat'];
% Path to weights:
weights_Hu_tblFile{1} = [savefig_path filesep 'Macrophages\human_composition_sum_to_arc_AM_'...
    num2str(human_weight_th) '.csv'];
weights_Ms_tblFile{1} =  [savefig_path filesep 'Macrophages\AM_sum_to_arc_mice.csv'];
weights_Hu_tblFile{2} =  [savefig_path filesep 'Fibroblasts\human_composition_sum_to_arc_AF_'...
    num2str(human_weight_th) '.csv'];
weights_Ms_tblFile{2} =  [savefig_path filesep 'Fibroblasts\AF_sum_to_arc_mice.csv'];
weights_Hu_tblFile{3} =  [savefig_path filesep 'EndothelialCells\human_composition_sum_to_arc_AE_'...
    num2str(human_weight_th) '.csv'];
weights_Ms_tblFile{3} =  [savefig_path filesep 'EndothelialCells\AE_sum_to_arc_mice.csv'];

% Figure data table with archetype functions: 
archetype_function_table = [savefig_path filesep 'Compare_Human_Mouse_Figure_Data.xlsx'];
%% Create save locations:
if ~exist(savefig_path,'dir')
    mkdir(savefig_path);
end
% Params:
nSim = 10000; % number of gene permutations for p-val
pval_th = 0.05; % for significance
gini_th = 0.31; % universal / tissue specific threshold
sheetNames = {'gene_match','match_pval','Ms_weight_03','Hu_weight_03','Ms_archetype_function','Hu_archetype_function'};
tissue_legend = {'Bladder','Fat','Heart',...
    'Kidney','Large intestine','Liver',...
    'Lung','Pancreas','Skin',...
    'Small intestine','Spleen','Thymus',...
    'Tongue','Trachea','Marrow',...
    'Limb muscle','Mammary gland',...
    'Muscle','Stomach'};
tissue_legend_clr = {'#3953A4','#ED2024','#6ABD45',...
    '#0F1031','#E03F97','#1F5429',...
    '#FED304','#4996D2','#9A4D42',...
    '#DF7D26','#714EA0','#1E9698',...
    '#E3AFD1','#A93493','#A781BA',...
    '#18BDC2','#C29A2D',...
    '#5C8A00','#C49A6C'}; 

%% Read full gene lists used for ParTi 
Hu_Ms_cmp = struct();
for cl = 1:numel(celltypes)
    fprintf('%s\n',celltypes{cl})
    Hu_Ms_cmp(cl).celltypeName = celltypes{cl} ;
    geneNamesHu = load(geneNames_Hu_matFile{cl},'geneNames');
    geneNamesMs = load(geneNames_Ms_matFile{cl},'geneNames');
    % Find common genes:
    compare_name = ['Ms_Hu_' celltypes{cl}];
    geneNamesHuLower = reshape(lower(geneNamesHu.geneNames),[],1);
    geneNamesMsLower = reshape(lower(geneNamesMs.geneNames),[],1);
    jointGenes = geneNamesHuLower(ismember(geneNamesHuLower,geneNamesMsLower));
    % Remove genes not in jointGenes from archetype table:
    tableHu = readtable(geneEnrich_Hu_tblFile{cl});
    tableMs = readtable(geneEnrich_Ms_tblFile{cl});
    tableHu.FeatureNameLower = lower(tableHu.FeatureName);
    tableMs.FeatureNameLower = lower(tableMs.FeatureName);
    tableHu(~ismember(tableHu.FeatureNameLower,jointGenes),:) = [];
    tableMs(~ismember(tableMs.FeatureNameLower,jointGenes),:) = [];
    % Compare matched genes per archetype
    arcsHu = unique(tableHu.archetype_);
    arcsMs = unique(tableMs.archetype_);
    genesMat = zeros(numel(arcsHu),numel(arcsMs));
    for arc1 = 1:numel(arcsHu)
        genes1 = lower(tableHu.FeatureName(tableHu.archetype_==arc1));
        for arc2 = 1:numel(arcsMs)
            genes2 = lower(tableMs.FeatureName(tableMs.archetype_==arc2));
            genesMat(arc1,arc2) = 100*sum(ismember(genes1,genes2))/numel(unique([genes1; genes2]));
        end
    end
    Hu_Ms_cmp(cl).gene_match = genesMat;
    % Plot
    figure(10+cl); clf; 
    imagesc(genesMat)
    xlabel('Mouse Archetypes')
    ylabel('Human Archetypes')
    clim([0 max(genesMat(:))])
    colorbar
    text(reshape(repmat(1:size(genesMat,2),[size(genesMat,1),1]),1,[])-0.25,repmat(1:size(genesMat,1),[1,size(genesMat,2)]),num2str(round(genesMat(:),2)))
    
    % Calculate significance
    if ~exist([savefig_path filesep 'Sim_' Hu_Ms_cmp(cl).celltypeName '.mat'],'file')
        genesMatRand = zeros(nSim,numel(arcsHu),numel(arcsMs));
        fprintf('Running permuations...')
        for n = 1:nSim
            randArcs = randperm(size(tableMs,1));
            arc2_rand = tableMs.archetype_(randArcs);
            for arc1 = 1:numel(arcsHu)
                genes1 = lower(tableHu.FeatureName(tableHu.archetype_==arc1));
                for arc2 = 1:numel(arcsMs)
                    genes2 = lower(tableMs.FeatureName(arc2_rand==arc2));
                    genesMatRand(n,arc1,arc2) =  100*sum(ismember(genes1,genes2))/numel(unique([genes1; genes2]));
                end
            end
        end
        save([savefig_path filesep 'Sim_' Hu_Ms_cmp(cl).celltypeName '.mat'],'genesMatRand')
        fprintf('Done\n')
    else
        fprintf('Loading permuations.\n')
        load([savefig_path filesep 'Sim_' Hu_Ms_cmp(cl).celltypeName '.mat'],'genesMatRand')
    end

    figure(20+cl); clf
    set(gcf,'Units','normalized','Position',[0.05 0.05 0.9 0.9])
    c = 1;
    pvals = zeros(numel(arcsHu),numel(arcsMs));
    for a1 = 1:numel(arcsHu)
        for a2 = 1:numel(arcsMs)
            subplot(numel(arcsHu),numel(arcsMs),c)
            histogram(squeeze(genesMatRand(:,a1,a2)),15)
            hold on
            xline(genesMat(a1,a2),'r','LineWidth',2);
            pvals(a1,a2) = sum(genesMatRand(:,a1,a2)>genesMat(a1,a2))/nSim;
            title(['p=' num2str(round(pvals(a1,a2),4))])
            c = c+1;
        end
    end
    Hu_Ms_cmp(cl).match_pval = pvals;
    % Read weights:
    Hu_Ms_cmp(cl).Hu_weight = readtable(weights_Hu_tblFile{cl});
    Hu_Ms_cmp(cl).Ms_weight = readtable(weights_Ms_tblFile{cl});
    Hu_Ms_cmp(cl).Mouse_weight_th = mouse_weight_th;
    % Read archetype functions:
    Hu_Ms_cmp(cl).Hu_archetype_function = readtable(archetype_function_table,'Sheet',[celltypes{cl} '_Hu_archetype_function']);
    Hu_Ms_cmp(cl).Ms_archetype_function = readtable(archetype_function_table,'Sheet',[celltypes{cl} '_Ms_archetype_function']);
end


%% Gini index for universal / tissue specific archetypes:
for cl = 1:numel(celltypes)
    arcNumsH = unique(Hu_Ms_cmp(cl).Hu_weight.archetype);
    arc_compH = table2array(Hu_Ms_cmp(cl).Hu_weight(:,2:end));
    tissuesH = Hu_Ms_cmp(cl).Hu_weight.Properties.VariableNames(2:end);
    for hi = numel(arcNumsH):-1:1
        arc_row = find(Hu_Ms_cmp(cl).Hu_weight.archetype==arcNumsH(hi));
        % Hu_Ms_cmp(cl).giniH(hi) = 1 - sum((Hu_Ms_cmp(cl).Hu_weight.Percentage(Hu_Ms_cmp(cl).Hu_weight.archetype==arcNumsH(hi))./100).^2);
        Hu_Ms_cmp(cl).giniH(hi) = 1 - sum((arc_compH(arc_row,:)./100).^2);
        % [~,maxTissueInd] = max(Hu_Ms_cmp(cl).Hu_weight.Percentage(Hu_Ms_cmp(cl).Hu_weight.archetype==arcNumsH(hi)));
        [~,maxTissueInd] = max(arc_compH(arc_row,:));
        % tissuesH = Hu_Ms_cmp(cl).Hu_weight.Tissue(Hu_Ms_cmp(cl).Hu_weight.archetype==arcNumsH(hi));
        Hu_Ms_cmp(cl).maxTissueH{hi} = tissuesH{maxTissueInd};
    end
    
    arcNumsM = unique(Hu_Ms_cmp(cl).Ms_weight.archetype); % includes "0" -> "generalists"
    arc_compM = table2array(Hu_Ms_cmp(cl).Ms_weight(:,2:end));
    tissuesM = Hu_Ms_cmp(cl).Ms_weight.Properties.VariableNames(2:end);

    for mi = numel(arcNumsM):-1:1
        arc_row = find(Hu_Ms_cmp(cl).Ms_weight.archetype==arcNumsM(mi));
        Hu_Ms_cmp(cl).giniM(mi) = 1 - sum((arc_compM(arc_row,:)./100).^2);

        % Hu_Ms_cmp(cl).giniM(mi) = 1 - sum((Hu_Ms_cmp(cl).Ms_weight.Percentage(Hu_Ms_cmp(cl).Ms_weight.archetype==arcNumsM(mi))./100).^2);
        % [~,maxTissueInd] = max(Hu_Ms_cmp(cl).Ms_weight.Percentage(Hu_Ms_cmp(cl).Ms_weight.archetype==arcNumsM(mi)));
        [~,maxTissueInd] = max(arc_compM(arc_row,:));
        % tissuesM = Hu_Ms_cmp(cl).Ms_weight.Tissue(Hu_Ms_cmp(cl).Ms_weight.archetype==arcNumsM(mi));
        Hu_Ms_cmp(cl).maxTissueM{mi} = tissuesM{maxTissueInd};
    end
end

%% Plot archetypes
uniColor = [0.9 0.9 0.9];
% l = 1;
for cl = 1:numel(celltypes)
    match_lw =[];
    l = 1;
    matchTbl =  Hu_Ms_cmp(cl).gene_match;
    if all(matchTbl(:)<=1)
        matchTbl = matchTbl*100;
    end
    pvalTbl = Hu_Ms_cmp(cl).match_pval;
    [sigPvalsH,sigPvalsM] = ind2sub(size(pvalTbl),find(pvalTbl<pval_th));
    rect_bot = 0.05;
    rect_top = 0.9;
    rect_h = 0.12;
    rect_w = 0.322;
    rectX_h = 0;
    rectX_m = 0.65;
    rect_arch_vec = linspace(rect_bot,rect_top-rect_h,size(matchTbl,1)); 
    rect_arcm_vec = linspace(rect_bot,rect_top-rect_h,size(matchTbl,2)); 
    curv = 0.1;
    lw = 1.5;
    alpha_val = 0.4;
    figure(cl); clf;
    set(gcf,'units','normalized','Position',[0.01 0.1 0.8 0.8],'Color','white',...
        'Name',Hu_Ms_cmp(cl).celltypeName);
    ax = axes('units','normalized','Position',[0.01 0.1 0.95 0.9]);
    hold on
    set(ax,'Visible','off','Box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[]);
    hold on
    % Plot human archetypes
    for harcs = 1:size(matchTbl,1)
        funcs = split(Hu_Ms_cmp(cl).Hu_archetype_function.(['Hu' num2str(harcs)]),',');
        for fnc = 1:size(funcs,1)
            currFun = strtrim(funcs{fnc});
            % currFun(1) = upper(currFun(1));
            text(ax,rectX_h+rect_w/50,rect_arch_vec(harcs)+...
                rect_h/(size(funcs,1))*(fnc-0.5),strtrim(funcs{fnc}))
        end
        if harcs==1
            text(ax,rectX_h+rect_w/3,rect_arch_vec(harcs)-0.3*rect_h,'Human Archetypes')
        end
        if Hu_Ms_cmp(cl).giniH(harcs)>=gini_th
            rect_text = 'Universal';
            rect_clr = uniColor;
        else
            rect_text = Hu_Ms_cmp(cl).maxTissueH{harcs};
            rect_clr = hex2rgb([tissue_legend_clr{strcmpi(Hu_Ms_cmp(cl).maxTissueH{harcs},tissue_legend)}]);
        end

        rectangle(ax,'Position',[rectX_h rect_arch_vec(harcs) rect_w rect_h],...
            'Curvature',curv,'LineWidth',lw,'FaceColor',rect_clr,'FaceAlpha',alpha_val);
        text(ax,rectX_h+rect_w+rect_w*0.02,rect_arch_vec(harcs)+rect_h-rect_h*0.12,rect_text)

    end
    xlim([0 1]);
    ylim([0 1]);
    
    % Plot mouse archetypes
    for marcs = 1:size(matchTbl,2)
        funcs = split(Hu_Ms_cmp(cl).Ms_archetype_function.(['M' num2str(marcs)]),', ');
        for fnc = 1:size(funcs,1)
            currFun = strtrim(funcs{fnc});
             % currFun(1) = upper(currFun(1));
            text(ax,rectX_m+rect_w/50,rect_arcm_vec(marcs)+rect_h/(size(funcs,1))*(fnc-0.5),...
               currFun)
        end
        if marcs==1
            text(ax,rectX_m+rect_w/3,rect_arcm_vec(marcs)-0.3*rect_h,'Mouse Archetypes')
        end
        if Hu_Ms_cmp(cl).giniM(marcs)>=gini_th
            rect_text = 'Universal';
            rect_clr = uniColor;
        else
            rect_text = Hu_Ms_cmp(cl).maxTissueM{marcs};
            rect_text = replace(rect_text,rect_text(1),upper(rect_text(1)));
            rect_clr = hex2rgb([tissue_legend_clr{strcmpi(Hu_Ms_cmp(cl).maxTissueM{marcs},tissue_legend)}]);
        end
        rectangle(ax,'Position',[rectX_m rect_arcm_vec(marcs) rect_w rect_h],...
            'Curvature',curv,'LineWidth',lw,'FaceColor',rect_clr,'FaceAlpha',alpha_val);
        text(ax,rectX_m+rect_w+rect_w*0.02,rect_arcm_vec(marcs)+rect_h-rect_h*0.12,...
            rect_text)

    end
    
    % Plot connecting lines:
    linecolors = brewermap(256,'Greys');
    % linecolors(1:50,:) = []; % pval-grey
    % linecolors(130:end,:) = []; % pval-grey
    % linecolors(1:100,:) = []; % rand-grey
    % linecolors(41:end,:) = []; % rand-grey
    linecolors(1:70,:) = []; % rand-grey
    linecolors(100:end,:) = []; % rand-grey
    xq1 = rectX_h+rect_w:0.001:rectX_m;
    for lines = 1:numel(sigPvalsH)
        x = linspace(rectX_h+rect_w,rectX_m,7);
        startArc = rect_arch_vec(sigPvalsH(lines));
        endArc = rect_arcm_vec(sigPvalsM(lines));
        y = [startArc startArc startArc startArc+(endArc-startArc)/2 ...
            endArc endArc endArc];
        linefun = makima(x,y,xq1)+rect_h/2;
        if sum(ismember(sigPvalsH,sigPvalsH(lines)))>1
            y_delta = 1;
        else
            y_delta = 0;
        end
        match_lw(l) = round(matchTbl(sigPvalsH(lines),sigPvalsM(lines)));
        % pval_lc(l) = round((abs(log10(pvalTbl(sigPvalsH(lines),sigPvalsM(lines))+1/nSim))-...
        %     abs(log10(pval_th)))*...
        %     size(linecolors,1)/(abs(log10(1/nSim))-abs(log10(pval_th))));
        plot(ax,xq1,linefun+y_delta*randn*0.005,...
            'LineWidth',match_lw(l),'Color',linecolors(randi(size(linecolors,1)),:));  % rand-grey
        % plot(ax,xq1,linefun+y_delta*randn*0.005,...
        %     'LineWidth',match_lw(l),'Color',linecolors(pval_lc(l),:)); % pval-grey
        l = l+1;
    end
    
    set(gca,'Visible','off','YDir','reverse')
    
    lineScaleBarVec = min(match_lw):2:max(match_lw);
    axScale = axes('units','normalized','Position',[0.31 0.01 0.32 0.1]);
    hold on
    for ln = 1:numel(lineScaleBarVec)
        plot(axScale,[ln ln+1],[1 1],...
                'LineWidth',lineScaleBarVec(ln),'Color',linecolors(round(size(linecolors,1)/2),:)); 
        text(axScale,ln+0.2,0.1,[num2str(lineScaleBarVec(ln)) '%'],'fontsize',10)
    end
    text(axScale,numel(lineScaleBarVec)/2-1,2,'Enriched gene overlap')
    set(axScale,'Visible','off','Box','off');
    exportgraphics(figure(cl),[savefig_path filesep Hu_Ms_cmp(cl).celltypeName '_' num2str(mouse_weight_th) '.pdf'],...
        'ContentType','vector'); % rand-grey
    % exportgraphics(figure(cl),[savefig_path filesep Hu_Ms_cmp(cl).celltypeName '_greyPval.pdf'],'ContentType','vector'); % pval-grey
end

%% Plot heatmaps for supplementary: 
clrmap = brewermap(100,'Oranges');
sq_sz = 1.5; %cm
for cl = 1:numel(celltypes)
    matchTbl = Hu_Ms_cmp(cl).gene_match;
    if any(matchTbl(:)>1)
        matchTbl = matchTbl/100;
    end
    pValTbl = Hu_Ms_cmp(cl).match_pval;
    figure(10+cl); clf
    set(gcf,'units','normalized','Position',[0.1 0.1 0.5 0.6],'Name',['Heatmap ' Hu_Ms_cmp(cl).celltypeName]);
    ax = axes('Units','centimeters','Position',[2 2 sq_sz*size(matchTbl,2)+2 sq_sz*size(matchTbl,1)]);
    imagesc(ax,matchTbl*100);
    xlabel('Mouse archetypes')
    ylabel('Human archetypes')
    xlim(0.5+[0 size(matchTbl,2)])
    xticks(1:size(matchTbl,2))
    ylim(0.5+[0 size(matchTbl,1)])
    yticks(1:size(matchTbl,1))
    colormap(clrmap)
    cb = colorbar;
    cb.Label.String = '% Enriched gene overlap';
    clim([0 round(max(matchTbl(:))*100)])
    [i,j] = ind2sub(size(pValTbl),find(pValTbl<pval_th));
    for p = 1:numel(i)
        text(j-0.1,i+0.1,'*','FontSize',20)
    end
    exportgraphics(figure(10+cl),[savefig_path filesep 'heatmap_' Hu_Ms_cmp(cl).celltypeName '.pdf'],'ContentType','vector');

end

