% Applying the Pareto Archetype Inference to the Macrophage scRNA-seq

% This implementation is based on the method and codebase provided
% by the Alon Lab at the Weizmann Institute, available at:
% https://github.com/AlonLabWIS/ParTI

% === Setup Paths ===
script_dir = fileparts(mfilename('fullpath'));
project_root = fullfile(script_dir, '..', '..');
addpath(genpath(fullfile(project_root, 'Code', 'particode')));

% === Locate CSV file ===
% Input file generated from Step 01 preprocessing
csv_filename = fullfile(project_root, 'Data', '02_Archetypes', 'Macrophage', 'Human_macro_forParti.csv');
if ~exist(csv_filename, 'file')
    error('CSV file not found: %s\nRun preprocessing step first.', csv_filename);
end

prefix = 'human_macro_5arc';
fprintf('Processing: %s\n', csv_filename);
fprintf('Output prefix: %s\n', prefix);

% === Start logging ===
timestamp = datetime('now', 'yyyymmdd_HHMM');
diary_filename = sprintf('%s_log_%s.txt', prefix, timestamp);
diary(diary_filename);

% === Load gene expression data from CSV ===
T = readtable(csv_filename, 'ReadRowNames', true);
geneNames = T.Properties.VariableNames;
geneExpression = table2array(T);   % [cells × genes]

fprintf('Loaded data: %d cells × %d genes\n', size(geneExpression, 1), size(geneExpression, 2));
fprintf('First 10 genes: ');
disp(geneNames(1:min(10, length(geneNames))));

% === Run ParTI ===
% This is the primary analysis step, executing the ParTI algorithm
%  to identify biological archetypes.
% Please note that this process takes several hours and is interactive:
%  the script will prompt the user at the start to select
%  the number of archetypes before proceeding (5 in this case)

% Run ParTI with 5 archetypes (k=5)
% The number of archetypes (5) was determined based on the elbow method
% and biological interpretability in the preliminary analysis.
[arc, arcOrig, pc, errs, pval] = ParTI( ...
    geneExpression, 5, 8, [], [], 0, [], [], [], 0.05, prefix);

save([prefix '_arc.mat'],     'arc');
save([prefix '_arcOrig.mat'], 'arcOrig');
save([prefix '_errs.mat'],    'errs');
save([prefix '_pval.mat'],    'pval');

% === Gene enrichment analysi ===

ParTI_lite(...
    geneExpression, 5, 8, [], [], 0, ...
    geneNames, geneExpression, [], 0.05, prefix, arcOrig);

% === Build GO matrix and run ParTI_lite on GO features ===
msigdb_path = fullfile(project_root, 'Code', 'particode', 'MSigDB');
if ~exist(msigdb_path, 'dir')
    warning('MSigDB directory not found at %s. GO analysis might fail.', msigdb_path);
end

[GOExpression, GONames, ~, GOcat2Genes] = MakeGOMatrix(...
    geneExpression, upper(geneNames), ...
    {fullfile(msigdb_path, 'c2.cp.v4.0.symbols.gmt'), fullfile(msigdb_path, 'c5.all.v4.0.symbols.gmt')}, ...
    20);

ParTI_lite(...
    geneExpression, 5, 8, [], [], 0, ...
    GONames, GOExpression, [], 0.05, [prefix '_GO'], arcOrig);

% === Finish ===
fprintf('\n=== Analysis Complete ===\n');
fprintf('Output prefix: %s\n', prefix);
fprintf('Check results in current directory\n');
diary off;
fprintf('Done.\n');
exit;
