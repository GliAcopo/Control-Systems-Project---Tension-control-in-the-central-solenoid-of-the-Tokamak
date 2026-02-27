function filename = SaveFigure(titleStr, figFolder, addTitle, addGrid)
% SaveFigure Saves the current figure (gcf) as a vector PDF.
%   filename = SaveFigure(Title, figFolder, addTitle, addGrid)
%
% - Title is used for the output filename.
% - If addTitle is true, Title is also applied to the current axes.
% - If addGrid is true, grid is turned on for the current axes.

    if nargin < 2 || isempty(figFolder)
        error('SaveFigure:MissingFolder', 'figFolder is required.');
    end
    if nargin < 3 || isempty(addTitle), addTitle = true; end
    if nargin < 4 || isempty(addGrid),  addGrid  = true; end

    % Ensure strings
    titleStr  = string(titleStr);
    figFolder = string(figFolder);

    % Ensure folder exists
    if ~isfolder(figFolder)
        disp("!!! WARNING: Folder"+string(figFolder)+"not found, creating folder")
        mkdir(figFolder);
    end

    fig = gcf;

    % Apply style to current axes (basic behavior)
    ax = gca;
    if addTitle
        title(ax, titleStr, 'Interpreter', 'none');
    end
    if addGrid
        grid(ax, 'on');
    end

    % Build a safe filename (Windows-safe too)
    safeName = strtrim(titleStr);
    safeName = regexprep(safeName, '[<>:"/\\|?*\x00-\x1F]', '_'); % invalid filename chars
    safeName = regexprep(safeName, '\s+', '_');                   % spaces -> underscore
    safeName = regexprep(safeName, '_+', '_');                    % collapse ___
    if strlength(safeName) == 0
        disp("!!! WARNING: string name invalid, defaulting to 'figure.pdf'")
        safeName = "figure";
    end

    filename = fullfile(figFolder, safeName + ".pdf");

    exportgraphics(fig, filename, 'ContentType', 'vector');

    fprintf("Saved: %s\n", filename);
end


% USAGE
%{
figure(1);
rlocus(L);
SaveFigure("Rlocus - Loop function of base system", figFolder, true, true);

figure(2);
step(L);
SaveFigure("Step response", figFolder, true, true);

SaveFigure("Same title but no grid", figFolder, true, false);
SaveFigure("Use only as filename", figFolder, false, true);

%}

%% BAD FUNCTION
%{
function SaveFigure(figNum, sys, titleStr, figFolder)
    % SAVE_RLOCUS_PDF  Create, style and export root locus as vector PDF
    %
    %   save_rlocus_pdf(1, L, "Rlocus - Base system", figFolder)
    %
    % Inputs:
    %   figNum     = figure number to use (1, 2, ...)
    %   sys        = the LTI system for rlocus (L in your case)
    %   titleStr   = title text (will also become filename base)
    %   figFolder  = folder path (your fullfile('..','Figures') thing)

    figure(figNum);
    clf;                       % clear anything old in this figure

    rlocus(sys);
    title(titleStr);
    grid on;

    % Build clean filename (replace forbidden / dangerous chars if needed)
    safeName = titleStr;
    safeName = regexprep(safeName, '[^a-zA-Z0-9 -:]', '_');  % keep letters, numbers, space, hyphen, colon
    safeName = strrep(safeName, ' ', '_');                   % optional: underscores instead of spaces
    filename = fullfile(figFolder, [safeName '.pdf']);

    exportgraphics(gcf, filename, 'ContentType', 'vector');

    fprintf('Saved: %s\n', filename);
end
%}

%% USAGE
%{
% Your main script
% Define relative path in wich to save file
figFolder = fullfile('..', 'Figures');

% One plot
SaveFigure(1, L, "Rlocus - Loop function of base system", figFolder);

% Another plot later
SaveFigure(2, L_new, "Rlocus - With compensator", figFolder);
%}


%{
function [outputArg1,outputArg2] = SaveFigure(inputArg1,inputArg2)
%SAVEFIGURE Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    inputArg1
    inputArg2
end

arguments (Output)
    outputArg1
    outputArg2
end

outputArg1 = inputArg1;
outputArg2 = inputArg2;
end
%}
%{
    tit = "Rlocus - Loop function of base system";
    title(tit);
    grid on;
    filename  = fullfile(figFolder, tit+".pdf");
    exportgraphics(gcf, filename, 'ContentType', 'vector');
%}