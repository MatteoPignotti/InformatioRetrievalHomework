%
% Copyright 2018-2019 University of Padua, Italy
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%     http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%
% Author: Nicola Ferro (ferro@dei.unipd.it)
%load('ap_C09.mat')

x = [1., 0.5, 0.3, 0.3, 0.2, 0., 1., 0.5, 0., 0.4, 0.5, 0.2, 0.1, 0.8, 1., 0.8, 0.3, 1., 0.2, 0., 0., 0.5, 0.6, 0.4, 0.7, 0.2, 0.5, 0., 0.1, 0.2, 0.1, 0.8, 0.2, 0.4, 0.9, 0.1, 0.2, 0., 0.5, 0.7, 0.2, 0.9, 0.5, 0., 0.1, 1., 0.9, 0., 0.4, 0.7];
y = [1., 0.4, 0.3, 0.3, 0.3, 0., 0.9, 0.5, 0., 0.5, 0.5, 0.2, 0.1, 0.8, 1., 0.6, 0.4, 0.9, 0.2, 0., 0., 0.5, 0.6, 0.5, 0.7, 0.1, 0.5, 0., 0.1, 0.2, 0.1, 0.8, 0.3, 0.4, 0.9, 0.1, 0.2, 0., 0.5, 0.7, 0.9, 0.9, 0.5, 0., 0., 1., 0.9, 0., 0.4, 0.8];
z = [1., 0.5, 0.3, 0.3, 0.2, 0., 1., 0.5, 0.1, 0.4, 0.5, 0.2, 0.1, 0.8, 1., 0.8, 0.3, 1., 0.2, 0., 0., 0.4, 0.6, 0.4, 0.7, 0.2, 0.5, 0., 0.1, 0.2, 0.1, 0.8, 0.2, 0.4, 0.9, 0.1, 0.2, 0., 0.5, 0.7, 0.2, 0.9, 0.5, 0.1, 0.1, 1., 0.9, 0., 0.4, 0.7];
w = [1., 0.4, 0.7, 0.4, 0.3, 0., 0.8, 0., 0.1, 0.9, 0.4, 0.2, 0.1, 0.8, 1., 0.6, 0.4, 0.8, 0.2, 0.1, 0., 0.4, 0.8, 0.4, 0.7, 0.2, 0.3, 0., 0.5, 0.2, 0.1, 0.6, 0.1, 0.3, 0.8, 0., 0.3, 0.1, 0.3, 0.8, 0.9, 0.3, 0.6, 0., 0., 1., 0.4, 0., 0.2, 0.8];

measure = [x.' y.' z.' w.'];
runID = ["BM25_0", "TFIDF_1", "BM25_2","TFIDF_3"];
runID = strtrim(cellstr(runID));
topicID = 351:400;
topicID = arrayfun(@num2str,topicID,'UniformOutput',false);

% the mean for each run across the topics
% Note that if the measure is AP (Average Precision), 
% this is exactly MAP (Mean Average Precision) for each run
m = mean(measure);

% sort in descending order of mean score
[~, idx] = sort(m, 'descend');

% re-order runs by ascending mean of the measure
% needed to have a more nice looking box plot
measure = measure(:, idx);
runID = runID(idx);

% perform the ANOVA
[~, tbl, sts] = anova1(measure, runID, 'off');

% display the ANOVA table
tbl

% perform
c = multcompare(sts, 'Alpha', 0.05, 'Ctype', 'hsd'); 

% display the multiple comparisons
c

%% plots of the data

% get the Tukey HSD test figure
currentFigure = gcf;

    ax = gca;
    ax.FontSize = 20;
    ax.XLabel.String = 'P10';
    ax.YLabel.String = 'Run';

    currentFigure.PaperPositionMode = 'auto';
    currentFigure.PaperUnits = 'centimeters';
    currentFigure.PaperSize = [42 22];
    currentFigure.PaperPosition = [1 1 40 20];

print(currentFigure, '-dpdf', 'ap-tukey.pdf');

% box plot
currentFigure = figure;
    % need to reverse the order of the columns to have bloxplot displayed
    % as the Tukey HSD plot
    boxplot(measure(:, end:-1:1), 'Labels', runID(end:-1:1), ...
        'Orientation', 'horizontal', 'Notch','off', 'Symbol', 'ro')
    
    ax = gca;
    ax.FontSize = 20;
    ax.XLabel.String = 'P10';
    ax.YLabel.String = 'Run';
    
    currentFigure.PaperPositionMode = 'auto';
    currentFigure.PaperUnits = 'centimeters';
    currentFigure.PaperSize = [42 22];
    currentFigure.PaperPosition = [1 1 40 20];

print(currentFigure, '-dpdf', 'ap-boxplot.pdf');