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

x = [0.5208, 0.0772, 0.4262, 0.2853, 0.2444, 0.0, 0.3926, 0.3137, 0.0714, 0.3576, 0.4444, 0.2051, 0.125, 0.5143, 0.6857, 0.4141, 0.2011, 0.4098, 0.3077, 0.2292, 0.1176, 0.3265, 0.4545, 0.3775, 0.325, 0.049, 0.3077, 0.1122, 0.125, 0.2857, 0.0714, 0.5909, 0.089, 0.2549, 0.2442, 0.0526, 0.1647, 0.3529, 0.0412, 0.2313, 0.3202, 0.4476, 0.0986, 0.1176, 0.2207, 0.322, 0.3704, 0.2069, 0.1667, 0.456];
y = [0.3958, 0.0, 0.418, 0.0693, 0.1556, 0.0, 0.0, 0.3137, 0.1071, 0.2384, 0.4444, 0.2051, 0.0625, 0.4857, 0.6571, 0.2929, 0.1217, 0.377, 0.1538, 0.003, 0.0588, 0.2857, 0.0, 0.3137, 0.175, 0.0, 0.2821, 0.0, 0.0625, 0.1429, 0.0357, 0.0455, 0.0959, 0.0, 0.2209, 0.0526, 0.1529, 0.3333, 0.0361, 0.209, 0.0, 0.0286, 0.0845, 0.0, 0.0516, 0.339, 0.3704, 0.0552, 0.1569, 0.416];
z = [0.5208, 0.0772, 0.4344, 0.2715, 0.2444, 0.0, 0.4, 0.3333, 0.1071, 0.3444, 0.4444, 0.2051, 0.125, 0.5143, 0.6857, 0.3939, 0.2063, 0.4098, 0.3077, 0.2083, 0.1176, 0.3265, 0.4545, 0.3775, 0.325, 0.0686, 0.3077, 0.1122, 0.125, 0.2857, 0.1071, 0.5909, 0.0822, 0.2353, 0.2674, 0.0526, 0.1882, 0.3333, 0.0412, 0.2388, 0.3427, 0.4381, 0.0986, 0.1176, 0.2254, 0.339, 0.3704, 0.2069, 0.1569, 0.456];
w = [0.4583, 0.0691, 0.4508, 0.2216, 0.1333, 0.0, 0.4074, 0.3137, 0.1429, 0.5033, 0.4444, 0.2308, 0.0, 0.6, 0.7143, 0.2323, 0.1693, 0.459, 0.2308, 0.1756, 0.1765, 0.2857, 0.4545, 0.3529, 0.325, 0.0392, 0.2821, 0.102, 0.3125, 0.1429, 0.0357, 0.4091, 0.089, 0.2353, 0.186, 0.1053, 0.1529, 0.2549, 0.0155, 0.1642, 0.4438, 0.2381, 0.1127, 0.1176, 0.2066, 0.3559, 0.1481, 0.1931, 0.098, 0.432];

measure = [x.' y.' z.' w.'];
runID = ["BM25_1", "BM25_2", "TF_IDF_0","TF_IDF_3"];
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
    ax.XLabel.String = 'Rprec';
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
    ax.XLabel.String = 'Rprec';
    ax.YLabel.String = 'Run';
    
    currentFigure.PaperPositionMode = 'auto';
    currentFigure.PaperUnits = 'centimeters';
    currentFigure.PaperSize = [42 22];
    currentFigure.PaperPosition = [1 1 40 20];

print(currentFigure, '-dpdf', 'ap-boxplot.pdf');