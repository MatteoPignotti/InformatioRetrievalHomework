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


x = [0.4181, 0.0373, 0.2185, 0.075, 0.0765, 0.0104, 0.2657, 0.2394, 0.0155, 0.2238, 0.3065, 0.0739, 0.0437, 0.4758, 0.6838, 0.3936, 0.0774, 0.4264, 0.2703, 0.0968, 0.0021, 0.137, 0.3783, 0.2292, 0.1793, 0.0137, 0.3113, 0.0047, 0.0704, 0.1909, 0.0303, 0.5199, 0.0235, 0.184, 0.2648, 0.0326, 0.1153, 0.0188, 0.038, 0.0824, 0.094, 0.4517, 0.0586, 0.0084, 0.1061, 0.2734, 0.4059, 0.0115, 0.1012, 0.3728];
y = [0.411, 0.0349, 0.2198, 0.0801, 0.0857, 0.0105, 0.2688, 0.2311, 0.0151, 0.2233, 0.3038, 0.0776, 0.0431, 0.4778, 0.6682, 0.4407, 0.0764, 0.3813, 0.2657, 0.0927, 0.0019, 0.1575, 0.3752, 0.2211, 0.1735, 0.0126, 0.3125, 0.0041, 0.0702, 0.1882, 0.0246, 0.5153, 0.0231, 0.1888, 0.2677, 0.0311, 0.1128, 0.0167, 0.0383, 0.082, 0.2819, 0.4444, 0.0555, 0.0065, 0.1044, 0.277, 0.4044, 0.0118, 0.0955, 0.3651];
z = [0.4163, 0.037, 0.2163, 0.0729, 0.0768, 0.0105, 0.2688, 0.2398, 0.0162, 0.2184, 0.3065, 0.0729, 0.0432, 0.4758, 0.6834, 0.3676, 0.0774, 0.4267, 0.2703, 0.0922, 0.0021, 0.1297, 0.3765, 0.2283, 0.1814, 0.0165, 0.3106, 0.0048, 0.0704, 0.1903, 0.0301, 0.5194, 0.0237, 0.1848, 0.2685, 0.0325, 0.1154, 0.0197, 0.0371, 0.0844, 0.0931, 0.4517, 0.0588, 0.0092, 0.1061, 0.2815, 0.4064, 0.0127, 0.101, 0.371];
w = [0.3349, 0.0347, 0.2221, 0.1026, 0.0702, 0.012, 0.2732, 0.1326, 0.0204, 0.391, 0.2564, 0.0906, 0.1594, 0.5314, 0.6964, 0.2415, 0.0767, 0.4257, 0.2657, 0.0868, 0.0022, 0.1008, 0.3986, 0.2812, 0.1737, 0.0178, 0.157, 0.0044, 0.3523, 0.146, 0.0258, 0.354, 0.0149, 0.1833, 0.2146, 0.0141, 0.1135, 0.0477, 0.0134, 0.078, 0.1531, 0.1124, 0.082, 0.014, 0.1044, 0.2792, 0.1975, 0.0076, 0.0589, 0.3379];

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
    ax.XLabel.String = 'Average Precision (AP)';
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
    ax.XLabel.String = 'Average Precision (AP)';
    ax.YLabel.String = 'Run';
    
    currentFigure.PaperPositionMode = 'auto';
    currentFigure.PaperUnits = 'centimeters';
    currentFigure.PaperSize = [42 22];
    currentFigure.PaperPosition = [1 1 40 20];

print(currentFigure, '-dpdf', 'ap-boxplot.pdf');