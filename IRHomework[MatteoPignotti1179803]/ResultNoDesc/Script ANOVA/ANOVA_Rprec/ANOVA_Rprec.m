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

x = [0.4167, 0.0894, 0.459, 0.1884, 0.1333, 0., 0.337, 0.2941, 0.0357, 0.3576, 0.4444, 0.2308, 0.125, 0.4571, 0.6571, 0.4747, 0.164, 0.4262, 0.2308, 0.247, 0., 0.1837, 0.4848, 0.3725, 0.275, 0.0392, 0.2821, 0.0306, 0.0625, 0.2857, 0.0714, 0.5455, 0.0753, 0.2157, 0.3023, 0.0526, 0.1176, 0.0784, 0.0928, 0.1045, 0.2079, 0.4476, 0.1127, 0.0588, 0.216, 0.3051, 0.4444, 0.0828, 0.2157, 0.424, ];
y = [0.3958, 0.0894, 0.4672, 0.1994, 0.1556, 0., 0.3296, 0.2941, 0.0357, 0.351, 0.4444, 0.2308, 0.125, 0.4571, 0.6571, 0.4646, 0.164, 0.4098, 0.2308, 0.247, 0., 0.1837, 0.4848, 0.3725, 0.275, 0.0392, 0.3333, 0.0306, 0.0625, 0.2857, 0.0714, 0.5455, 0.0753, 0.2157, 0.3023, 0.0526, 0.1176, 0.0784, 0.0928, 0.1045, 0.3371, 0.4381, 0.0986, 0.0588, 0.2207, 0.3051, 0.4074, 0.069, 0.2059, 0.416, ];
z = [0.4167, 0.0854, 0.459, 0.1884, 0.1333, 0., 0.337, 0.2941, 0.0357, 0.3444, 0.4444, 0.2308, 0.125, 0.4571, 0.6571, 0.4545, 0.164, 0.4262, 0.2308, 0.25, 0., 0.1633, 0.4848, 0.3725, 0.275, 0.049, 0.2821, 0.0408, 0.0625, 0.2857, 0.0714, 0.5455, 0.0822, 0.2353, 0.3023, 0.0526, 0.1176, 0.0784, 0.0876, 0.1045, 0.2022, 0.4476, 0.1127, 0.0588, 0.216, 0.322, 0.4444, 0.0828, 0.2157, 0.424, ];
w = [0.3333, 0.0854, 0.3361, 0.2271, 0.1111, 0., 0.3556, 0.1961, 0.0357, 0.4901, 0.3333, 0.2308, 0.1875, 0.4857, 0.6571, 0.3232, 0.164, 0.4754, 0.2308, 0.2351, 0., 0.102, 0.4545, 0.3431, 0.275, 0.0588, 0.3333, 0.0306, 0.4375, 0.2857, 0.0714, 0.4545, 0.0822, 0.2353, 0.2442, 0.1053, 0.1176, 0.1176, 0.0619, 0.097, 0.1966, 0.2381, 0.169, 0.0588, 0.2207, 0.339, 0.2222, 0.0621, 0.1569, 0.384, ];

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