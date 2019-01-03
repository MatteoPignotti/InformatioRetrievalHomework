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

% set the random number generator for reproducibility
rng(6, 'twister');

% the number of samples to generate
n = 10;

% generate the data
x = [0.5229, 0.0237, 0.2362, 0.1677, 0.1332, 0.0106, 0.2956, 0.2997, 0.0462, 0.2238, 0.3065, 0.1122, 0.0454, 0.5543, 0.7150, 0.3440, 0.1248, 0.4079, 0.2746, 0.1254, 0.0574, 0.2488, 0.3859, 0.2389, 0.2589, 0.0052, 0.2968, 0.0401, 0.0725, 0.1695, 0.0569, 0.6056, 0.0257, 0.2380, 0.2177, 0.0223, 0.1382, 0.2662, 0.0145, 0.2022, 0.2668, 0.4335, 0.0479, 0.0235, 0.1290, 0.2880, 0.3125, 0.0975, 0.0800, 0.4169];
y = [0.5191, 0.0225, 0.2330, 0.1549, 0.1309, 0.0107, 0.3028, 0.3001, 0.0533, 0.2184, 0.3073, 0.1129, 0.0447, 0.5528, 0.7204, 0.3193, 0.1184, 0.4067, 0.2761, 0.1190, 0.0569, 0.2499, 0.3874, 0.2370, 0.2617, 0.0099, 0.2967, 0.0380, 0.0796, 0.1612, 0.0575, 0.5975, 0.0261, 0.2411, 0.2286, 0.0221, 0.1474, 0.2458, 0.0150, 0.1956, 0.2839, 0.4340, 0.0473, 0.0353, 0.1295, 0.3037, 0.3068, 0.1010, 0.0772, 0.4163];
z = [0.3805, 0.0000, 0.2238, 0.0080, 0.0716, 0.0060, 0.0000, 0.2865, 0.0325, 0.1742, 0.2827, 0.0858, 0.0546, 0.5145, 0.6599, 0.1596, 0.0498, 0.3362, 0.1608, 0.0001, 0.0141, 0.1843, 0.0050, 0.1847, 0.0817, 0.0004, 0.1814, 0.0000, 0.0514, 0.1467, 0.0360, 0.0345, 0.0215, 0.0002, 0.1614, 0.0133, 0.1174, 0.2232, 0.0134, 0.1628, 0.0000, 0.0165, 0.0271, 0.0025, 0.0065, 0.2988, 0.3290, 0.0062, 0.0591, 0.3582];
w = [0.4381, 0.0272, 0.3074, 0.1019, 0.0807, 0.0105, 0.2955, 0.2114, 0.0579, 0.4318, 0.2837, 0.0956, 0.0280, 0.5497, 0.7318, 0.1627, 0.1002, 0.4231, 0.2571, 0.0765, 0.0590, 0.2005, 0.3834, 0.2933, 0.2624, 0.0067, 0.1594, 0.0299, 0.2195, 0.1220, 0.0484, 0.4490, 0.0232, 0.2166, 0.1299, 0.0118, 0.1257, 0.1284, 0.0017, 0.1186, 0.4116, 0.1272, 0.0452, 0.0593, 0.1086, 0.3406, 0.1326, 0.0959, 0.0198, 0.3811];

data = [x.' y.' z.' w.'];

% the number of rows and columns
[p, q] = size(data);


% the estimated grand mean
hat_mu = mean(mean(data));
hat_mu

% the estimated column mean
hat_mu_col = mean(data);

% the predicted score, which is equal the the column mean,
% i.e. repeat the column mean for each row
hat_y = repmat(hat_mu_col, p, 1);

% the prediction error
hat_eps = data - hat_y;

% the total effects
ss_total = sum(sum( (data - hat_mu).^2 ));
df_total = p*q - 1;

% the column effects
% note that the column mean is the same for all the rows
% therefore we can multiply p (rows) times the inner summation
ss_col = p * sum( (hat_mu_col - hat_mu).^2 );
df_col = q - 1;
ms_col = ss_col ./ df_col;

% the error effects
ss_error = sum(sum( (data - repmat(hat_mu_col, p, 1)).^2  ));
df_error = q * (p - 1);
ms_error = ss_error ./ df_error;


% the significance level
alpha = 0.05;

% the F-statistic
Fstat = ms_col ./ ms_error;

% the F critical value under H0
Fcrit = finv(1 - alpha, df_col, df_error);

% the p-value
pval = (1 - fcdf(Fstat, df_col, df_error));

% fill in the ANOVA table
tbl = {'Source',    'SS',       'DF',       'MS',       'F',    'Prob>F'; ...
       'Columns',   ss_col,     df_col,     ms_col,     Fstat,  pval; ...
       'Error',     ss_error,   df_error,   ms_error,   [],     []; ...
       'Total',     ss_total,   df_total,   [],         [],     []}

   
% compute the p-value according to the Tukey HSD test for all the pairwise
% comparisons
t = abs(hat_mu_col(1) - hat_mu_col(2)) ./ sqrt(ms_error./p);
ptuk = internal.stats.stdrcdf(t, q*(p-1), q, 'upper');   
c(1, :) =  [1 2 ptuk];

t = abs(hat_mu_col(1) - hat_mu_col(3)) ./ sqrt(ms_error./p);
ptuk = internal.stats.stdrcdf(t, q*(p-1), q, 'upper');   
c(2, :) =  [1 3 ptuk];

t = abs(hat_mu_col(2) - hat_mu_col(3)) ./ sqrt(ms_error./p);
ptuk = internal.stats.stdrcdf(t, q*(p-1), q, 'upper');   
c(3, :) =  [2 3 ptuk];

% display the multiple comparisons
c


%% plots of the data

% box plot
currentFigure = figure;
    boxplot([x.' y.' z.' w.'], 'Labels', {'x', 'y', 'z','w'}, 'Notch','off', 'Symbol', 'ro')
    
    ax = gca;
    ax.FontSize = 20;

    currentFigure.PaperPositionMode = 'auto';
    currentFigure.PaperUnits = 'centimeters';
    currentFigure.PaperSize = [22 22];
    currentFigure.PaperPosition = [1 1 20 20];

print(currentFigure, '-dpdf', 'toy-boxplot.pdf');



% F distribution
a = 0:.01:6;
f = fpdf(a, df_col, df_error);

% find where we are above and below tCrit
idx = a >= Fcrit;

currentFigure = figure;

    % plot the F distribution
    plot(a, f, 'LineWidth', 4, 'Color', 'k')
    hold on
    
    % plot a vertical line corresponding to Fcrit
    h(1) = plot([Fcrit Fcrit], get(gca,'ylim'), 'Color', 'r', 'LineWidth', 3, 'LineStyle', '--');
    
    % color the area under the t distribution above Fcrit
    area(a(idx), f(idx), 'FaceColor', 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none')
    
    % plot a vertical lines corresponding to Fstat
    h(2) = plot([Fstat Fstat], get(gca,'ylim'), 'Color', 'b', 'LineWidth', 3, 'LineStyle', '--');
    
    
    legend(h, {'$F_{crit}$', '$F_{stat}$'}, 'Interpreter', 'Latex', 'FontSize', 36);
    
    % adjust the axes
    ax = gca;
    ax.FontSize = 36;
    ax.YLim = [0 ax.YLim(2) + 0.05];

    % print the figure
    currentFigure.PaperPositionMode = 'auto';
    currentFigure.PaperUnits = 'centimeters';
    currentFigure.PaperSize = [47 22];
    currentFigure.PaperPosition = [1 1 45 20];
    print(currentFigure, '-dpdf', 'toy-anova.pdf');