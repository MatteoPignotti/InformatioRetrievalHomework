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
x = [0.5208, 0.0772, 0.4262, 0.2853, 0.2444, 0.0, 0.3926, 0.3137, 0.0714, 0.3576, 0.4444, 0.2051, 0.125, 0.5143, 0.6857, 0.4141, 0.2011, 0.4098, 0.3077, 0.2292, 0.1176, 0.3265, 0.4545, 0.3775, 0.325, 0.049, 0.3077, 0.1122, 0.125, 0.2857, 0.0714, 0.5909, 0.089, 0.2549, 0.2442, 0.0526, 0.1647, 0.3529, 0.0412, 0.2313, 0.3202, 0.4476, 0.0986, 0.1176, 0.2207, 0.322, 0.3704, 0.2069, 0.1667, 0.456];
y = [0.3958, 0.0, 0.418, 0.0693, 0.1556, 0.0, 0.0, 0.3137, 0.1071, 0.2384, 0.4444, 0.2051, 0.0625, 0.4857, 0.6571, 0.2929, 0.1217, 0.377, 0.1538, 0.003, 0.0588, 0.2857, 0.0, 0.3137, 0.175, 0.0, 0.2821, 0.0, 0.0625, 0.1429, 0.0357, 0.0455, 0.0959, 0.0, 0.2209, 0.0526, 0.1529, 0.3333, 0.0361, 0.209, 0.0, 0.0286, 0.0845, 0.0, 0.0516, 0.339, 0.3704, 0.0552, 0.1569, 0.416];
z = [0.5208, 0.0772, 0.4344, 0.2715, 0.2444, 0.0, 0.4, 0.3333, 0.1071, 0.3444, 0.4444, 0.2051, 0.125, 0.5143, 0.6857, 0.3939, 0.2063, 0.4098, 0.3077, 0.2083, 0.1176, 0.3265, 0.4545, 0.3775, 0.325, 0.0686, 0.3077, 0.1122, 0.125, 0.2857, 0.1071, 0.5909, 0.0822, 0.2353, 0.2674, 0.0526, 0.1882, 0.3333, 0.0412, 0.2388, 0.3427, 0.4381, 0.0986, 0.1176, 0.2254, 0.339, 0.3704, 0.2069, 0.1569, 0.456];
w = [0.4583, 0.0691, 0.4508, 0.2216, 0.1333, 0.0, 0.4074, 0.3137, 0.1429, 0.5033, 0.4444, 0.2308, 0.0, 0.6, 0.7143, 0.2323, 0.1693, 0.459, 0.2308, 0.1756, 0.1765, 0.2857, 0.4545, 0.3529, 0.325, 0.0392, 0.2821, 0.102, 0.3125, 0.1429, 0.0357, 0.4091, 0.089, 0.2353, 0.186, 0.1053, 0.1529, 0.2549, 0.0155, 0.1642, 0.4438, 0.2381, 0.1127, 0.1176, 0.2066, 0.3559, 0.1481, 0.1931, 0.098, 0.432];

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