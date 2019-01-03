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
x = [1., 0.5, 0.3, 0.3, 0.2, 0., 1., 0.5, 0., 0.4, 0.5, 0.2, 0.1, 0.8, 1., 0.8, 0.3, 1., 0.2, 0., 0., 0.5, 0.6, 0.4, 0.7, 0.2, 0.5, 0., 0.1, 0.2, 0.1, 0.8, 0.2, 0.4, 0.9, 0.1, 0.2, 0., 0.5, 0.7, 0.2, 0.9, 0.5, 0., 0.1, 1., 0.9, 0., 0.4, 0.7];
y = [1., 0.4, 0.3, 0.3, 0.3, 0., 0.9, 0.5, 0., 0.5, 0.5, 0.2, 0.1, 0.8, 1., 0.6, 0.4, 0.9, 0.2, 0., 0., 0.5, 0.6, 0.5, 0.7, 0.1, 0.5, 0., 0.1, 0.2, 0.1, 0.8, 0.3, 0.4, 0.9, 0.1, 0.2, 0., 0.5, 0.7, 0.9, 0.9, 0.5, 0., 0., 1., 0.9, 0., 0.4, 0.8];
z = [1., 0.5, 0.3, 0.3, 0.2, 0., 1., 0.5, 0.1, 0.4, 0.5, 0.2, 0.1, 0.8, 1., 0.8, 0.3, 1., 0.2, 0., 0., 0.4, 0.6, 0.4, 0.7, 0.2, 0.5, 0., 0.1, 0.2, 0.1, 0.8, 0.2, 0.4, 0.9, 0.1, 0.2, 0., 0.5, 0.7, 0.2, 0.9, 0.5, 0.1, 0.1, 1., 0.9, 0., 0.4, 0.7];
w = [1., 0.4, 0.7, 0.4, 0.3, 0., 0.8, 0., 0.1, 0.9, 0.4, 0.2, 0.1, 0.8, 1., 0.6, 0.4, 0.8, 0.2, 0.1, 0., 0.4, 0.8, 0.4, 0.7, 0.2, 0.3, 0., 0.5, 0.2, 0.1, 0.6, 0.1, 0.3, 0.8, 0., 0.3, 0.1, 0.3, 0.8, 0.9, 0.3, 0.6, 0., 0., 1., 0.4, 0., 0.2, 0.8];

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