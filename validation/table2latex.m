function table2latex(T,fmt)
% TABLE2LATEX  Format table T as a LaTeX tabular (booktabs style).
%   fmt (optional) is a format string for numeric values, e.g. '%.4f'.

if nargin < 2 || isempty(fmt)
    fmt = '%.4f';
end

vars = T.Properties.VariableNames;
rowNames = T.Properties.RowNames;

if isempty(rowNames)
    rowNames = cell(height(T), 1);
    for r = 1:height(T)
        rowNames{r} = sprintf('row_%d', r);
    end
end

nVars = numel(vars);
nRows = height(T);

header = cell(1, nVars + 1);
header{1} = '';
for c = 1:nVars
    header{c + 1} = latexEscape(vars{c});
end

rows = cell(nRows, nVars + 1);
for r = 1:nRows
    rows{r, 1} = latexEscape(rowNames{r});
    for c = 1:nVars
        rows{r, c + 1} = formatCell(T{r, c}, fmt);
    end
end

widths = zeros(1, nVars + 1);
for c = 1:(nVars + 1)
    widths(c) = strlength(header{c});
    for r = 1:nRows
        widths(c) = max(widths(c), strlength(rows{r, c}));
    end
end

colSpec = ['l', repmat('r', 1, nVars)];
fprintf('\\begin{tabular}{%s}\n', colSpec);
fprintf('\\toprule\n');
fprintf('%s \\\\\n', joinPaddedRow(header, widths));
fprintf('\\midrule\n');
for r = 1:nRows
    fprintf('%s \\\\\n', joinPaddedRow(rows(r, :), widths));
end
fprintf('\\bottomrule\n');
fprintf('\\end{tabular}\n');
end

function s = formatCell(v, fmt)
if isnumeric(v) && isscalar(v)
    if isnan(v)
        s = '';
    else
        s = sprintf(fmt, v);
    end
elseif isstring(v) || ischar(v)
    s = char(v);
elseif iscategorical(v)
    s = char(v);
else
    s = char(string(v));
end
s = latexEscape(s);
end

function s = latexEscape(s)
s = strrep(s, '&', '\&');
s = strrep(s, '%', '\%');
s = strrep(s, '_', '\_');
s = strrep(s, '#', '\#');
s = strrep(s, '{', '\{');
s = strrep(s, '}', '\}');
end

function line = joinPaddedRow(cells, widths)
n = numel(cells);
parts = cell(1, n);
for c = 1:n
    pad = widths(c) - strlength(cells{c});
    parts{c} = [cells{c}, repmat(' ', 1, pad)];
end
line = strjoin(parts, ' & ');
end
