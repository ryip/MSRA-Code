% Diagonal Down-Left DTV operator
%
% Author: Richard Yip
% Email: rbyip1@gmail.com

function op = linop_TV_DL( sz, variation, action )

error(nargchk(1,3,nargin));
if nargin < 2 || isempty(variation), variation = 'regular'; end
if nargin < 3 || isempty(action), action = 'handle'; end

CALCULATE_TV = false;
if numel(sz) > 4 
    CALCULATE_TV = true;
    X   = sz;
    sz  = size(X);
end

if iscell(sz)
    n1 = sz{1};
    n2 = sz{2};
else
    n1 = sz(1);
    n2 = sz(2);
end

% Setup the Total-Variation operators
mat = @(x) reshape(x,n1,n2);

if strcmpi(action,'matrix') || strcmpi(action,'cvx')
    disp('take care of this');
    switch lower(variation)
        case 'regular'
            e = ones(max(n1,n2),1);
            e2 = e;
            e2(n1:end) = 0;
            J = spdiags([-e2,e], 0:1,n1,n1);
            I = eye(n2);
            Dl = kron(I,J);  % vertical differences, sparse matrix
        case 'circular'
            e = ones(max(n1,n2),1);
            e2 = e;
%             e2(n1:end) = 0;
            J = spdiags([-e2,e], 0:1,n1,n1);
            J(end,1) = 1;
            I = eye(n2);
            Dl = kron(I,J);  % vertical differences, sparse matrix
    end
    if strcmpi(action,'matrix')
        op = Dl;
    else
        % "norms" is a CVX function
        op = @(X) sum( norms( [ Dl*X(:)]' ) );
    end
    return;
end

switch lower(variation)
    case 'regular'
        Dl     = @(X) vec(DLdiffs(X)); %CHANGED
        diff_l = @(X)[[zeros(1,n2-1);X(1:end-1,2:end)],zeros(n1,1)] - [[zeros(n1-1,1),X(1:end-1,2:end)];zeros(1,n2)]; %CHANGED
    case 'circular'
        % For circular version, 2 x 2 case is special.
%         error('not yet implemented');
        disp('GOD PLEASE HELP US');
        Dl     = @(X) vec( [diff(X,1,1); X(1,:) - X(end,:) ] );
        % diff_v needs to be checked
        diff_l = @(X) [X(end,:);X(1:end-1,:)] - X;

    otherwise
        error('Bad variation parameter');
end
if iscell(sz)
    Dl_transpose = @(X)      diff_l(mat(X))  ;
else
    Dl_transpose = @(X) vec( diff_l(mat(X)) );
end

TV  = @(x) ( Dl(mat(x)) );     % real to complex
TVt = @(z) ( Dl_transpose((z)) );

if CALCULATE_TV
    op = norm( TV(X), 1 );
    return;
end

if strcmpi(action,'norm')
    % to compute max eigenvalue, I use a vector
    % that is very likely to be the max eigenvector:
    %  matrix with every entry alternating -1 and 1
    even = @(n) ~( n - 2*round(n/2) );  % returns 1 if even, 0 if odd
    Y = zeros( n1 + even(n1), n2 + even(n2) );
    nn = numel(Y);
    Y(:) = (-1).^(1:nn);
    Y = Y(1:n1,1:n2);
    
    % Nearly equivalent to:
    % norm(full( [real(tv); imag(tv)] ) ) 
    op = norm( TV(Y) )/norm(Y(:));
    % where tv is the matrix form
else
    if iscell(sz)
        szW = { [n1,n2], [n1*n2,1] };
    else
        szW = sz(1)*sz(2); % n1 * n2
        szW = [szW,szW];
    end
    op = @(x,mode)linop_tv_r2c(szW,TV,TVt,x,mode);
end

function y = linop_tv_r2c( sz, TV, TVt, x, mode )
switch mode,
    case 0, y = sz;
    case 1, y = TV( realcheck( x ) );
    case 2, y = realcheck( TVt( x ) );
end

function y = realcheck( y )
if ~isreal( y ), 
    error( 'Unexpected complex value in linear operation.' );
end

% TFOCS v1.1 by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2012 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.
