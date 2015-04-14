% Modified projection onto the l_infinity l_2 ball that incorporates a diagonal component for our DTV
% Original code done by C. Fernandez-Granda
% Modifications done by Richard Yip

% Projection onto the l_infinity l_2 ball
%
% Author: Carlos Fernandez-Granda
% Email: cfgranda@stanford.edu

function op = proj_linfl2_DL_mod_diag( param )

if nargin == 0,
    q = 1;
    param=[1,64,64];
else
    q=param(1);
    if ~isnumeric( q ) || ~isreal( q ) || numel( q ) ~= 1 || q <= 0,
        error( 'Argument must be positive.' );
    end
end
op = @(varargin)proj_linfl2_q( param,varargin{:} );

function [ v, x ] = proj_linfl2_q( param, x, t )
q=param(1);
n1=param(2);
n2=param(3);
m1=param(4);
m2=param(5);
pow=param(6);
if m1*m2 == length(x)
    n1=m1;
    n2=m2;
end
v = 0;
x=reshape(x,n1,n2);
aux = sqrt(DLsums(x));
switch nargin,
    case 2,
        if nargout == 2,
            error( 'This function is not differentiable.' );
        elseif norm( aux(:), Inf ) > q,
            v = Inf;
        end
    case 3,
        %CHANGED TO FIX OVERCOMPENSATION FOR DIAGONAL
        %  x = x ./ (2*max(1,DLrepmat(aux)/q));
        x = x ./ max(1,DLrepmat(aux)/q);
        %pow = 0.5; %THIS IS THE POWER TO WHICH I RAISE X
        x = x./(max(abs(x),0.001).^(1-pow));
        x=x(:);
    otherwise,
        error( 'Not enough arguments.' );
end
