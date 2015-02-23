% Improvement on original DTV solver by implementing vertical, horizontal,
% and diagonal encouragement norms.
% Modifications done by Richard Yip


% DTV solver implementing transform-invariant directional total-variation together
% with l2-norm penalty to perform image super-resolution For more information
% see "Super-resolution via Transform-invariant Group-sparse Regularization"
% by C. Fernandez-Granda.
%
% Author: Carlos Fernandez-Granda 
% Email: cfgranda@stanford.edu

function varargout = solver_transf_DTV_L2_diag( A, b, trans_op_x,trans_op_y, trans_op_DR, trans_op_DL, lambda, mu, x0, z0,L,n1,n2,m1,m2, opts, varargin )

nonneg = false;
if isfield(opts,'nonneg')
    nonneg  = opts.nonneg;
    opts = rmfield(opts,'nonneg');
end
if isfield(opts,'nonNeg')
    nonneg  = opts.nonNeg;
    opts = rmfield(opts,'nonNeg');
end

if nonneg       
    % -- case: x >= 0 constraints
    prox    = proj_Rplus;
else
    % -- case: no x >= 0 constraint
    prox    = [];
end

%W1= linop_compose(trans_op,linop_TV_x( [m1,m2] ));
%W2= linop_compose(trans_op,linop_TV_y( [m1,m2] ));

%This looks like what I need
W1 = trans_op_x;
W2 = trans_op_y;
W3 = trans_op_DR;
W4 = trans_op_DL;

% Need to estimate the norms of A*A' and W*W' in order to be most efficient
if isfield( opts, 'noscale' ) && opts.noscale,
    normA2 = 1; normW12 = 1; normW22 = 1;
else
    normA2 = []; normW12 = []; normW22 = [];
    if isfield( opts, 'normA2'  )
        normA2 = opts.normA2;
        opts = rmfield( opts, 'normA2' );
    end
end
if isempty( normA2 ),
    normA2 = linop_normest( A ).^2;
end

%This looks like the thing in the numerator... 
%So I need to "make up" a normW32 (down right) and normW42 (down left)

normW12 = linop_normest_mod( W1,[L,1],'r2r',1e-5,1e3).^2; %Does this mean linop_normest_mod returns a list of values?
normW22 = linop_normest_mod( W2,[L,1],'r2r',1e-5,1e3).^2;
normW32 = linop_normest_mod( W3,[L,1],'r2r',1e-5,1e3).^2;
normW42 = linop_normest_mod( W4,[L,1],'r2r',1e-5,1e3).^2;

%This right here looks promising
%proxScale1 corresponds with the x-differential

proxScale1 = sqrt( normW12 / normA2 );
proxScale2 = sqrt( normW22 / normA2 );
proxScale3 = sqrt( normW32 / normA2 );
proxScale4 = sqrt( normW42 / normA2 ); 

%Must alter prox_2 so that it takes additional inputs (namely stuff with
%proxScale3-4

%finish the sentence: proxScale2 takes the form of __________

prox_2       = { smooth_quad, ...
               proj_linfl2_T_mod( [proxScale1 *lambda ,m1,m2,n1,n2] ),...
               proj_linfl2_mod( [proxScale2 *lambda ,m1,m2,n1,n2] ) };
W1         = linop_compose( W1, 1 / proxScale1 );
W2         = linop_compose( W2, 1 / proxScale2 );
W3         = linop_compose( W3, 1 / proxScale3 );
W4         = linop_compose( W4, 1 / proxScale4 );


[varargout{1:max(nargout,1)}] = ...
    tfocs_SCD( prox, { A, -b; W1, 0; W2, 0  }, prox_2, mu, x0, z0, opts, varargin{:} );


