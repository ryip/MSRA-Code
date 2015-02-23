% Improvement on original DTV solver by implementing vertical, horizontal,
% and diagonal encouragement norms.
% Modifications done by Richard Yip

% DTV solver implementing transform-invariant directional total-variation together
% with l2-norm penalty and total-variation penalties to perform image
% super-resolution For more information see "Super-resolution via 
% Transform-invariant Group-sparse Regularization" by C. Fernandez-Granda.
%
% Author: Carlos Fernandez-Granda 
% Email: cfgranda@stanford.edu

function varargout = solver_transf_DTV_L2_TV_diag( A, b, trans_op_x,trans_op_y,trans_op_DR,trans_op_DL,aux_op,lambda,TV_par, mu, x0, z0,L,n1,n2,m1,m2, opts, varargin )

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
W1 = trans_op_x;
W2 = trans_op_y;
W3   = linop_compose(linop_TV( [n1,n2] ),aux_op);

%DIAGONAL COMPONENTS:
W4 = trans_op_DR;
W5 = trans_op_DL;

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
normW12 = linop_normest_mod( W1,[L,1],'r2r',1e-5,1e3).^2;
normW22 = linop_normest_mod( W2,[L,1],'r2r',1e-5,1e3).^2;
normW42 = linop_normest_mod( W4,[L,1],'r2r',1e-5,1e3).^2;
normW52 = linop_normest_mod( W5,[L,1],'r2r',1e-5,1e3).^2;

pow = 2;

proxScale1 = ( normW12 / normA2 ).^(pow);
proxScale2 = ( normW22 / normA2 ).^(pow);
proxScale4 = ( normW42 / normA2 ).^(pow);
proxScale5 = ( normW52 / normA2 ).^(pow);
normW      = linop_TV( [n1,n2], [], 'norm' );
proxScale  = sqrt( normW / normA2 );
W3 = linop_compose( W3, 1 / proxScale );

power = 0.9;

prox_2       = { smooth_quad, ...
               proj_linfl2_T_mod_diag( [proxScale1 *lambda ,m1,m2,n1,n2,power] ),...
               proj_linfl2_mod_diag( [proxScale2 *lambda ,m1,m2,n1,n2,power] ) ,...
               proj_linfl2_DR_mod_diag( [proxScale4 *lambda ,m1,m2,n1,n2,power] ) ,...
               proj_linfl2_DL_mod_diag( [proxScale5 *lambda ,m1,m2,n1,n2,power] ) ,...
               proj_linf(TV_par*proxScale )};
W1         = linop_compose( W1, 1 / proxScale1 );
W2         = linop_compose( W2, 1 / proxScale2 );
W4         = linop_compose( W4, 1 / proxScale4 );
W5         = linop_compose( W5, 1 / proxScale5 );

%Change proxscale so that it takes inputs to raise powers


[varargout{1:max(nargout,1)}] = ...
    tfocs_SCD( prox, { A, -b; W1, 0; W2, 0  ; W4,0 ;W5,0;W3,0 }, prox_2, mu, x0, z0, opts, varargin{:} );


