%% introduction to linear algebra and matrix analysis

% mikexcohen@gmail.com

%% vectors as lines

% 2-dimensional vector
v2 = [ 3 -2 ];

% 3-dimensional vector
v3 = [ 4 -3 2 ];


% plot them
figure(1), clf
subplot(211)
plot([0 v2(1)],[0 v2(2)],'linew',2)
axis square
axis([ -4 4 -4 4 ])
hold on
plot(get(gca,'xlim'),[0 0],'k--')
plot([0 0],get(gca,'ylim'),'k--')
xlabel('X_1 dimension')
ylabel('X_2 dimension')


subplot(212)
plot3([0 v3(1)],[0 v3(2)],[0 v3(3)],'linew',2)
axis square
axis([ -4 4 -4 4 -4 4 ])
hold on, grid on
plot3(get(gca,'xlim'),[0 0],[0 0],'k--')
plot3([0 0],get(gca,'ylim'),[0 0],'k--')
plot3([0 0],[0 0],get(gca,'zlim'),'k--')
xlabel('X_1 dimension')
ylabel('X_2 dimension')
zlabel('X_3 dimension')

% might be easier to see when rotated
rotate3d on 

%% vector-vector addition

v1 = [3 2]';
v2 = [1; 2];
v3 = v1+v2;

figure(2), clf
subplot(121)
plot([0 v1(1)],[0 v1(2)],'k','linew',2)
axis square
axis([ -6 6 -6 6 ])
hold on
plot([0 v2(1)],[0 v2(2)],'b','linew',2)
plot([0 v3(1)],[0 v3(2)],'r','linew',2)

plot(get(gca,'xlim'),[0 0],'k--')
plot([0 0],get(gca,'ylim'),'k--')
xlabel('X_1 dimension')
ylabel('X_2 dimension')

legend({'v1';'v2';'v1+v2'},'location','northwest') % also specify the location of the legend

%% vector-scalar multiplication

v2m = v2 * 1.5;

% plot them
subplot(122)
plot([0 v2(1)],[0 v2(2)],'k','linew',2)
axis square
axis([ -6 6 -6 6 ])
hold on
plot([0 v2m(1)],[0 v2m(2)],'r--','linew',2)

plot(get(gca,'xlim'),[0 0],'k--')
plot([0 0],get(gca,'ylim'),'k--')
xlabel('X_1 dimension')
ylabel('X_2 dimension')

legend({'v2';'1.5v2'})

%% inner/outer products

% inner product (aka dot product)
v1'*v2
v1(:)'*v2(:)
sum(v1.*v2)
dot(v1,v2)

% outer product
v1*v2'

%% matrix addition

aMat = [2 3; 5 7];

aMat
aMat + 5
aMat * 5

%% Three important matrices

N = 4; % for easy viewing and inspection

% square matrix
sqr = round( 10*rand(N) );

% symmetric matrix
% Use the rule that a matrix times its transpose is symmetric.
% and shrink down by 10 to make the numbers smaller for the figure ;)
sym = round( (sqr'*sqr)/10 ); 



% identity matrix
idt = eye(N);

figure(3), clf
subplot(131), imagesc(sqr), axis square, title('Square')
subplot(132), imagesc(sym), axis square, title('Symmetric')
subplot(133), imagesc(idt), axis square, title('Identity')

%% finding your way around matrices

% create a matrix by reshaping a vector
A = reshape( 1:12, 3,4)

% matrix indexing: state the element you want in the dimensions you want
A(3,2)
A(2:3,3:4)
A([1 3],[1 3 4])

% linear indexing: use one number to access elements (see figure 10.6)
A(4)
A(11)
A(3:5)
A([2 4 10 12])


% the same principle applies in higher dimensions
B = randn(4,2,6,5,3);

% please please please do yourself a favor and 
% use matrix indexing whenever possible!
B(1,2,3,2,1)

% linear indexing can get quite confusing. 
% Where in the matrix is this element??
B(582)

%% matrix multiplication

M1 = 10;
N1 =  5;
M2 = 10;
N2 =  5;

mat1 = randn(M1,N1);
mat2 = randn(M2,N2);

size( mat1*mat2' )

%%


%% transpose

v = [1 2 4 3];
v'
v''
transpose(v')

m = [ 1 3 2; 4 1 2 ];
m'
m''

%% rank

A = [1 1 3; 2 3 -2 ];

rank(A)

%% matrix inverse

A = [2 3; 1 5];
Ainv = inv(A);

figure(5), clf
subplot(131)
imagesc(A)
axis off, axis square
title('A')

subplot(132)
imagesc(Ainv)
axis off, axis square
title('A^-^1')

subplot(133)
imagesc(A*Ainv)
axis off, axis square
title('AA^-^1 = I')

%%

