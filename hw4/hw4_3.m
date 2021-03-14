%%
hold on
plot3(  ones(3,:),   ones(2,:),   ones(5,:),'ko')
plot3(    zs(3,:),     zs(2,:),     zs(5,:),'ro')
plot3(  twos(3,:),   twos(2,:),   twos(5,:),'go')
plot3(threes(3,:), threes(2,:), threes(5,:),'yo')
plot3( fours(3,:),  fours(2,:),  fours(5,:),'co')
plot3( fives(3,:),  fives(2,:),  fives(5,:),'bo')
plot3(  sixs(3,:),   sixs(2,:),   sixs(5,:),'o', 'Color', [.61 .51 .74])
plot3(sevens(3,:), sevens(2,:), sevens(5,:),'wo')
plot3(eights(3,:), eights(2,:), eights(5,:),'o', 'Color', [.31 .41 .92])
plot3( nines(3,:),  nines(2,:),  nines(5,:),'o', 'Color', [.87 .65 .11])
% proj = S*V(:,1)';


% seven and nine stick together
% pro
%% Build LDA for two digit classification
numSev = size(sevenInd,2);
numNin = size(nineInd,2);
meanSev = mean(sevens,2);
meanNin = mean(nines,2);

wn = 0; % within variances
for k = 1:numSev
   wn = wn + (sevens(:,k)-meanSev)*(sevens(:,k)-meanSev)'; 
end

for k = 1:numNin
    wn = wn + (nines(:,k) - meanNin)*(nines(:,k)-meanNin)';
end

bw = (meanSev-meanNin)*(meanSev-meanNin)';

% finding proj using LDA
[V2, D] = eig(bw, wn); 
[lambda, ind] = max(abs(diag(D)));
w = V2(:,ind);
w = w/norm(w,2);

projSev = w'*sevens;
projNin = w'*nines;

if mean(projSev) > mean(projNin)
    w = -w;
    projSev = -projSev;
    projNin = -projNin;
end

% nines below 0
% sevens above 0
figure
plot(projSev,sevens(length(projSev)),'ko')
hold on
plot(projNin,nines(length(projNin)),'bo')

% there is pretty big overlap between the values
