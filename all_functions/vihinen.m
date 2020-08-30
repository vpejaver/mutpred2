function [prediction] = vihinen (sequence)

% X is added and the average residue flexibility is added for it
% (mean(FlexScales(1:20))
AminoAcids = 'ACDEFGHIKLMNPQRSTVWYXBZ';

% DS
%FlexScales = [-0.605 -0.693 -0.279 -0.160 -0.719 -0.537 -0.662 -0.682 -0.043 -0.631 -0.626 -0.381 -0.271 -0.369 -0.448 -0.423 -0.525 -0.669 -0.727 -0.721];
% V
FlexScales =[0.984 0.906 1.068 1.094 0.915 1.031 0.950 0.927 1.102 0.935 0.952 1.048 1.049 1.037 1.008 1.046 0.997 0.931 0.904 0.929 0.9906 1.068 1.094];

%DS
%window = [0.2 0.4 0.6 0.8 1 0.8 0.6 0.4 0.2];
%V
window = [0.25 0.4375 0.625 0.8125 1 0.8125 0.625 0.4275 0.25];

%load DData; % use for David Smith's data
%load DDataRaw; % use for Vihinen's data

% normalize raw data (only for Vihinen)
% for i = 1 : 290
%     mn = mean(Data{i}.Bfactors(4 : length(Data{i}.Bfactors) - 3));
%     st = std(Data{i}.Bfactors(4 : length(Data{i}.Bfactors) - 3));
%     Data{i}.Bfactors = (Data{i}.Bfactors - mn) / st;
% end


prediction = zeros(1, length(sequence));

for j = 1 : length(sequence)
    m = 1;
    vect = [];
    for k = max(1, j - 4) : min(j + 4, length(sequence))
        q = find(AminoAcids == sequence(k));
        vect(m) = FlexScales(q);
        m = m + 1;
    end
    
    if j == 1
        wind = window(5 : length(window));
        wind = wind * sum(window) / sum(wind);
    elseif j == 2
        wind = window(4 : length(window));
        wind = wind * sum(window) / sum(wind);
    elseif j == 3
        wind = window(3 : length(window));
        wind = wind * sum(window) / sum(wind);
    elseif j == 4
        wind = window(2 : length(window));
        wind = wind * sum(window) / sum(wind);
    elseif j == length(sequence) - 3
        wind = window(1 : length(window) - 1);
        wind = wind * sum(window) / sum(wind);
    elseif j == length(sequence) - 2
        wind = window(1 : length(window) - 2);
        wind = wind * sum(window) / sum(wind);
    elseif j == length(sequence) - 1
        wind = window(1 : length(window) - 3);
        wind = wind * sum(window) / sum(wind);
    elseif j == length(sequence)
        wind = window(1 : length(window) - 4);
        wind = wind * sum(window) / sum(wind);
    else
        wind = window;
    end
       
    prediction(j) = wind * vect';
end


% for j = 5 : length(sequence) - 4
%     m = 1;
%     for k = j - 4 : j + 4
%         q = find(AminoAcids == sequence(k));
%         vect(m) = FlexScales(q);
%         m = m + 1;
%     end
%     prediction(j) = window * vect';
% end


return
