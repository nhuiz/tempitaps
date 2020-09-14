function pattern = make_rhythm_v8_2(num_measures,num_rhythms)

% function pattern = make_rhythm(num_measures,num_rhythms)
% V8 CHECKS TO MAKE SURE NON OF THE RHYTHMS ARE REPEATS OF OTHER RHYTHMS IN
% THE SET
% V8_2 doesn't make a complex rhythm, only simple. 
%
% Inputs:
%   num_measures = number of measures, self explanatory
%   num_rhythms = number of unique patterns to generate
%
% Output:
%   pattern = struct with two fields (simple, complex), each is a cell array containing as
%   many rhythms as are specified in num_patterns
%
%
% Make sure randomization is fresh.
% rng('shuffle')

% There are various things hard-coded in this function:
meter    = 4;
num_reps = 1; % Max number of repetitions allowed (max. number of measures that can happen in a row with the same structure)
benergy_max  = .33; % Proportion of beats that are allowed to get events for metric complex versions
simpbar_max = 2; % Max number of consecutive MS bars in a MC rhythm


% Make the rhythms:
for ii = 1:num_rhythms
%     comp = [];
%     while isempty(comp) || ireps == 1 || irepc == 1
        
        repcnt = 0; % check for repetitions
        for h = 1:num_measures
            m{h} = make_measure_v2(h);
            if h > 1
                if length(m{h}) == length(m{h-1})
                    if sum(m{h} == m{h-1}) == length(m{h})
                        repcnt = repcnt + 1;
                    end
                end
                while repcnt >= num_reps
                    m{h} = make_measure_v2(h);
                    if length(m{h}) == length(m{h-1})
                        if sum(m{h} == m{h-1}) == length(m{h})
                            repcnt = repcnt + 1;
                        else
                            repcnt = 0;
                        end
                    else
                        repcnt = 0;
                    end
                end
            end
        end
        
        % once non-repeating rhythms exists...
        rhy = [];
        for k = 1:length(m)
            rhy = [rhy m{k}];
        end
        
%         % check that simple rhythm isn't a repetition of another in the set
%         ireps = 0;
%         if ii > 1
%             for jj = 1:ii-1
%                 if length(rhy) == length(pattern.simple{jj,1})
%                     if sum(rhy == pattern.simple{jj,1}) == length(rhy)
%                         ireps = 1;
%                     end
%                 end
%             end
%         end
        pattern.simple{ii,1} = rhy;
        
%         % for metrical complex versions of the same rhythms
%         comp =  make_complex_rhythm(pattern.simple{ii,1},benergy_max,simpbar_max);
%         if ~isempty(comp)
%             % check that simple rhythm isn't a repetition of another in the set
%             irepc = 0;
%             if ii > 1
%                 for jj = 1:ii-1
%                     if length(comp) == length(pattern.complex{jj,1})
%                         if sum(comp == pattern.complex{jj,1}) == length(comp)
%                             irepc = 1;
%                         end
%                     end
%                 end
%             end
%             pattern.complex{ii,1} = comp;
%         end
    end
    %     keyboard;
end

% end


% subroutine for making individual measures
function y = make_measure_v2(h)

% individual measures and their probabilities:
% template{1} = [4];
% template{2} = [3 1];
% template{3} = [1 3];
% template{4} = [2 2];
% template{5} = [2 1 1];
% template{6} = [1 2 1];
% template{7} = [1 1 2];
template{8} = [1 1 1 1];

P = [1];

Ep = [];
if h < 2
    for ii = [1 2 4 5 7 8]
        Ep = [Ep; repmat(ii,P(ii)*100,1)];
    end
else
    for ii = [1 2 3 4 5 7 8]
        Ep = [Ep; repmat(ii,P(ii)*100,1)];
    end
end

[~,ix] = sort(rand(size(Ep)));
y = template{Ep(ix(1))};

end


% subroutine for permuting measures
function n = permute_measures(m)

num_measures = length(m);
x = 0;
while x < num_measures/2
    [~,p] = sort(rand(num_measures,1));
    n = m(p);
    for jj = 1:num_measures/2
        if length(n{jj}) == 2
            if (sum(n{jj}==[3 1]) == 2) || (sum(n{jj}==[2 2]) == 2)
                x = x + 1;
            elseif (sum(n{jj}==[1 3]) == 2) && (rand(1,1) <= .2)
                x = x + 1;
            end
        elseif length(n{jj}) == 3
            if (sum(n{jj}==[1 1 2]) == 3) || (sum(n{jj}==[2 1 1]) == 3)
                x = x + 1;
            elseif (sum(n{jj}==[1 2 1]) == 3) && (rand(1,1) <= .2)
                x = x + 1;
            end
        end
    end
end

end