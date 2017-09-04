% Copyright (c) 2017.
% All rights reserved. Please read the 'license.txt' for license terms.
% 
% Developers: Zhen Zhang, Dr Fumio Motegi, Dr Pakorn Kanchanawong
% Contact: 
% Dr Pakorn Kanchanawong (biekp@nus.edu.sg)
% Dr Fumio Motegi (fmotegi@tll.org.sg)

function thresh = SubregionOtsu(counts)
total = sum(counts); % '''total''' is the number of pixels in the given image.
S1 = 0;
W1 = 0;
M = 0.00;
S2 = dot( (0:255), counts);
i = 1;
while i<=256
    W1 = W1 + counts(i);
    if W1 == 0
        continue;
    end
    W2 = total - W1;
    if W2 == 0
        break;
    end
    S1 = S1 +  (i-1) * counts(i);
    bwn = W1 * W2 * ((S1/W1) - ((S2-S1)/W2)) * ((S1/W1) - ((S2-S1)/W2));
    if ( bwn >= M )
        thresh = i;
        M = bwn;
    end
    i = i + 1;
end
end