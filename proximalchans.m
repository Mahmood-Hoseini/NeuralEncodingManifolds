function [prox] = proximalchans(probetype,chan)
% Given a tetrode probe type ('64D' or '64F' and channel number)
% returns the N nearest neighbors
% Note all channels are shifted up one relative to INTAN usage 
% e.g. Chan 0 (INTAN) is Chan 1 here; grabbed from amplifier_data(1,:)
% USEAGE
%{
    [prox] = proximalchans('64D',30)
%}

if strcmp(probetype,'64D')
    arrangement = [ 1, 3, 6, 24, 29, 30;... % chan 0
                    0, 62, 24, 9, 31, 47;...% chan 1
                    5, 3, 7, 25, 27, 28;...% chan 2
                    2, 0, 25, 6, 28, 29;...% chan 3
                    26, 27, 5, 7, 28, 2;...% chan 4
                    4, 26, 27, 2, 7, 28;...% chan 5
                    25, 29, 30, 24, 0, 3;...% chan 6
                    26, 27, 5, 28, 2, 25;...% chan 7
                    23, 46, 45, 22, 50, 53;...% chan 8
                    24, 31, 47, 23, 62, 1;...% chan 9
                    21, 42, 41, 20, 52, 49;...% chan 10
                    22, 44, 43, 21, 48, 51;...% chan 11
                    19, 38, 37, 18, 56, 57;...% chan 12
                    20, 40, 39, 19, 54, 55;...% chan 13
                    17, 34, 33, 16, 60, 61;...% chan 14
                    18, 36, 35, 17, 58, 59;...% chan 15
                    14, 33, 32, 63, 60, 17;...% chan 16
                    15, 35, 34, 14, 61, 58;...% chan 17
                    12, 37, 36, 15, 59, 56;...% chan 18
                    13, 39, 38, 12, 57, 54;...% chan 19
                    10, 41, 40, 13, 55, 52;...% chan 20
                    11, 43, 42, 10, 49, 48;...% chan 21
                    8, 45, 44, 11, 51, 50;...% chan 22
                    9, 47, 46, 8, 53, 62;...% chan 23
                    6, 30, 31, 9, 1, 0;...% chan 24
                    7, 28, 29, 6, 3, 2;...% chan 25
                    4, 27, 5, 7, 28, 2;...% chan 26
                    4, 26, 5, 7, 28, 25;...% chan 27
                    27, 7, 2, 25, 29, 3;...% chan 28
                    28, 25, 3, 6, 30, 31;...% chan 29
                    29, 6, 0, 24, 31, 47;...% chan 30
                    30, 24, 1, 9, 47, 46;...% chan 31
                    33, 16, 63, 34, 14, 60;...% chan 32
                    34, 14, 60, 16, 63, 32;...% chan 33
                    35, 17, 61, 14, 33, 32;...% chan 34
                    36, 15, 58, 17, 34, 33;...% chan 35
                    37, 18, 59, 15, 35, 34;...% chan 36
                    38, 12, 56, 18, 36, 35;...% chan 37
                    39, 19, 57, 12, 37, 36;...% chan 38
                    40, 13, 54, 19, 38, 37;...% chan 39
                    41, 20, 55, 13, 39, 38;...% chan 40
                    42, 10, 52, 20, 40, 39;...% chan 41
                    43, 21, 49, 10, 41, 40;...% chan 42
                    44, 11, 48, 21, 42, 41;...% chan 43
                    45, 22, 51, 11, 43, 42;...% chan 44
                    46, 8, 50, 22, 44, 43;...% chan 45
                    47, 23, 53, 8, 45, 44;...% chan 46
                    31, 9, 62, 23, 46, 45;...% chan 47
                    51, 11, 43, 21, 49, 52;...% chan 48
                    48, 21, 42, 10, 52, 55;...% chan 49
                    53, 8, 45, 22, 51, 48;...% chan 50
                    50, 22, 44, 11, 48, 49;...% chan 51
                    49, 10, 41, 20, 55, 54;...% chan 52
                    62, 23, 46, 8, 50, 51;...% chan 53
                    55, 13, 39, 19, 57, 56;...% chan 54
                    52, 20, 40, 13, 54, 57;...% chan 55
                    57, 12, 37, 18, 59, 58;...% chan 56
                    54, 19, 38, 12, 56, 59;...% chan 57
                    59, 15, 35, 17, 61, 60;...% chan 58
                    56, 18, 36, 15, 58, 61;...% chan 59
                    61, 14, 33, 16, 63, 32;...% chan 60
                    58, 17, 34, 14, 60, 63;...% chan 61
                    1, 9, 47, 23, 53, 50;...% chan 62
                    60, 16, 32, 61, 14, 33]; % chan 63
                    
    prox = arrangement(chan,:) + 1; % add 1 to get back into normal numbering
    
else
    
    disp('Do not have probe 64F coded up yet... You need to do it!')
    
end

end