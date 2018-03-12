% ggplot2 color scheme, i.e. equidistant colors in HSV scale
%
function C = ggplotcolors(n)
        switch n
                case 1
                        C = [0 0 0];
                case 2
                        C = hex2rgb({'#F8766D' '#00BFC4'});
                case 3
                        C = hex2rgb({'#F8766D' '#00BA38' '#619CFF'});
                case 4
                        C = hex2rgb({'#F8766D' '#7CAE00' '#00BFC4' '#C77CFF'});
                case 5
                        C = hex2rgb({'#F8766D' '#A3A500' '#00BF7D' '#00B0F6' '#E76BF3'});
                case 6
                        C = hex2rgb({'#F8766D' '#B79F00' '#00BA38' '#00BFC4' '#619CFF' '#F564E3'});
                case 7
                        C = hex2rgb({'#F8766D' '#C49A00' '#53B400' '#00C094' '#00B6EB' '#A58AFF' '#FB61D7'});
                case 8
                        C = hex2rgb({'#F8766D' '#CD9600' '#7CAE00' '#00BE67' '#00BFC4' '#00A9FF' '#C77CFF' '#FF61CC'});
                case 9
                        C = hex2rgb({'#F8766D' '#D39200' '#93AA00' '#00BA38' '#00C19F' '#00B9E3' '#619CFF' '#DB72FB' '#FF61C3'});
                case 10
                        C = hex2rgb({'#F8766D' '#D89000' '#A3A500' '#39B600' '#00BF7D' '#00BFC4' '#00B0F6' '#9590FF' '#E76BF3' '#FF62BC'});
                case 11
                        C = hex2rgb({'#F8766D' '#DB8E00' '#AEA200' '#64B200' '#00BD5C' '#00C1A7' '#00BADE' '#00A6FF' '#B385FF' '#EF67EB' '#FF63B6'});
        end
end

