function out = w_syllabletyper(in, num)


figure(1); clf; 
    specgram(in(num).duet, 512, in(num).Fs, [], 500);
    colormap('hot');
    ylim([500 4500]);
    caxis([-50 30]);   

    mintim = -in(num).tim(1);
    
    numsyls = length(in(num).syl);
    if num == 3
        numsyls = 13; % I HATE THAT
    end

    hold on;
    for j=1:numsyls
        plot([in(num).syl(j).tim(1)+mintim, in(num).syl(j).tim(1)+mintim], [500, 4500], 'g-');
        plot([in(num).syl(j).tim(2)+mintim, in(num).syl(j).tim(2)+mintim], [500, 4500], 'r-');
        if in(num).sylsex(j) == 1 % Male
            text(in(num).syl(j).tim(1)+mintim, 4000, num2str(j), 'FontSize', 18, 'Color', 'c');
        end
        if in(num).sylsex(j) == 2 % Female
            text(in(num).syl(j).tim(1)+mintim, 4000, num2str(j), 'FontSize', 18, 'Color', 'm');
        end
    end
    
    f = input('Female syl idents? \n');
    fidx = find(in(num).sylsex == 2);
    if length(f) == length(fidx)
        for j=1:length(fidx)
            out(fidx(j)) = f(j) + 20;
        end
    else
        fprintf('Fucked that right up, did you not?\n');
    end
    
    
    m = input('Male syl idents? \n');
    midx = find(in(num).sylsex == 1);
    if length(m) == length(midx)
        for j=1:length(midx)
            out(midx(j)) = m(j) + 10;
        end
    else
        fprintf('Fucked that right up man, did you not?\n');
    end
    