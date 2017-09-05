%% calculate average score and choose reliable score into the ensemble

if( hogScore > hog_low_thres * hogAver )
    hogAver  = (hogAver*hogScorenum + hogScore)/(hogScorenum + 1);
    hogScorenum = hogScorenum + 1;
end

if ( colorScore > color_low_thres * colorAver )
    colorAver  = (colorAver*colorScorenum + colorScore)/(colorScorenum + 1);
    colorScorenum = colorScorenum + 1;
end

%% adaptive update
if (colorScore > 0.7 * colorAver)&&(hogScore > 0.7 * hogAver)   % looser thresholds are used to contain more useful samples  
   p.learning_rate_pwp = 0.02;
   p.learning_rate_cf = 0.01;
else
   p.learning_rate_pwp = 0;                                   % we want pure color model, just discard occluded samples
   p.learning_rate_cf = 0.8*(hogScore/hogAver)^3*0.01;  % penalize the sample with low hog_score, which still remains a little useful info
end

