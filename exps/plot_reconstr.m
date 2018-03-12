function plot_reconstr(Y,Y_rec,obs_idx,miss_idx)
%PLOT_RECONSTR Visualizes the reconstructions along each dimension.
%
%   INPUT
%       Y - original training data
%       Y_rec - reconstructed data
%       obs_idx - indices of the observed data
%       miss_idx - indices of the missing part

%%
figure;
for y = 1:25
    subplot(5,5,y);
    hold on
    plot(obs_idx,Y(obs_idx,y),'bo','markersize',3)
    plot(miss_idx,Y(miss_idx,y),'o','markersize',3,'color',[0,0.73,0.22])
    plot(Y_rec(:,y),'r-','linewidth',2)
    if y==1, legend('obs','miss','rec'), end
end
figure;
for y = 1:25
    subplot(5,5,y);
    hold on
    plot(obs_idx,Y(obs_idx,y+25),'bo','markersize',3)
    plot(miss_idx,Y(miss_idx,y+25),'o','markersize',3,'color',[0,0.73,0.22])
    plot(Y_rec(:,y+25),'r-','linewidth',2)
    if y==1, legend('obs','miss','rec'), end
end

end

