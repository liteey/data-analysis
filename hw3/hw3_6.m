
prin_comp = U' * fourthMat;
plot(1:trimLen, prin_comp(3:5,:),'Linewidth',1)
xlabel('Time (frames)'); ylabel('Displacement (pixels)')
legend('3rd Principal Component', '4th Principal Component', ...
       '5th Principal Component', 'location', 'northeast')
title('Case 4: Three to Five Principal Component Projections')