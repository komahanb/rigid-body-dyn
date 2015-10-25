function [] = plotstates(pos,tt)
subplot(2,2,1)
plot(tt, pos(:,1), tt, pos(:,2), tt, pos(:,3));title('position r');
xlabel('time');

subplot(2,2,2)
plot(tt, pos(:,4), tt, pos(:,5), tt, pos(:,6));title('orientation theta');
xlabel('time');

subplot(2,2,3)
plot(tt, pos(:,7), tt, pos(:,8), tt, pos(:,9));title('linear velocity v');
xlabel('time');

subplot(2,2,4)
plot(tt, pos(:,10), tt, pos(:,11), tt, pos(:,12)); title('angular velocity omega');
xlabel('time');
end