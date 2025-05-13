function tracedArms = manualArmTracing(img)
    figure, imshow(img, []);
    title('Click along each spiral arm. Press Enter when done with each arm. Press Enter twice to finish.');
    hold on;
    tracedArms = {};
    armCount = 0;
    while true
        [x, y] = getpts;  % User clicks along one arm, Enter to finish
        if isempty(x)
            break;  % Double Enter to finish all arms
        end
        armCount = armCount + 1;
        tracedArms{armCount} = [x, y];
        plot(x, y, 'LineWidth', 2, 'DisplayName', ['Arm ' num2str(armCount)]);
        scatter(x, y, 25, 'filled');
    end
    legend show;
    hold off;
    disp(['Total arms traced: ' num2str(armCount)]);
end
