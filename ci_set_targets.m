function target = ci_set_targets(Mat, thresh)
%SETTARGETS creates cyclic invariant matrices to become penalty terms for
%the next run of ci_cp_dgn, it uses the previous iteration.

    target = Mat;
    % Cap large values
    target = min(target,1);
    target = max(target,-1);

    % Set small entries to zero
    mask = abs(target) > thresh;
    target = target .* mask;

    target = roundWithThreshold(target, thresh);
end