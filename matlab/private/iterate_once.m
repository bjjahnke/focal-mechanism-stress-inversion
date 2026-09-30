function [choice, result, iterations, status] = iterate_once(planes, choice, friction, max_iterations)
%ITERATE_ONCE One iterative stress inversion from a starting plane choice.
%   status is 'converged', 'cycling' (choice flips between states; the state
%   with the lowest mean misfit is kept) or 'max_iterations'.
keys = {};
history = {};
for iteration = 1:max_iterations
    [normals, slips] = select_planes(planes, choice);
    result = invert_stress(normals, slips);
    keys{end + 1} = encode_choice(choice); %#ok<AGROW>
    history{end + 1} = struct('choice', choice, 'result', result); %#ok<AGROW>

    inst1 = fault_instability(planes.normal1, result, friction);
    inst2 = fault_instability(planes.normal2, result, friction);
    new_choice = 1 + (inst1 < inst2);
    if isequal(new_choice, choice)
        iterations = iteration; status = 'converged';
        return
    end
    seen = find(strcmp(keys, encode_choice(new_choice)), 1);
    if ~isempty(seen)
        cycle = history(seen:end);
        misfits = cellfun(@(h) mean(h.result.misfit_deg), cycle);
        [~, best] = min(misfits);
        choice = cycle{best}.choice; result = cycle{best}.result;
        iterations = iteration; status = 'cycling';
        return
    end
    choice = new_choice;
end
choice = history{end}.choice; result = history{end}.result;
iterations = max_iterations; status = 'max_iterations';
end
