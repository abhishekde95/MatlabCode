% This function maps events (specific codes dropped by REX) to behavior that
% sets up the next stimulus and updates the QUEST functions based on observer
% performance.

function p = parse_trial_codes(p, codes)
global gl

code_func_map = {
    codes.CORRECTCD @fCorrect
    codes.ERRCD @fWrong
    codes.EOTCD @fEot
    };

for codenum = 1:size(code_func_map,1)
    if any(p.events == code_func_map{codenum,1})
        feval(code_func_map{codenum,2});
        codeidx = find(p.events == code_func_map{codenum,1}, 1, 'last');
        p.lastprocessed_t = p.times(codeidx);
    end
end

    function fCorrect
        gl.conds_remaining(gl.cond_idx) = gl.conds_remaining(gl.cond_idx)-1;
        update_QUEST(1);
    end

    function fWrong
        gl.conds_remaining(gl.cond_idx) = gl.conds_remaining(gl.cond_idx)-1;
        update_QUEST(0);
    end

    function fEot
        pick_next_cond();
        message_REX('CONDIDXUPDATED'); % this handshake is crucial
        % without it, we may skip over stimuli and show duplicates
    end

    function pick_next_cond()
        candidates = find(gl.conds_remaining > 0);
        if isempty(candidates)
            gl.blocks_remaining = gl.blocks_remaining-1;
            if gl.blocks_remaining > 0
                gl.conds_remaining(:) = gl.trials_per_block/2;
                gl.cond_idx = randi(length(gl.conds_remaining));
            else
                return
            end
        else
            gl.cond_idx = candidates(randi(length(candidates)));
        end
    end
end
