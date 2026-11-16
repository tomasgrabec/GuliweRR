classdef ParallelPolicy
%%% ParallelPolicy — shared parallel-execution policy for GuliweRR sweep classes
%%%
%%% Two decisions, made in one place for every sweep class: whether a loop of
%%% independent iterations runs in parallel ("auto"/"always"/"never"), and
%%% acquisition of a MATLAB THREAD pool. Threads, not processes, because a
%%% thread pool shares memory instead of copying broadcast data (elastic
%%% tensors, CalcParams, rawByDepthCell) into every worker.
%%%
%%% Stateless static class.
%%%
%%% GuliweRR v6.0

    properties (Constant)
        AUTO_THRESHOLD = 20  % "auto" goes parallel once nItems >= this; a single
                              % benchmarked threshold, not scaled by problem size
    end

    methods (Static)
        function tf = shouldUseParallel(useParallelSetting, nItems)
        %%% Decide whether a loop of nItems independent iterations should run in parallel
        %%% >> tf = ParallelPolicy.shouldUseParallel(useParallelSetting, nItems)
        %%%
        %%% INPUT:  1. useParallelSetting .. "auto" | "always" | "never"
        %%%         2. nItems             .. number of independent loop iterations
        %%% OUTPUT: tf .. true/false

            useParallelSetting = string(useParallelSetting);
            validateattributes(nItems, {'numeric'}, {'scalar', 'real', 'nonnegative'}, ...
                'ParallelPolicy.shouldUseParallel', 'nItems');

            switch useParallelSetting
                case "never"
                    tf = false;
                case "always"
                    tf = true;
                case "auto"
                    tf = nItems >= ParallelPolicy.AUTO_THRESHOLD;
                otherwise
                    error('ParallelPolicy:shouldUseParallel:badSetting', ...
                        "useParallelSetting must be 'auto', 'always', or 'never'; got '%s'.", ...
                        useParallelSetting);
            end
        end

        function tf = isOnWorker()
        %%% True when running inside a parallel worker (thread or process pool)
        %%% >> tf = ParallelPolicy.isOnWorker()
        %%%
        %%% OUTPUT: tf .. true on a pool worker, false on the client (or without the
        %%%               Parallel Computing Toolbox)
        %%%
        %%% A loop that goes parallel on its own (e.g. GreenScore.reference) must stay
        %%% serial on a worker: there is no pool there (gcp is empty), and
        %%% ensureThreadPool would try to start one. getCurrentTask does not work in a
        %%% thread pool (it errors), so thread workers are recognised by MATLAB's
        %%% internal parallel.internal.pool.isPoolThreadWorker (present in R2025b) and
        %%% process workers by getCurrentTask. Any failure means "client".
            tf = false;
            try
                tf = parallel.internal.pool.isPoolThreadWorker();
            catch
            end
            if ~tf
                try
                    tf = ~isempty(getCurrentTask());
                catch
                end
            end
        end

        function ok = ensureThreadPool()
        %%% Make a MATLAB thread pool the active parallel pool, starting one if needed
        %%% and swapping out any active process pool
        %%% >> ok = ParallelPolicy.ensureThreadPool()
        %%%
        %%% OUTPUT: ok .. true if a thread pool is ready; false (with a warning, not an
        %%%               error) if parallel execution is unavailable here — callers
        %%%               must treat false as "run serial"
        %%%
        %%% An active process pool is deleted to make room (MATLAB allows one pool),
        %%% with the warning ParallelPolicy:ensureThreadPool:replacedPool.

            ok = true;
            try
                pool = gcp('nocreate');
                if isempty(pool) || ~isa(pool, 'parallel.ThreadPool')
                    if ~isempty(pool)
                        warning('ParallelPolicy:ensureThreadPool:replacedPool', ...
                            ['Deleting the active %s to start the thread pool GuliweRR''s ' ...
                             'parallel sweeps use.'], class(pool));
                        delete(pool);   % swap a process pool for a thread pool
                    end
                    parpool('Threads');
                end
            catch ME
                warning('ParallelPolicy:ensureThreadPool:unavailable', ...
                    ['Could not start a MATLAB thread pool (%s) -- falling back to ' ...
                     'serial execution. Requires Parallel Computing Toolbox.'], ME.message);
                ok = false;
            end
        end
    end
end
