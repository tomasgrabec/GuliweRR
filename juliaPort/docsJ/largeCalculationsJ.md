# GuliwerrJ — large calculations

How to use every core of a machine when a calculation runs for minutes to days: many
orientations, materials or parameter points, each one an ordinary GuliwerrJ sweep. For a
single sweep nothing here is needed — `-t auto` (`quickStartJ.md`) already spreads its angles
over the threads.

## 1. Choose the kind of parallelism

| Workload | Use | Why |
|---|---|---|
| One sweep (tens to hundreds of angles) | `julia -t auto`, the sweep's own threads (`useParallel = :auto`) | No set-up cost; all angles share one copy of the structure and the Ritz precompute |
| Many independent sweeps, a few minutes in total | Threads over the sweeps, each sweep `useParallel = :never` (section 5) | Starting worker processes would cost more than it saves |
| Many independent sweeps, longer than ~1 minute | **Worker processes** (`Distributed`, section 3), one thread each | Scale much better than threads (section 2) |

"Independent" means a job needs nothing from another job: one crystal orientation, one
material, one point of a parameter grid.

## 2. Why processes beat threads for batches

All threads of one Julia process share one memory heap. Every garbage collection — Julia's
automatic freeing of memory that is no longer used — stops all threads together, and GuliwerrJ's
solvers allocate a lot of short-lived memory. On macOS, Apple's Accelerate library (which
GuliwerrJ uses for its linear algebra, `matlabVsJuliaJ.md` §3.4) also slows down when several threads
of one process make many small LAPACK calls. Separate processes have separate heaps and do not
share these costs.

Measured on an Apple M4 Pro (8 performance + 4 efficiency cores), Ni half-space,
`CalcParams(relDepth = 20, legOrder = 50)`, 180 directions per orientation:

| Workload | 8 threads | 8 processes |
|---|---|---|
| `angSweepHS` `calcMode = :full` (Ritz + Green + classification), 32 orientations | 4.1× faster than one thread | 5.9× |
| 200-point Green grid + peak refinement per direction, 5005 orientations | 385 s | 203 s (incl. worker start-up) |

Memory bandwidth was not the limit in either case. No start-up option changed the thread
result by more than ~5 % (garbage-collector threads, concurrent sweeping, heap-size hint,
OpenBLAS instead of Accelerate, Accelerate forced single-threaded).

The price of processes:

- **RAM.** Each worker is a complete Julia process with GuliwerrJ loaded: 0.8–1.5 GB each in
  the cases above (8 workers 6–12 GB, against ~2 GB for 8 threads).
- **Start-up.** Starting the workers and compiling the code on each takes ~10 s, once per run.
- **Results are copied** back to the main process. Return what you need (velocities, scores,
  a few numbers) — never whole sweep objects, which can be hundreds of MB each.

## 3. Recipe: a batch on worker processes

Save as e.g. `myBatch.jl` and run with `julia --project=juliaPort/guliwerrJ myBatch.jl outDir`
(no `-t` needed: each worker is started with one thread).

```julia
using Distributed, Serialization
addprocs(8; exeflags = ["--project=$(Base.active_project())", "-t 1"])   # number of workers: section 4

@everywhere begin                        # defined on the main process and on every worker
    using GuliwerrJ, Serialization
    function runJob(job)                 # ONE independent calculation
        P, q = job
        mat = Material(elConsts = ElastConsts([247.0, 153.0, 122.0], 2), rho = 8.912,
                       dirMat = Euler2DirMat(0.0, P, q))
        str = Structure(materIn = mat, calcParamsIn = CalcParams(relDepth = 20, legOrder = 50),
                        name = "Ni", nMaxModes = 500)
        sweep = calc(angSweepHS(structure = str, angles = collect(0.0:1.0:179.0), calcMode = :full,
                                retentionLevel = :none, useParallel = :never))
        return sweep.sawVelMap           # small result: SAW velocities per angle
    end
end

outDir = ARGS[1]; mkpath(outDir)
jobs = [(P, q) for P in 0.0:5.0:50.0 for q in 0.0:5.0:90.0]
results = pmap(eachindex(jobs)) do i
    file = joinpath(outDir, "job_$i.jls")
    isfile(file) && return deserialize(file)         # done in an earlier, interrupted run
    r = with_blas_threads(1) do
        with_task_blas_single(() -> runJob(jobs[i]))
    end
    serialize(file, r)                               # each result saved as soon as it exists
    return r
end
rmprocs(workers())
serialize(joinpath(outDir, "results.jls"), (jobs = jobs, results = results))
```

What each part does:

- `addprocs` starts the worker processes with the same environment (`--project`) as the main
  one. `@everywhere` runs the block on all of them, so each knows `runJob`.
- `pmap` hands the jobs out one at a time to whichever worker is free, so fast and slow jobs
  balance themselves.
- `with_blas_threads(1)` and `with_task_blas_single` keep the linear algebra on one thread per
  job (`pitfallsJ.md` jNUM1); several BLAS threads per worker would only compete with the other
  workers.
- `useParallel = :never` inside the job: the worker has one thread anyway; it states the intent.
- The per-job files make the batch **restartable**: run the same command again after an
  interruption and the finished jobs are read back instead of recomputed. Delete `outDir` to
  start afresh. Serialised files are for the same Julia version only; store final results in
  your own format if you need them long term.

A job that throws an error stops `pmap`. For long unattended runs, catch errors inside the
`pmap` body (`try ... catch e; return sprint(showerror, e); end`) and check the results
afterwards.

## 4. How many workers

- **Start from the number of physical cores.** Logical "threads" from SMT (simultaneous
  multithreading, two logical CPUs per core on most x86 processors) add little to this kind of
  work.
- **Cores of unequal speed.** On processors with slower efficiency cores, a worker on an
  efficiency core finishes its job later; with many jobs `pmap` balances this and every core
  still adds throughput. On the M4 Pro above, its 4 efficiency cores together added ~10 % (about
  one performance core). With only a few long jobs, the slow workers become the tail: use only
  the fast cores.
- **RAM.** Workers × (1–1.5 GB) plus the main process must fit comfortably in memory; check with
  the first few jobs (`Sys.maxrss()` on a worker, or Activity Monitor / `top`).
- Leave a core free if you want to keep using the machine.

## 5. Threads over several sweeps (short batches)

For a batch of a few minutes, use threads over the sweeps instead (`julia -t auto`):

```julia
results = Vector{Any}(undef, length(jobs))
with_blas_threads(1) do
    Threads.@threads :dynamic for i in eachindex(jobs)
        results[i] = with_task_blas_single(() -> runJob(jobs[i]))   # runJob as above
    end
end
```

Each task writes only its own slot of `results`; never accumulate into a shared variable from
inside the loop (`pitfallsJ.md` jNUM1).

## 6. Long runs

- Test the job on two or three points first and time it; multiply up before starting the
  full batch.
- Keep the machine awake: on macOS start the run under `caffeinate -i julia ...`; on Linux run it
  inside `tmux` or with `nohup` so closing the terminal does not stop it.
- On a laptop, run on mains power: power-saving settings on battery can slow the processor down.
- Other heavy programs running at the same time slow every worker; for timings, check the
  machine is otherwise idle.
