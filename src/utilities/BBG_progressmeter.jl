#UTILITY PROGRESSMETER, REPORTS WORKER NUMBER AND CURRENT ITERATE
mutable struct ProgressHMM{T<:Real} <: ProgressMeter.AbstractProgress
    thresh::T
    dt::AbstractFloat
    val::T
    counter::Integer
    triggered::Bool
    tfirst::AbstractFloat
    tlast::AbstractFloat
    printed::Bool        # true if we have issued at least one status update
    desc::AbstractString # prefix to the percentage, e.g.  "Computing..."
    color::Symbol        # default to green
    output::IO           # output stream into which the progress is written
    numprintedvalues::Integer   # num values printed below progress in last iteration
    offset::Integer             # position offset of progress bar (default is 0)
    steptime::AbstractFloat
    function ProgressHMM{T}(thresh;
                               dt::Real=0.1,
                               desc::AbstractString="Progress: ",
                               color::Symbol=:green,
                               output::IO=stderr,
                               offset::Integer=0,
                               start_it::Integer=1) where T
        tfirst = tlast = time()
        printed = false
        new{T}(thresh, dt, typemax(T), start_it, false, tfirst, tlast, printed, desc, color, output, 0, offset, 0.0)
    end
end

ProgressHMM(thresh::Real, dt::Real=0.1, desc::AbstractString="Progress: ",
         color::Symbol=:green, output::IO=stderr;
         offset::Integer=0, start_it::Integer=1) = ProgressHMM{typeof(thresh)}(thresh, dt=dt, desc=desc, color=color, output=output, offset=offset, start_it=start_it)

ProgressHMM(thresh::Real, desc::AbstractString, offset::Integer=0, start_it::Integer=1) = ProgressHMM{typeof(thresh)}(thresh, desc=desc, offset=offset, start_it=start_it)

function update!(p::ProgressHMM, val, steptime; options...)
    p.val = val
    p.counter += 1
    p.steptime = steptime
    updateProgress!(p; options...)
end

function updateProgress!(p::ProgressHMM; showvalues = Any[], valuecolor = :blue, offset::Integer = p.offset, keep = (offset == 0))
    p.offset = offset
    t = time()
    if p.val <= p.thresh && !p.triggered
        p.triggered = true
        if p.printed
            p.triggered = true
            dur = ProgressMeter.durationstring(t-p.tfirst)
            msg = @sprintf "%s Time: %s (%d iterations)" p.desc dur p.counter
            print(p.output, "\n" ^ (p.offset + p.numprintedvalues))
            ProgressMeter.move_cursor_up_while_clearing_lines(p.output, p.numprintedvalues)
            ProgressMeter.printover(p.output, msg, p.color)
            ProgressMeter.printvalues!(p, showvalues; color = valuecolor)
            if keep
                println(p.output)
            else
                print(p.output, "\r\u1b[A" ^ (p.offset + p.numprintedvalues))
            end
        end
        return
    end

    if t > p.tlast+p.dt && !p.triggered
        elapsed_time = t - p.tfirst
        msg = @sprintf "%s (convergence: %g => %g, iterate: %g, step time: %s)" p.desc p.val p.thresh p.counter hmss(p.steptime)
        print(p.output, "\n" ^ (p.offset + p.numprintedvalues))
        ProgressMeter.move_cursor_up_while_clearing_lines(p.output, p.numprintedvalues)
        ProgressMeter.printover(p.output, msg, p.color)
        ProgressMeter.printvalues!(p, showvalues; color = valuecolor)
        print(p.output, "\r\u1b[A" ^ (p.offset + p.numprintedvalues))
        # Compensate for any overhead of printing. This can be
        # especially important if you're running over a slow network
        # connection.
        p.tlast = t + 2*(time()-t)
        p.printed = true
    end
end

            function hmss(dt)
                isnan(dt) && return "NaN"
                (h,r) = divrem(dt,60*60)
                (m,r) = divrem(r, 60)
                (isnan(h)||isnan(m)||isnan(r)) && return "NaN"
                string(Int(h),":",Int(m),":",Int(ceil(r)))
            end