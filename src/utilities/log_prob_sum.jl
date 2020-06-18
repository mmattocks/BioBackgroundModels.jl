                #subfuncs to handle sums of log probabilities that may include -Inf (ie p=0), returning -Inf in this case rather than NaNs
                function lps(adjuvants)
                    prob = sum(adjuvants) ; isnan(prob) ? - Inf : prob
                end
                
                function lps(base, adjuvants...)
                    prob = base+sum(adjuvants) ; isnan(prob) ? -Inf : prob
                end

