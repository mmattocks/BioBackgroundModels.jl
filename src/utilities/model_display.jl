cs_dict=Dict('A'=>:green,'C'=>:blue, 'G'=>:yellow, 'T'=>:red)

function print_emitters(state::Categorical) #print informative symbols from the state's emission distribution, if any
    order=Integer(log(4,length(state.p))-1)
    alphabet=CompoundAlphabet(ACGT,order)
    bits=log(2,length(state.p)) #higher order hmms express more information per symbol
    dummy_p=state.p.+10^-99 #prevent -Inf or NaN from zero probability symbols
    infoscore=(bits+sum(x*log(2,x) for x in dummy_p))
    infovec=[x*infoscore for x in dummy_p]
    infovec./=bits
    print("<<< ")
    for (symbol,symbol_prob) in enumerate(infovec)
        if symbol_prob >= .05
            str=string(alphabet.integers[symbol])
            if symbol_prob >= .7
                for char in str
                    printstyled(stdout, char; bold=true, color=cs_dict[char])
                end
            elseif .7 > symbol_prob >= .25
                for char in str
                    printstyled(stdout, char; color=cs_dict[char])
                end
            elseif .25 > symbol_prob
                for char in str
                    printstyled(stdout, lowercase(char); color=cs_dict[char])
                end
            end
            print(" ")
        end
    end
    println(">>>")
end