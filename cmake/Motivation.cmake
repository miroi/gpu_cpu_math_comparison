

function(get_motivational_message in_compiled_by_dev out_msg)
    set(msg "good job!")
    if(${in_compiled_by_dev})
        # get a random integer in the range 0-9
        string(RANDOM LENGTH 1 ALPHABET 0123456789 random_integer)
        math(EXPR number "${random_integer} + 0") # remove extra leading 0s
        
        set(msg_list
            "good job!"
            "well done!"
            "another job well done!"
            "perhaps time for a coffee?"
            "nice!"
            "excellent!"
            "you are a beautiful person writing a beautiful code." # thanks to Yesod
            "keep going, you have nearly fixed the bug!"           # thanks to Yesod
            "check your posture, don’t lean over too much."        # thanks to Yesod
            "get a glass of water, keep hydrated."                 # thanks to Yesod
            )
        
        list(GET msg_list ${random_integer} msg)
    endif()
    set(${out_msg} ${msg} PARENT_SCOPE)
endfunction()
